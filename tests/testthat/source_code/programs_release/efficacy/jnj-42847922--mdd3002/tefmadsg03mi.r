###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefmadsg03mi.r
## R Version:                 4.5.2
## junco Version:             0.1.6
## Short Description:         Program to create tefmadsg03:
##                            [Montgomery-Asberg Depression Rating Scale (MADRS)
##                            Total Score]: Change From Baseline[(DB)] Over Time:
##                            MMRM Observed Case by [Subgroup][- Double-blind Phase]
## Note:                      This version is with multiple imputation. Please see tefmadsg03.r for
##                            the version without multiple imputation.
## Author:                    C&SP Methodology
## Date:                      14 May 2026
## Input:                     adsl, admadrsi
## Output:                    tefmadsg03mi.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:
## Modified By:
## Reporting effort:
## Date:
## Description:
################################################################################

# Environment ----

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)

library(dplyr)
library(rtables)
library(junco)
library(haven)

# Parameters ----

# Test run or real production run? Please set to FALSE for use in the study analysis.
testrun <- TRUE

# Define output ID and file location.
tblid <- "TEFMADSG03mi"
fileid <- write_path(opath, tblid)

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Placebo"

# Define population flags used.
popfl <- "FAS1FL"

# Define subgroup variable used, alongside the label to be used and the order of the
# categories.
subgroup <- "AGEGR1"
subgrlbl <- "Age Group"
subgrorder <- c("18-34 years", "35-54 years", "55-64 years", ">=65 years")

# Define the covariates to use in the MMRM (imputation and analysis).
# Note: AVISIT and trtvar are always included already.
covariates <- c("BASE", "RRSIWRS", "COUNTRY")

# Interaction terms based on covariates, AVISIT and trtvar
# to be used in the MMRM (imputation and analysis).
interactions <- c(paste(trtvar, "* AVISIT"))

# Define analysis domain flags to be used (for ADMADRS).
domfl <- quote(ANL03FL == "Y" & DTYPE != "ENDPOINT")

# Define visit flags to be used.
visitfl <- quote(!(AVISIT %in% c("Baseline (DB)", "Endpoint (DBFU)")))

# Define response parameter to be used (for ADMADRS).
resppar <- "MADES1"

# Define between one and two sided p-values.
pval_sided <- "2"
# one of:
# "2" - two sided
# "1" - one sided, greater
# "-1" - one sided, less

# Settings for the imputation method.
library(rbmi)
set.seed(3469) # Needed because method_bayes generates a random seed.
imputation_method <- rbmi::method_bayes(
  control = rbmi::control_bayes(
    warmup = ifelse(testrun, 10, 200),
    thin = ifelse(testrun, 1, 5)
  ),
  n_samples = ifelse(testrun, 10, 500) # Number of imputed datasets to be used.
)

# Reference based imputation method to be used.
reference_method <- "CR"

# References to be used during imputation.
reference_spec <- c(
  "Placebo" = "Placebo",
  "Seltorexant 20 mg" = "Placebo"
)

# Double blind labeling in the output or not.
doubleblind <- TRUE

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived formats specifications.
formats <- list(
  lsmean = jjcsformat_xx("xx.x"),
  mse = jjcsformat_xx("xx.x"),
  df = jjcsformat_xx("xx."),
  lsmean_diff = jjcsformat_xx("xx.x"),
  se = jjcsformat_xx("xx.xx"),
  ci = jjcsformat_xx("(xx.xx, xx.xx)"),
  pval = jjcsformat_pval_fct(alpha)
)
dblabel <- if (doubleblind) " (DB)" else ""

set_default_na_str("NE")

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  mutate(!!trtvar := forcats::fct_relevel(!!rlang::sym(trtvar), ctrlab)) |>
  mutate(!!subgroup := factor(.data[[subgroup]], levels = subgrorder)) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(subgroup), any_of(covariates))

## ADMADRSI ----

admadrsi <- haven::read_sas(read_path(a_in, "admadrsi.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(!!domfl) |>
  filter(!!visitfl) |>
  filter(PARAMCD == resppar) |>
  select(USUBJID, AVISIT, AVAL, BASE, CHG) |>
  mutate(
    USUBJID = factor(USUBJID),
    AVISIT = factor(paste0(AVISIT, dblabel))
  )

# Move BASE to adsl because it is a time-constant, baseline covariate.
base_vals <- admadrsi |>
  select(USUBJID, BASE) |>
  distinct()
adsl <- adsl |>
  inner_join(base_vals, by = "USUBJID")
admadrsi <- admadrsi |>
  select(-BASE)

## Analysis ----

# For multiple imputation.
ana_mi <- admadrsi |>
  # Complete patient/visits grid.
  expand(
    USUBJID = levels(admadrsi$USUBJID),
    AVISIT = levels(admadrsi$AVISIT)
  ) |>
  inner_join(adsl, by = "USUBJID") |>
  arrange(USUBJID, AVISIT)

# Split by subgroup.
# Note that:
# - In contrast to tefmadsg03 without multiple imputation there is no option here
#   to have joint MMRM fitted across subgroups. This is a deliberate choice
#   based on the statistical methodology.
# - We only here make the subject variable a factor to avoid
# levels with no observations in the subgroups.
ana_mi_by_subgr <- split(ana_mi, ana_mi[[subgroup]], drop = TRUE) |>
  lapply(function(x) mutate(x, USUBJID = factor(USUBJID)))

# Multiple Imputation ----

## Intercurrent events ----

# This could look more complicated, when different strategies are used
# for different events. Here we follow the approach that
# - intermittent missing values are imputed using hypothetical strategy (MAR)
# - final missing values are imputed using the copy-reference (CR) strategy
#   and therefore in `ana_ice` the first such missing value is defined.
# All complete patients, or patients who only have intermittent missing values
# are therefore not included in `ana_ice`.
ana_ice_by_subgr <- lapply(
  ana_mi_by_subgr,
  function(x) {
    x |>
      select(USUBJID, AVISIT, CHG) |>
      arrange(USUBJID, AVISIT) |>
      group_by(USUBJID) |>
      summarize(AVISIT = find_missing_chg_after_avisit(pick(AVISIT, CHG))) |>
      ungroup() |>
      filter(!is.na(AVISIT)) |>
      mutate(AVISIT = factor(AVISIT, levels = levels(ana_mi$AVISIT))) |>
      mutate(strategy = reference_method)
  }
)

## Bayesian imputation ----

vars_draw <- rbmi::set_vars(
  outcome = "CHG",
  visit = "AVISIT",
  subjid = "USUBJID",
  group = trtvar,
  covariates = c(covariates, interactions)
)
vars_draw$arm <- trtvar # Needed for table building at the end.

set.seed(12345)
wrapper <- if (testrun) suppressWarnings else identity
draw_obj_by_subgr <- mapply(
  function(x, y) {
    wrapper(rbmi::draws(
      data = x,
      data_ice = y, # Here we pass the intercurrent events information already.
      vars = vars_draw,
      method = imputation_method,
      quiet = TRUE
    ))
  },
  x = ana_mi_by_subgr,
  y = ana_ice_by_subgr,
  SIMPLIFY = FALSE
)
impute_obj_by_subgr <- lapply(
  draw_obj_by_subgr,
  rbmi::impute,
  references = reference_spec
)

## Analysis step ----

vars_analyse <- rbmi::set_vars(
  subjid = "USUBJID",
  outcome = "CHG",
  visit = "AVISIT",
  group = trtvar,
  covariates = c(covariates, interactions)
)

ana_obj_by_subgr <- lapply(
  impute_obj_by_subgr,
  rbmi_analyse,
  fun = rbmi_mmrm,
  vars = vars_analyse,
  cov_struct = "us",
  weights = "equal",
  reml = TRUE,
  method = "Kenward-Roger",
  vcov = "Kenward-Roger-Linear"
)

## Pooling ----

pool_obj_by_subgr <- lapply(
  ana_obj_by_subgr,
  junco:::rbmi_pool,
  conf.level = conflvl,
  alternative = switch(
    pval_sided,
    "2" = "two.sided",
    "-1" = "less",
    "1" = "greater"
  )
)

## Results data ----

n_data <- ana_mi |>
  group_by(AVISIT, !!rlang::sym(trtvar), !!rlang::sym(subgroup)) |>
  summarize(n = sum(!is.na(CHG)))

ana_mi_results <- lapply(pool_obj_by_subgr, broom::tidy, visits = levels(ana_mi$AVISIT)) |>
  bind_rows(.id = subgroup) |>
  mutate(p_value_less = p_value, p_value_greater = p_value) |>
  rename(
    !!trtvar := group,
    AVISIT = visit,
    estimate_est = est,
    estimate_contr = est_contr
  ) |>
  inner_join(n_data) |>
  relocate(!!rlang::sym(subgroup), AVISIT, !!rlang::sym(trtvar))

# Layout ----

lyt <- basic_table() |>
  split_cols_by(
    vars_draw$group,
    split_fun = lsmeans_wide_first_split_fun_fct(include_variance = FALSE)
  ) |>
  split_cols_by(
    vars_draw$group,
    split_fun = lsmeans_wide_second_split_fun_fct(
      include_pval = TRUE,
      pval_sided = pval_sided,
      conf_level = conflvl
    )
  ) |>
  split_rows_by(
    subgroup,
    page_by = TRUE,
    split_fun = drop_split_levels,
    split_label = subgrlbl,
    section_div = " "
  ) |>
  # Note: This line is currently needed to keep the above page split labels in the rtf,
  # as a workaround for an issue in rtables.
  summarize_row_groups(cfun = ac_blank_line) |>
  split_rows_by(vars_draw$visit, section_div = "") |>
  summarize_row_groups(
    var = vars_draw$group,
    cfun = lsmeans_wide_cfun,
    extra_args = list(
      variables = vars_draw,
      ref_level = levels(ana_mi[[trtvar]])[1],
      treatment_levels = levels(ana_mi[[trtvar]])[-1],
      pval_sided = pval_sided,
      conf_level = conflvl,
      formats = formats
    )
  )

# Output ----

result <- build_table(lyt, ana_mi_results, round_type = "sas")

# Add blank lines after each subgroup page split label.
section_div_at_path(result, c(subgroup, "*", "@content", "*")) <- " "

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
