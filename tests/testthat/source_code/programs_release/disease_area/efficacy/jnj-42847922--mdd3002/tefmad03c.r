###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefmad03c.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefmad03c: Change From Baseline[(DB)] to [Time Point]:
##                            Copy Increment From Reference Multiple Imputation Analysis[ - Double-blind Phase];
##                            Full Analysis Set Analysis Set (Study mdd)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADMADRSI
## Output:                    TEFMAD03c.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:                     1
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
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
tblid <- "TEFMAD03c"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Placebo"

# Define population flags used.
popfl <- "FAS2FL"

# Define the covariates to use in the MMRM (imputation and analysis).
# Note: AVISIT and trtvar are always included already.
covariates <- c("BASE", "RRSIWRS", "COUNTRY", "AGEGR2")

# Interaction terms based on covariates, AVISIT and trtvar
# to be used in the MMRM (imputation and analysis).
interactions <- c(paste(trtvar, "* AVISIT"))

# Define analysis domain flags to be used (for ADMADRSI).
domfl <- quote(ANL03FL == "Y" & DTYPE != "ENDPOINT")

# Define response parameter to be used (for ADMADRSI).
resppar <- "MADES2"

# Define time point to be used.
timept <- "Day 43"

# Define baseline time point label.
basetpt <- "Baseline (DB)"

# Define between one and two sided p-values.
pval_sided <- "2"
# one of:
# "2" - two sided
# "-1" - one sided, less
# "1" - one sided, greater

# Settings for the imputation method.
library(rbmi)
set.seed(3469) # Needed because method_bayes generates a random seed.
imputation_method <- rbmi::method_bayes(
  control = rbmi::control_bayes(
    warmup = 200,
    thin = ifelse(testrun, 1, 5)
  ),
  n_samples = ifelse(testrun, 10, 500) # Number of imputed datasets to be used.
)

# Reference based imputation method to be used.
reference_method <- "CIR"

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
  mean_sd = jjcsformat_xx("xx.x (xx.xx)"),
  median = jjcsformat_xx("xx.x"),
  range = jjcsformat_xx("xx, xx"),
  pval = jjcsformat_pval_fct(alpha),
  diff_mean_se = jjcsformat_xx("xx.x (xx.xx)"),
  diff_mean_ci = jjcsformat_xx("(xx.xx, xx.xx)")
)
dblabel <- if (doubleblind) " (DB)" else ""

set_default_na_str("NE")

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  mutate(!!trtvar := forcats::fct_relevel(!!rlang::sym(trtvar), ctrlab)) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl), any_of(covariates))

adsl$colspan_trt <- factor(
  ifelse(adsl[[trtvar]] == ctrlab, " ", "Active Study Agent"),
  levels = c("Active Study Agent", " ")
)

colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = ctrlab,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

ref_path <- c("colspan_trt", " ", trtvar, ctrlab)

## ADMADRSI ----

admadrsi <- haven::read_sas(read_path(a_in, "admadrsi.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(PARAMCD == resppar) |>
  filter(!!domfl) |>
  select(USUBJID, AVISIT, AVAL, BASE, CHG) |>
  mutate(
    USUBJID = factor(USUBJID),
    AVISIT = factor(AVISIT)
  )

## Analysis ----

# For descriptive statistics.
ana_desc <- admadrsi |>
  filter(AVISIT %in% c(basetpt, timept)) |>
  mutate(
    AVISIT = stringr::str_replace(AVISIT, timept, paste0(timept, dblabel))
  ) |>
  inner_join(adsl, by = "USUBJID")

# For multiple imputation.
ana_mi <- admadrsi |>
  inner_join(adsl, by = "USUBJID") |>
  # Complete patient/visits grid.
  expand_locf(
    USUBJID = levels(admadrsi$USUBJID),
    AVISIT = levels(admadrsi$AVISIT),
    vars = c("STUDYID", popfl, trtvar, covariates),
    group = c("USUBJID"),
    order = c("USUBJID", "AVISIT")
  ) |>
  # Drop baseline visit (because there is no change to model yet).
  filter(AVISIT != basetpt) |>
  mutate(AVISIT = droplevels(AVISIT))

# Multiple Imputation ----

## Intercurrent events ----

# This could look more complicated, when different strategies are used
# for different events. Here we follow the approach that
# - intermittent missing values are imputed using hypothetical strategy (MAR)
# - final missing values are imputed using the copy-reference (CR) strategy
#   and therefore in `ana_ice` the first such missing value is defined.
# All complete patients, or patients who only have intermittent missing values
# are therefore not included in `ana_ice`.
ana_ice <- ana_mi |>
  select(USUBJID, AVISIT, CHG) |>
  arrange(USUBJID, AVISIT) |>
  group_by(USUBJID) |>
  summarize(AVISIT = find_missing_chg_after_avisit(pick(AVISIT, CHG))) |>
  ungroup() |>
  filter(!is.na(AVISIT)) |>
  mutate(AVISIT = factor(AVISIT, levels = levels(ana_mi$AVISIT))) |>
  mutate(strategy = reference_method)

## Bayesian imputation ----

vars <- rbmi::set_vars(
  outcome = "CHG",
  visit = "AVISIT",
  subjid = "USUBJID",
  group = trtvar,
  covariates = c(covariates, interactions)
)

set.seed(12345)
wrapper <- if (testrun) suppressWarnings else identity
draw_obj <- wrapper(rbmi::draws(
  data = ana_mi,
  data_ice = ana_ice, # Here we pass the intercurrent events information already.
  vars = vars,
  method = imputation_method,
  quiet = TRUE
))
impute_obj <- rbmi::impute(
  draw_obj,
  references = reference_spec
)

## Analysis step ----

ana_obj <- rbmi::analyse(
  impute_obj,
  rbmi_mmrm,
  vars = set_vars(
    subjid = "USUBJID",
    outcome = "CHG",
    visit = "AVISIT",
    group = trtvar,
    covariates = c(covariates, interactions)
  ),
  cov_struct = "us",
  weights = "equal",
  reml = TRUE,
  method = "Kenward-Roger",
  vcov = "Kenward-Roger-Linear"
)

## Pooling ----

pool_obj <- junco:::rbmi_pool(
  ana_obj,
  conf.level = conflvl,
  alternative = switch(
    pval_sided,
    "2" = "two.sided",
    "-1" = "less",
    "1" = "greater"
  )
)
ana_mi_results <- broom::tidy(pool_obj, visits = levels(ana_mi$AVISIT)) |>
  filter(visit == timept) |>
  rename(!!trtvar := group) |>
  mutate(
    colspan_trt = factor(
      ifelse(!!rlang::sym(trtvar) == ctrlab, " ", "Active Study Agent"),
      levels = c("Active Study Agent", " ")
    )
  )

# Layout ----

## Descriptive analysis ----

lyt_desc <- rtables::basic_table(
  top_level_section_div = " ",
  show_colcounts = TRUE,
  colcount_format = "N=xx"
) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  split_rows_by("AVISIT", section_div = "") |>
  analyze_values("AVAL", formats = formats) |>
  analyze_values(
    "CHG",
    var_labels = paste0("Change from baseline", dblabel),
    formats = formats,
    show_labels = "visible",
    nested = FALSE
  )

## Multiple imputation results ----

lyt_mi <- rtables::basic_table(
  top_level_section_div = " ",
  show_colcounts = TRUE,
  colcount_format = "N=xx"
) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  analyze(
    trtvar,
    afun = a_rbmi_lsmeans,
    show_labels = "hidden",
    na_str = default_na_str(),
    extra_args = list(
      .stats = c(
        "p_value",
        "diff_mean_se",
        "diff_mean_ci"
      ),
      .labels = c(
        p_value = paste0(
          pval_sided,
          "-sided p-value (minus ",
          ctrlab,
          ")~[super a,b]"
        ),
        diff_mean_se = "Diff. of LS Means (SE)",
        diff_mean_ci = tern::f_conf_level(conflvl)
      ),
      .formats = c(
        p_value = formats$pval,
        diff_mean_se = formats$diff_mean_se,
        diff_mean_ci = formats$diff_mean_ci
      ),
      .indent_mods = c(diff_mean_se = 1L),
      ref_path = ref_path
    )
  )

# Output ----

result_desc <- build_table(lyt_desc, df = ana_desc, alt_counts_df = adsl)
result_mi <- build_table(lyt_mi, ana_mi_results, alt_counts_df = adsl)
result <- rbind(result_desc, result_mi)


# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
