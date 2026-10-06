###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefmad04_mi.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefmad04:
##                            Change From Baseline[(DB)] Over Time: ANCOVA [Imputation
##                            Method][ - Double-blind Phase]; Full Analysis Set Analysis
##                            Set (Study mdd).
## Note:                      This version is with multiple imputation. Please see tefmad04.r for
##                            the version without multiple imputation.
## Author:                    C&SP Methodology
## Date:                      23 Jan 2025
## Input:                     ADSL, ADMADRSI
## Output:                    TEFMAD04mi.rtf
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
tblid <- "TEFMAD04mi"
fileid <- write_path(opath, tblid)

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Placebo"

# Define population flags used.
popfl <- "FAS1FL"

# Define the covariates to use in the MMRM (imputation) and ANCOVA (analysis).
# Note: AVISIT (MMRM) and trtvar (MMRM and ANCOVA) are always included already.
covariates <- c("BASE", "RRSIWRS", "COUNTRY", "AGEGR2")

# Interaction terms based on covariates, AVISIT and trtvar
# to be used in the MMRM (imputation).
interactions <- c("BASE * AVISIT", paste(trtvar, "* AVISIT"))

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
    warmup = 200,
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
  select(STUDYID, USUBJID, all_of(trtvar), any_of(covariates))

## ADMADRSI ----

admadrsi <- haven::read_sas(read_path(a_in, "admadrsi.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(!!domfl) |>
  filter(PARAMCD == resppar) |>
  select(USUBJID, AVISIT, AVAL, BASE, CHG) |>
  mutate(
    USUBJID = factor(USUBJID),
    AVISIT = factor(AVISIT)
  )

## Analysis ----

# For multiple imputation.
ana_mi <- admadrsi |>
  inner_join(adsl, by = "USUBJID") |>
  # Complete patient/visits grid.
  expand_locf(
    USUBJID = levels(admadrsi$USUBJID),
    AVISIT = levels(admadrsi$AVISIT),
    vars = c("STUDYID", trtvar, covariates),
    group = c("USUBJID"),
    order = c("USUBJID", "AVISIT")
  ) |>
  filter(!!visitfl) |>
  mutate(AVISIT = factor(paste0(AVISIT, dblabel)))

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
  rbmi_ancova,
  vars = set_vars(
    subjid = "USUBJID",
    outcome = "CHG",
    visit = "AVISIT",
    group = trtvar,
    covariates = covariates
  ),
  weights = "equal"
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

## Results data ----

n_data <- ana_mi |>
  group_by(AVISIT, !!rlang::sym(trtvar)) |>
  summarize(n = sum(!is.na(CHG)))

ana_mi_results <- broom::tidy(pool_obj, visits = levels(ana_mi$AVISIT)) |>
  mutate(p_value_less = p_value, p_value_greater = p_value) |>
  rename(
    !!trtvar := group,
    AVISIT = visit,
    estimate_est = est,
    estimate_contr = est_contr
  ) |>
  inner_join(n_data) |>
  relocate(AVISIT, !!rlang::sym(trtvar))

# Layout ----

variables <- list(
  response = "CHG",
  covariates = covariates,
  id = "USUBJID",
  arm = trtvar,
  visit = "AVISIT"
)
trtlab <- levels(ana_mi[[trtvar]])[-1]
lyt <- basic_table() |>
  summarize_lsmeans_wide(
    variables = variables,
    pval_sided = pval_sided,
    conf_level = conflvl,
    ref_level = ctrlab,
    treatment_levels = trtlab,
    formats = formats
  )

# Output ----

result <- build_table(lyt, ana_mi_results)


# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
