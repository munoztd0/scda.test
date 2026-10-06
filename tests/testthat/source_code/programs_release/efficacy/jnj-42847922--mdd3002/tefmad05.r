###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefmad05.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefmad05:
##                            Change From Baseline[(DB)] Over Time: MMRM [Observed Case] Analysis
##                            [ - Double-blind Phase]; Full Analysis Set Analysis Set (Study mdd)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADMADRSI
## Output:                    TEFMAD05.rtf
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

# Define output ID and file location.
tblid <- "TEFMAD05"
fileid <- write_path(opath, tblid)


# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Placebo"

# Define population flags used.
popfl <- "FAS1FL"

# Define the covariates to use in the MMRM.
# Note: The main effects and interaction of AVISIT and trtvar are always
# included already.
covariates <- c("BASE", "RRSIWRS", "COUNTRY", "AGEGR2")

# Define analysis domain flags to be used (for ADMADRSI).
domfl <- quote(ANL03FL == "Y" & DTYPE != "ENDPOINT")

# Define visit flags to be used.
visitfl <- quote(!(AVISIT %in% c("Baseline (DB)", "Endpoint (DBFU)")))

# Define response parameter to be used (for ADMADRSI).
resppar <- "MADES1"

# Define between one and two sided p-values.
pval_sided <- "2"
# one of:
# "2" - two sided
# "1" - one sided, greater
# "-1" - one sided, less

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
  mse = jjcsformat_xx("xx.xx"),
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
  filter(!!visitfl) |>
  filter(PARAMCD == resppar) |>
  select(USUBJID, AVISIT, AVAL, BASE, CHG) |>
  mutate(
    USUBJID = factor(USUBJID),
    AVISIT = factor(AVISIT)
  )

## Analysis ----

ana <- admadrsi |>
  mutate(AVISIT = factor(paste0(AVISIT, dblabel))) |>
  inner_join(adsl, by = "USUBJID")

# MMRM analysis ----

variables <- list(
  response = "CHG",
  covariates = covariates,
  id = "USUBJID",
  arm = trtvar,
  visit = "AVISIT"
)
model_obj <- fit_mmrm_j(
  vars = variables,
  data = ana,
  cor_struct = "unstructured",
  conf_level = conflvl,
  method = "Kenward-Roger",
  vcov = "Kenward-Roger-Linear",
  weights_emmeans = "equal"
)
ana_results <- broom::tidy(model_obj)

# Layout ----

lyt <- basic_table() |>
  summarize_lsmeans_wide(
    variables = variables,
    pval_sided = pval_sided,
    conf_level = conflvl,
    ref_level = model_obj$ref_level,
    treatment_levels = model_obj$treatment_levels,
    formats = formats
  )

# Output ----

result <- build_table(lyt, ana_results)


# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
