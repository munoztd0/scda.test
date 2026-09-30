###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefmadsg02.r
## R Version:                 4.4.2
## junco Version:             0.1.6 (tbc)
## Short Description:         Program to create tefmadsg02:
##                            [Montgomery-Asberg Depression Rating Scale (MADRS) Total Score]
##                            Change From Baseline[(DB)] Over Time:
##                            ANCOVA by [Subgroup][- Double-blind Phase];
##                            Full Analysis Set (Study mdd)
## Author:                    Technology Solutions
## Date:                      2026-09-302026
## Input:                     ADSL, ADMADRSI
## Output:                    TEFMADSG02.rtf
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
tblid <- "TEFMADSG02"
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

# Define the covariates to use in the ANCOVA.
# Note: trtvar is always included already.
covariates <- c("BASE", "RRSIWRS", "COUNTRY")

# Define analysis domain flags to be used (for ADMADRS).
domfl <- quote(
  ANL03FL == "Y" &
    DTYPE != "ENDPOINT" &
    !(AVISIT %in% c("Baseline (DB)", "Endpoint (DBFU)"))
)

# Define response parameter to be used (for ADMADRS).
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

# ANCOVA analysis ----

variables <- list(
  response = "CHG",
  covariates = covariates,
  id = "USUBJID",
  arm = trtvar,
  visit = "AVISIT"
)

ana_by_subgr <- split(ana, ana[[subgroup]], drop = TRUE)
model_by_subgr <- lapply(
  ana_by_subgr,
  FUN = fit_ancova,
  vars = variables,
  conf_level = conflvl,
  weights_emmeans = "equal"
)
results_by_subgr <- lapply(model_by_subgr, broom::tidy)
ana_results <- bind_rows(results_by_subgr, .id = subgroup)

# Layout ----

lyt <- basic_table() |>
  split_cols_by(
    variables$arm,
    split_fun = lsmeans_wide_first_split_fun_fct(include_variance = FALSE)
  ) |>
  split_cols_by(
    variables$arm,
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
  split_rows_by(variables$visit, section_div = "") |>
  summarize_row_groups(
    var = variables$arm,
    cfun = lsmeans_wide_cfun,
    extra_args = list(
      variables = variables,
      ref_level = model_by_subgr[[1]]$ref_level,
      treatment_levels = model_by_subgr[[1]]$treatment_levels,
      pval_sided = pval_sided,
      conf_level = conflvl,
      formats = formats
    )
  )


# Output ----

result <- build_table(lyt, ana_results)

# Add blank lines after each subgroup page split label.
section_div_at_path(result, c(subgroup, "*", "@content", "*")) <- " "

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
