###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefos03.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefos03: Overall Survival (Multivariate Analysis)
##                            Proportional Hazards Model; Full Analysis Set (Study mmy, bc, lc)
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADTTEEF
## Output:                    TEFOS03.rtf
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
tblid <- "TEFOS03"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
# Warning: Title file should contain exactly one title record per Table ID
# Therefore need to make sure we only have one title record:
tab_titles$title <- tab_titles$title[1]
tab_titles$main_footer <- tab_titles$main_footer[1]


# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Dummy B - DPd"

# Define population flag used (here take measurable disease at baseline flag).
popfl <- "ITTFL"

# Define time-to-event parameter to be used.
ttepar <- "OS"

# Define additional covariates (incl. stratification variables)
# to be used in the multivariate model.
covariates <- c("AGE", "SEX", "RACEGR1", "STRAT02A", "STRAT03")

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived formats specifications.
formats <- list(
  coef_se = jjcsformat_xx("xx.xx (xx.xx)"),
  hr_est = jjcsformat_xx("xx.xx"),
  hr_ci = jjcsformat_xx("(xx.xx, xx.xx)"),
  pval = jjcsformat_pval_fct(alpha)
)

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(!!trtvar := forcats::fct_relevel(.data[[trtvar]], ctrlab)) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl), all_of(covariates))

## ADTTEEF ----

adtteef <- haven::read_sas(read_path(a_in, "adtteef.sas7bdat")) |>
  select(STUDYID, USUBJID, PARAMCD, CNSR, AVAL, all_of(trtvar))

adtteef <- adtteef |>
  select(USUBJID, PARAMCD, CNSR, AVAL) |>
  filter(PARAMCD == ttepar) |>
  select(-PARAMCD) |>
  mutate(
    USUBJID = factor(USUBJID),
    EVENT = 1 - CNSR,
    is_event = (EVENT == "Event")
  )

## Analysis ----

ana <- adtteef |>
  right_join(adsl, by = "USUBJID")

# Layout ----

variables <- list(
  time = "AVAL",
  event = "EVENT",
  arm = trtvar,
  covariates = covariates
)
control <- tern::control_coxreg(
  ties = "breslow",
  conf_level = conflvl
)
lyt <- basic_table() |>
  summarize_coxreg_multivar(
    var = "USUBJID",
    variables = variables,
    control = control,
    formats = formats
  )

# Output ----

result <- build_table(lyt, ana)

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
