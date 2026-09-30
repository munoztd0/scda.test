###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              teforr01a.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create teforr01a:  Objective Response Rate Based on
##                            RECIST [Version 1.1] Criteria in Subjects With Measurable Disease at Baseline -
##                            [Stratified/Unstratified] Analysis; Full Analysis Set (Study mmy, bc, lc)
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADEFF
## Output:                    TEFORR01a.rtf
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
tblid <- "TEFORR01a"
fileid <- write_path(opath, tblid)

# Please select titles and footnotes as appropriate.
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
tab_titles$title <- tab_titles$title[1]
tab_titles$main_footer <- tab_titles$main_footer[c(1, 4)]

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "CP"

# Define treatment group label to look at here (only one treatment group allowed).
trtlab <- "ACP"

# Define population flag used (here take measurable disease at baseline flag).
popfl <- "MDBIRCFL"

# Define the strata variables to use below.
strata <- c("STRATA1", "STRATA2", "STRATA3")
# For unstratified analyses, please use:
# strata <- NULL

# Define response parameter to be used.
resppar <- "BORCIRC"

# Define method for response rate confidence intervals.
# See ?tern::s_proportion for possible options.
cimethod <- "clopper-pearson"

# Define method for comparing response.
method <- "or_cmh"
# one of:
# "rr": relative risk
# "or_logistic": odds ratio by logistic regression
# "or_cmh": odds ratio by Cochran-Mantel-Haenszel (CMH)

# If unstratified, then this choice of the test matters:
testmethod <- "fisher"
# one of:
# "fisher": Fisher's exact test
# "chisq": Chi-Square test

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived formats and methods specifications.
formats <- list(
  comp_stat_ci = jjcsformat_xx("xx.xx (xx.xx, xx.xx)"),
  pval = jjcsformat_pval_fct(alpha),
  prop_ci = jjcsformat_xx("(xx.x%, xx.x%)")
)
set_default_na_str("NE")
methods <- list(
  comp_stat_ci = method,
  pval = testmethod,
  prop_ci = cimethod
)

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(
    !!rlang::sym(popfl) == "Y",
    !!rlang::sym(trtvar) %in% c(ctrlab, trtlab)
  ) |>
  mutate(!!trtvar := factor(.data[[trtvar]], levels = c(trtlab, ctrlab))) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl), all_of(strata))

## ADEFF ----

adeff <- haven::read_sas(read_path(a_in, "adeff.sas7bdat")) |>
  select(STUDYID, USUBJID, PARAMCD, AVALC, all_of(trtvar))

adeff <- adeff |>
  select(USUBJID, PARAMCD, AVALC) |>
  filter(PARAMCD == resppar) |>
  select(-PARAMCD) |>
  mutate(
    USUBJID = factor(USUBJID),
    AVALC = forcats::fct_recode(
      AVALC,
      "CR" = "Complete Response (CR)",
      "PR" = "Partial Response (PR)",
      "SD" = "Stable Disease (SD)",
      "PD" = "Progressive Disease (PD)",
      "NE" = "Not Evaluable/Unknown",
      "NE" = "Non-CR/Non-PD"
    )
  )

## Analysis ----

ana <- adeff |>
  mutate(
    rsp_lab = jjcs_lung_rsp_label(AVALC),
    is_rsp_any = AVALC %in% c("CR", "PR"),
    is_dcr = AVALC %in% c("CR", "PR", "SD")
  ) |>
  right_join(adsl, by = "USUBJID")

# Functions ----

analyze_response_logical <- function(lyt, vars, label, include_comp) {
  analyze(
    lyt = lyt,
    vars = vars,
    afun = resp01_acfun,
    show_labels = "hidden",
    extra_args = list(
      arm = trtvar,
      include_comp = include_comp,
      conf_level = conflvl,
      strata = strata,
      label = label,
      methods = methods,
      formats = formats
    )
  )
}

# Layout ----

lyt <- basic_table() |>
  split_cols_by(
    trtvar,
    show_colcounts = TRUE,
    # Note: It needs to stay "Overall" here (this won't be shown anyway).
    split_fun = add_overall_level("Overall", label = "", first = FALSE)
  ) |>
  split_cols_by(
    "STUDYID", # The particular choice of this is irrelevant.
    split_fun = resp01_split_fun_fct(method = method, conf_level = conflvl)
  ) |>
  ## Response category ----
  analyze(
    vars = "rsp_lab",
    afun = resp01_acfun,
    show_labels = "visible",
    var_labels = "Response category",
    extra_args = list(
      arm = trtvar,
      include_comp = "Complete response (CR)",
      conf_level = conflvl,
      strata = strata,
      formats = formats,
      methods = methods
    )
  ) |>
  insert_blank_line() |>
  ## Overall response ----
  analyze_response_logical(
    vars = "is_rsp_any",
    label = "Overall response (CR+PR)",
    include_comp = TRUE
  ) |>
  ## Disease control rate ----
  analyze_response_logical(
    vars = "is_dcr",
    label = "Disease control rate (CR+PR+SD)",
    include_comp = FALSE
  ) |>
  append_topleft("Response")

# Output ----

result <- build_table(lyt, ana, alt_counts_df = adsl)

# Post-hoc we can suppress the column count for the Overall column.
colcount_visible(result, c(trtvar, "Overall")) <- FALSE


# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
