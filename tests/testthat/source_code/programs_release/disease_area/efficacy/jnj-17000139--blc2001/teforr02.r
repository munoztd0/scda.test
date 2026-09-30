###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

## Original Reporting Effort: Standards
## Program Name:              teforr02.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create teforr02: Overall Best [Confirmed] Response
##                            Rate Based on RECIST [Version 1.1] Criteria in Subjects With Measurable Disease
##                            at Baseline by [Subgroup] - [Stratified/Unstratified] Analysis; Full
##                            Analysis Set (Study mmy, bc, lc)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADEFF
## Output:                    TEFORR02.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
##  Rev #:                    1
##  Modified By:
##  Reporting Effort:
##  Date:
##  Description:

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
tblid <- "TEFORR02"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
tab_titles$title <- tab_titles$title[1]
tab_titles$main_footer <- tab_titles$main_footer[c(1, 2, 5)]

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define treatment groups to be included.
trtlab <- c("")

# Define control group label used in the treatment variable.
ctrlab <- "CETRELIMAB"

# Define population flag used.
popfl <- "FASFL"

# Define subgroup variable used, alongside the label to be used.
subgroup <- "AGEGR1"
subgrlbl <- "Age Group"

# Define the strata variables to use below.
strata <- c("REGION1")
# For unstratified analyses, please use:
# strata <- NULL

# Define response parameter to be used.
resppar <- "BORSPIRC"

# Define method for comparing response.
method <- "rr"
# one of:
# "rr": relative risk
# "or_logistic": odds ratio by logistic regression
# "or_cmh": odds ratio by Cochran-Mantel-Haenszel (CMH)

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived format specifications.
formats <- list(
  pval = jjcsformat_pval_fct(alpha),
  ci = jjcsformat_xx("xx.xx (xx.xx, xx.xx)"),
  est_ci = jjcsformat_xx("xx.xx (xx.xx, xx.xx)")
)

options(tern_default_na_str = rep("NE", 10))

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(
    !!trtvar := as.factor(.data[[trtvar]]),
    !!subgroup := as.factor(.data[[subgroup]])
  ) |>
  select(
    STUDYID,
    USUBJID,
    all_of(trtvar),
    all_of(popfl),
    all_of(strata),
    all_of(subgroup)
  )

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

## ADEFF ----

adeff <- haven::read_sas(read_path(a_in, "adeff.sas7bdat")) |>
  select(STUDYID, USUBJID, PARAMCD, AVALC, all_of(trtvar))

adeff <- adeff |>
  select(USUBJID, PARAMCD, AVALC) |>
  filter(PARAMCD == resppar) |>
  select(-PARAMCD) |>
  mutate(
    USUBJID = factor(USUBJID),
    AVALC = forcats::fct_recode(AVALC, "NE" = "NON-CR")
  )

## Analysis ----
ana <- adeff |>
  mutate(
    rsp_lab = jjcs_lung_rsp_label(AVALC),
    is_rsp_cr_pr = AVALC %in% c("CR", "PR"),
    is_rsp_cr = AVALC == "CR",
    is_rsp_cr_pr_sd = AVALC %in% c("CR", "PR", "SD")
  ) |>
  right_join(adsl, by = "USUBJID")

# Functions ----

# Layout generating function for the response variable analysis.
analyze_response <- function(
  lyt,
  vars,
  var_labels,
  method = c("rr", "or_logistic", "or_cmh")
) {
  method <- match.arg(method)

  ## Response proportion ----
  lyt <- lyt |>
    analyze(
      vars = vars,
      afun = s_proportion_logical,
      table_names = paste0(vars, "est_prop"),
      extra_args = list(label = var_labels),
      show_labels = "hidden"
    )

  lyt <- if (method == "rr") {
    ## Relative risk ----
    # Stratified CMH test for point estimate and p-value,
    # and Wald statistic for the confidence interval
    lyt |>
      analyze(
        vars = vars,
        afun = a_relative_risk,
        table_names = paste0(vars, "est_relrisk_strat"),
        show_labels = "hidden",
        na_str = default_na_str(),
        extra_args = list(
          variables = list(strata = strata, arm = trtvar),
          conf_level = conflvl,
          .stats = c("rel_risk_ci", "pval"),
          .formats = c(rel_risk_ci = formats$ci, pval = formats$pval),
          .labels = c(
            rel_risk_ci = paste0(
              "Relative Risk (",
              tern::f_conf_level(conflvl),
              ")~[super a]"
            ),
            pval = "p-value~[super b]"
          ),
          ref_path = ref_path
        )
      )
  } else if (method == "or_logistic") {
    ## Stratified logistic regr. ----
    # Variant 1:
    # Stratified logistic regression for odds ratio point estimate,
    # p-value and confidence interval.
    lyt |>
      analyze(
        vars = vars,
        afun = a_odds_ratio_j,
        table_names = paste0(vars, "est_or_strat"),
        show_labels = "hidden",
        na_str = default_na_str(),
        extra_args = list(
          conf_level = conflvl,
          variables = list(strata = strata, arm = trtvar),
          method = "exact",
          na_if_no_events = TRUE,
          .stats = c("or_ci", "pval"),
          .formats = c(or_ci = formats$ci, pval = formats$pval),
          .labels = c(
            or_ci = paste0(
              "Odds Ratio (",
              tern::f_conf_level(conflvl),
              ")~[super a]"
            ),
            pval = "p-value~[super b]"
          ),
          ref_path = ref_path
        )
      )
  } else {
    ## Stratified CMH ----
    # Variant 2:
    # Stratified CMH test for odds ratio point estimate and confidence interval
    # and p-value.
    lyt |>
      analyze(
        vars = vars,
        afun = a_odds_ratio_j,
        table_names = paste0(vars, "est_or_strat_cmh"),
        show_labels = "hidden",
        na_str = default_na_str(),
        extra_args = list(
          conf_level = conflvl,
          variables = list(strata = strata, arm = trtvar),
          method = "cmh",
          na_if_no_events = TRUE,
          .stats = c("or_ci", "pval"),
          .formats = c(or_ci = formats$est_ci, pval = formats$pval),
          .labels = c(
            or_ci = paste0(
              "Odds Ratio (",
              tern::f_conf_level(conflvl),
              ")~[super a]"
            ),
            pval = "p-value~[super b]"
          ),
          ref_path = ref_path
        )
      )
  }
  lyt
}

# Layout ----

lyt <- rtables::basic_table(
  top_level_section_div = " ",
  show_colcounts = TRUE,
  colcount_format = "N=xx"
) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  append_topleft("Response") |>
  split_rows_by(
    subgroup,
    page_by = TRUE,
    split_fun = drop_split_levels
  ) |>
  summarize_row_counts(label_fstr = paste0(subgrlbl, ": %s")) |>
  insert_blank_line() |>
  ## Response Category ----
  analyze(
    vars = "rsp_lab",
    afun = s_proportion_factor,
    show_labels = "visible",
    var_labels = "Response category"
  ) |>
  insert_blank_line() |>
  ## CR+PR ----
  analyze_response(
    vars = "is_rsp_cr_pr",
    var_labels = "Overall response (CR+PR)",
    method = method
  ) |>
  insert_blank_line() |>
  ## CR ----
  analyze_response(
    vars = "is_rsp_cr",
    var_labels = "Complete response (CR)",
    method = method
  ) |>
  insert_blank_line() |>
  ## DCR (CR+PR+SD) ----
  analyze_response(
    vars = "is_rsp_cr_pr_sd",
    var_labels = "Disease control rate (CR+PR+SD)",
    method = method
  )

# Output ----

result <- build_table(lyt, df = ana, alt_counts_df = adsl)


# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
