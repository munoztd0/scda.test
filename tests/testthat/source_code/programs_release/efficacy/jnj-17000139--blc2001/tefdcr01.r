###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

## Original Reporting Effort: Standards
## Program Name:              tefdcr01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefdcr01: Disease Control Rate With Stable Disease
##                            at 3 Months and 6 Months– Independent Radiographic Review -
##                            [Stratified/Unstratified] Analysis; Full Analysis Set (Study mmy, bc, lc)
## Author:                    Technology Solutions
## Date:                      2026-09-30
## Input:                     ADSL, ADEFF
## Output:                    TEFDCR01.rtf
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
tblid <- "TEFDCR01"
fileid <- write_path(opath, tblid)


tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Please choose appropriate title and footnotes.
tab_titles$title <- tab_titles$title[1]
tab_titles$main_footer <- tab_titles$main_footer[c(1, 2)]

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "CETRELIMAB"

# Define population flag used.
popfl <- "FASFL"

# Define the strata variables to use below.
# strata <- c("SEX")
# For unstratified analyses, please use:
strata <- NULL

# Define response parameter to be used.
resppar <- "BORSPIRC"

# Define method for comparing response.
method <- "or_logistic"
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
set_default_na_str("NE")

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl), all_of(strata))

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
    rsp_lab = forcats::fct_collapse(
      AVALC,
      "Responder" = "CR",
      "Non-responder" = c("NE", "NON-CR")
    ),
    is_resp = (rsp_lab == "Responder")
  )

## Analysis ----

ana <- adeff |>
  inner_join(adsl, by = "USUBJID")

# Functions ----

# Layout generating function for the response variable analysis.
analyze_response <- function(
  lyt,
  vars,
  method = c("rr", "or_logistic", "or_cmh")
) {
  method <- match.arg(method)

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
    ## Stratified (or unstratified) logistic regr. ----
    # Variant 1:
    # Stratified logistic regression for odds ratio point estimate,
    # p-value and confidence interval.
    lyt |>
      analyze(
        vars = vars,
        afun = a_odds_ratio_j,
        table_names = paste0(vars, "est_or_strat"),
        na_str = default_na_str(),
        show_labels = "hidden",
        extra_args = list(
          conf_level = conflvl,
          variables = list(strata = strata, arm = trtvar),
          method = "exact", # Note: Only used when strata are present.
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
          .indent_mods = c(
            or_ci = 0L,
            pval = 0L
          ),
          ref_path = ref_path
        )
      )
  } else {
    ## Stratified CMH ----
    # Variant 2:
    # Stratified CMH test for odds ratio point estimate and confidence interval, but p-value
    # from Fisher’s exact test.
    lyt |>
      analyze(
        vars = vars,
        afun = a_odds_ratio_j,
        table_names = paste0(vars, "est_or_strat_cmh"),
        na_str = default_na_str(),
        show_labels = "hidden",
        extra_args = list(
          conf_level = conflvl,
          variables = list(strata = strata, arm = trtvar),
          method = "cmh",
          na_if_no_events = TRUE,
          .stats = "or_ci",
          .formats = c(or_ci = formats$est_ci),
          .labels = c(
            or_ci = paste0(
              "Odds Ratio (",
              tern::f_conf_level(conflvl),
              ")~[super a]"
            )
          ),
          ref_path = ref_path
        )
      ) |>
      ## Fisher's exact test ----
      analyze(
        vars = vars,
        afun = a_test_proportion_diff,
        table_names = paste0(vars, "test_prop_strat"),
        show_labels = "hidden",
        na_str = default_na_str(),
        extra_args = list(
          variables = list(strata = strata),
          method = "fisher",
          ref_path = ref_path,
          .formats = c(pval = formats$pval)
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
  ## Response Category ----
  analyze(
    vars = "rsp_lab",
    afun = s_proportion_factor,
    show_labels = "hidden"
  ) |>
  insert_blank_line() |>
  ## Response Analysis ----
  analyze_response(
    vars = "is_resp",
    method = method
  )

# Output ----

result <- build_table(lyt, df = ana, alt_counts_df = adsl)


# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
