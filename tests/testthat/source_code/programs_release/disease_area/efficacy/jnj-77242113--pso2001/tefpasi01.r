###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefpasi01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefpasi01: PASI Response Analysis
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     adsl.sas7bdat, adparspi.sas7bdat
## Output:                    tefpasi01.rtf
## Remarks:
##
## Modification History:
## Rev #:
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
## Description:
################################################################################

################################################################################
# Prep environment:
################################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)


library(dplyr)
library(rtables)
library(junco)
library(haven)

################################################################################
# Define script level parameters:
################################################################################

tblid <- "TEFPASI01"
fileid <- write_path(opath, tblid)
popfl <- "FASFL"
trtvar <- "TRT01P"

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

timepoint <- "Week 16"
resp_paramcd <- "PASI75P"
resp_txt <- "PASI 75 responders"
stratvar <- "STRATWTG"

################################################################################
# Process data:
################################################################################

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(
    !!popfl := factor(!!rlang::sym(popfl)),
    !!trtvar := factor(
      !!rlang::sym(trtvar),
      levels = c(
        "JNJ-77242113 25 MG QD",
        "JNJ-77242113 50 MG QD",
        "JNJ-77242113 25 MG BID",
        "JNJ-77242113 100 MG QD",
        "JNJ-77242113 100 MG BID",
        "PLACEBO"
      )
    )
  ) |>
  create_colspan_var(
    non_active_grp = "PLACEBO",
    non_active_grp_span_lbl = " ",
    active_grp_span_lbl = "Active Study Agent",
    colspan_var = "colspan_trt",
    trt_var = trtvar
  ) |>
  select(
    USUBJID,
    !!rlang::sym(popfl),
    !!rlang::sym(trtvar),
    colspan_trt,
    !!rlang::sym(stratvar)
  )

adpasi <- haven::read_sas(read_path(a_in, "adparspi.sas7bdat")) |>
  filter(
    !!rlang::sym(popfl) == "Y" & AVISIT == timepoint & PARAMCD == resp_paramcd
  ) |>
  mutate(
    AVISIT = forcats::fct_reorder(factor(AVISIT), AVISITN),
    response = case_when(
      AVALC == "Y" ~ TRUE,
      TRUE ~ FALSE
    )
  ) |>
  select(USUBJID, AVISIT, PARAMCD, AVALC, response)

adpasi <- inner_join(x = adsl, y = adpasi, by = "USUBJID")

################################################################################
# Define layout and build table:
################################################################################

# Map each treatment group to appropriate columns spanning header.
colspan_trt_map <- create_colspan_map(
  df = adsl,
  trt_var = trtvar,
  colspan_var = "colspan_trt",
  non_active_grp = "PLACEBO",
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent"
)

# Specification for reference column
ref_path <- c("colspan_trt", " ", trtvar, "PLACEBO")

lyt <- basic_table(show_colcounts = TRUE, colcount_format = "N=xx") |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  estimate_proportion(
    vars = "response",
    table_names = "est_prop",
    .stats = c("n_prop"),
    .labels = c("n_prop" = resp_txt),
    .formats = c("n_prop" = jjcsformat_count_fraction)
  ) |>
  analyze(
    vars = "response",
    afun = a_proportion_diff_j,
    show_labels = "hidden",
    na_str = default_na_str(),
    table_names = "diff_est_ci",
    extra_args = list(
      .stats = c("diff_est_ci"),
      .labels = c(diff_est_ci = "Treatment difference (95% CI)"),
      .formats = c(diff_est_ci = jjcsformat_xx("xx.x (xx.x, xx.x)")),
      .indent_mods = 1,
      method = "cmh",
      variables = list(strata = stratvar),
      ref_path = ref_path
    )
  ) |>
  analyze(
    vars = "response",
    afun = a_test_proportion_diff,
    table_names = "pval",
    show_labels = "hidden",
    na_str = default_na_str(),
    extra_args = list(
      method = "cmh",
      variables = list(strata = stratvar),
      ref_path = ref_path,
      .labels = c(pval = "p-value")
    )
  )

result <- build_table(lyt, adpasi, alt_counts_df = adsl)

################################################################################
# Post-Processing:
################################################################################

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)
################################################################################
# Convert to tbl file and output table:
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
