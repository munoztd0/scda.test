###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              teffpga01.r
## R Version:                 4.5.2
## junco Version:             0.1.6
## Short Description:         Program to create teffpga01:
##                            Subjects Achieving [Response Description] at [Time Point]
##                            Among Subjects With [Baseline Criteria] at Baseline
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adbdc, adclrori
## Output:                    teffpga01.rtf
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

# Define output ID and file location.
tblid <- "TEFFPGA01"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Placebo to JNJ-77242113"

# Define the order of treatment group labels in the treatment variable.
trtlab <- c("JNJ-77242113")

# Define population flags used.
popfl <- "FASFL"

# Define the stratification variables to use for CMH.
stratvar <- c("STRAT02")

# Define visit at which to analyze the results.
timepoint <- "Week 16"

# Define response parameters to be used and analyzed separately,
# along with their labels to use in the table.
resppar <- c(
  "FPCLRP" ~ paste0("Subjects achieving an f-PGA score of 0 at ", timepoint),
  "FPMINP" ~ paste0("Subjects achieving an f-PGA score of 0 or 1 at ", timepoint)
)

# Define the baseline criteria.
basefl <- quote(PARAMCD == "FPGABL" & AVAL >= 2)

# The label to describe baseline criteria in the table.
baselbl <- "an f-PGA score >= 2"

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived formats specifications.
formats <- list(
  est_prop = jjcsformat_count_fraction,
  est_prop_diff = jjcsformat_xx("xx.x% (xx.x%, xx.x%)"),
  pval = jjcsformat_pval_fct(alpha)
)

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  mutate(!!trtvar := factor(.data[[trtvar]], levels = c(ctrlab, trtlab))) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(stratvar))

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

## ADBDC ----

adbdc <- haven::read_sas(read_path(a_in, "adbdc.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  # Only take patients fulfilling baseline criteria here.
  filter(!!basefl) |>
  select(STUDYID, USUBJID)

## ADCLRORI ----

adclrori <- haven::read_sas(read_path(a_in, "adclrori.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(
    PARAMCD %in% sapply(resppar, junco:::leftside),
    AVISIT == timepoint
  ) |>
  select(STUDYID, USUBJID, PARAMCD, AVALC) |>
  mutate(
    USUBJID = factor(USUBJID),
    PARAMLBL = case_match(
      PARAMCD,
      !!!resppar
    ),
    PARAMCD = factor(PARAMCD, levels = sapply(resppar, junco:::leftside))
  )

## Analysis ----

ana <- adbdc |>
  # Inner join here because we only want to look at patients fulfilling baseline
  # criteria.
  inner_join(adsl, by = c("STUDYID", "USUBJID")) |>
  # Left join here because missing patients would be counted as non-responders.
  left_join(adclrori, by = c("STUDYID", "USUBJID")) |>
  mutate(
    response = case_when(
      AVALC == "Y" ~ TRUE,
      .default = FALSE
    )
  )

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
  split_rows_by(
    "PARAMCD",
    section_div = " ",
    labels_var = "PARAMLBL"
  ) |>
  summarize_row_groups(
    "response",
    cfun = function(df, labelstr, .var) {
      in_rows(
        .list = list(
          nrow(df),
          s_proportion(df, .var = .var, conf_level = conflvl)$n_prop
        ),
        .labels = c(
          paste("Subjects with", baselbl, "at baseline"),
          labelstr
        ),
        .formats = list("xx", formats$est_prop),
        .indent_mods = list(0, 1)
      )
    }
  ) |>
  analyze(
    "response",
    afun = a_proportion_diff_j,
    na_str = default_na_str(),
    table_names = "est_prop_diff",
    show_labels = "hidden",
    extra_args = list(
      .stats = "diff_est_ci",
      .labels = c(
        diff_est_ci = paste0(
          "% Difference (",
          tern::f_conf_level(conflvl),
          ")"
        )
      ),
      .formats = c(diff_est_ci = formats$est_prop_diff),
      .indent_mods = 1,
      method = "cmh_mn",
      variables = list(strata = stratvar),
      conf_level = conflvl,
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
      .formats = c("pval" = formats$pval),
      .indent_mods = 2,
      ref_path = ref_path,
      .labels = c(pval = "p-value")
    )
  )

# Output ----

result <- build_table(lyt, ana, alt_counts_df = adsl, round_type = "sas")

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
