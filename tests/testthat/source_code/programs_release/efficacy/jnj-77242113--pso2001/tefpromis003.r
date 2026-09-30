###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefpromis03.R
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefpromis003:
##                            TEFPROMIS003: Subjects Achieving ≥5-point Improvement From Baseline in
##                            PROMIS-29 Physical and Mental Component Scores Through [Time Point];
##                            Full Analysis Set (Study psoriasis)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADPRRSPI
## Output:                    TEFPROMIS003.rtf
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
tblid <- "TEFPROMIS003"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "PLACEBO"

# Define the order of treatment group labels in the treatment variable.
trtlab <- c(
  "JNJ-77242113 25 MG QD",
  "JNJ-77242113 50 MG QD",
  "JNJ-77242113 25 MG BID",
  "JNJ-77242113 100 MG QD",
  "JNJ-77242113 100 MG BID"
)

# Define population flags used.
popfl <- "FASFL"

# Define the stratification variables to use for CMH.
strata <- c("STRATWTG")

# Define visits at which to analyze the results.
timepoints <- c("Week 8", "Week 16")

# Define response parameters to be used and analyzed separately,
# along with their labels to use in the table.
resppar <- c(
  "PRPCGE5P" ~ "PROMIS-29 - PCS (T-score)",
  "PRMCGE5P" ~ "PROMIS-29 - MCS (T-score)"
)

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
  select(STUDYID, USUBJID, all_of(trtvar), all_of(strata))

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

## ADPRRSPI ----

adprrspi <- haven::read_sas(read_path(a_in, "adprrspi.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(
    PARAMCD %in% sapply(resppar, junco:::leftside),
    AVISIT %in% timepoints
  ) |>
  select(USUBJID, AVISIT, PARAMCD, AVALC) |>
  mutate(
    USUBJID = factor(USUBJID),
    PARAMLBL = case_match(
      PARAMCD,
      !!!resppar
    ),
    PARAMCD = factor(PARAMCD, levels = sapply(resppar, junco:::leftside)),
    AVISIT = factor(AVISIT, levels = timepoints)
  )

## Analysis ----

ana <- adprrspi |>
  inner_join(adsl, by = "USUBJID") |>
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
  split_rows_by(
    "AVISIT",
    section_div = " "
  ) |>
  summarize_row_groups(
    "response",
    cfun = c_proportion_logical,
    extra_args = list(
      label_fstr = "Subjects achieving >=5-point improvement from baseline at %s",
      format = formats$est_prop
    )
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
          "Treatment difference (",
          tern::f_conf_level(conflvl),
          ")"
        )
      ),
      .formats = c(diff_est_ci = formats$est_prop_diff),
      method = "cmh",
      variables = list(strata = strata),
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
      variables = list(strata = strata),
      .formats = c("pval" = formats$pval),
      ref_path = ref_path,
      .labels = c(pval = "p-value")
    )
  )


# Output ----

result <- build_table(lyt, ana, alt_counts_df = adsl)

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = 
  result,
  file = fileid,
  fontspec = formatters::font_spec("Times", 8L, 1.2),
  orientation = "landscape"
)
