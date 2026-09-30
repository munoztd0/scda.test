###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefpasi02.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefpasi02: PASI Response Analysis by Subgroup
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     adsl.sas7bdat, adparspi.sas7bdat
## Output:                    tefpasi02.rtf
## Remarks:                   Note: PASI50 response category not available in sample
##                            data. Therefore, this program summarizes PASI 75, 90, 100.
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

tblid <- "TEFPASI02"
fileid <- write_path(opath, tblid)
popfl <- "FASFL"
trtvar <- "TRT01P"

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

timepoints <- c("Week 1", "Week 2", "Week 4", "Week 8", "Week 12", "Week 16")
subgrpvar <- "WGTGR1"
subgrptxt <- "Subgroup: Weight"

################################################################################
# Process data:
################################################################################

# Read in SAS dataset and convert to R dataframe.
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
    ),
    subgrp = factor(paste(subgrptxt, !!rlang::sym(subgrpvar)))
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
    subgrp
  )

# Read in SAS dataset and convert to R dataframe.
adpasi <- haven::read_sas(read_path(a_in, "adparspi.sas7bdat")) |>
  filter(
    !!rlang::sym(popfl) == "Y" &
      PARAMCD %in% c("PASI75P", "PASI90P", "PASI100P") &
      AVISIT %in% timepoints
  ) |>
  mutate(
    AVISIT = forcats::fct_reorder(factor(AVISIT), AVISITN),
    response = factor(
      case_when(
        PARAMCD == "PASI75P" &
          AVALC == "Y" ~
          ">= 75% improvement from baseline",
        PARAMCD == "PASI90P" &
          AVALC == "Y" ~
          ">= 90% improvement from baseline",
        PARAMCD == "PASI100P" & AVALC == "Y" ~ "100% improvement from baseline",
        TRUE ~ "Nonresponder"
      ),
      levels = c(
        "Nonresponder",
        ">= 75% improvement from baseline",
        ">= 90% improvement from baseline",
        "100% improvement from baseline"
      )
    )
  ) |>
  select(USUBJID, AVISIT, response)

adpasi <- inner_join(x = adsl, y = adpasi, by = "USUBJID")

################################################################################
# Define layout and build table:
################################################################################

# Map each treatment group to appropriate columns spanning header.
colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = "PLACEBO",
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

lyt <- basic_table(show_colcounts = TRUE, colcount_format = "N=xx") %>%
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) %>%
  split_cols_by(trtvar) %>%
  split_rows_by("subgrp", section_div = c(" ")) %>%
  summarize_row_groups(
    "subgrp",
    cfun = a_freq_j,
    extra_args = list(
      denom = "n_altdf",
      denom_by = "subgrp",
      riskdiff = FALSE,
      extrablankline = TRUE,
      .stats = c("n_altdf")
    )
  ) %>%
  split_rows_by(
    "AVISIT",
    split_fun = keep_split_levels(only = timepoints),
    section_div = " "
  ) %>%
  analyze(
    vars = "response",
    afun = a_freq_j,
    extra_args = list(
      denom = "n_df",
      .stats = c("n_df", "count_unique_fraction"),
      val = c(
        ">= 75% improvement from baseline",
        ">= 90% improvement from baseline",
        "100% improvement from baseline"
      )
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
