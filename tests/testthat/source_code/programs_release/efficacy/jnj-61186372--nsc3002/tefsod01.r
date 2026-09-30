###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefsod01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefsod01: Maximum Maximum Percent Change From Baseline
##                            in the Sum of Diameters for All Target Lesions; Full Analysis Set
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADRECIST
## Output:                    TEFSOD01.rtf
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
tblid <- "TEFSOD01"
fileid <- write_path(opath, tblid)

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
# Warning: Title file should contain exactly one title record per Table ID
# Therefore need to make sure we only have one title record:
tab_titles$title <- tab_titles$title[1]

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "CP"

# Define population flags used.
popfl <- c("FASFL", "MDBIRCFL")

# Define the PARAMCD value defining the percent change for analysis.
resppar <- "RTRGSUM"

# Define additional analysis flags to be used.
anlfl <- c("SAFFL", "ANL05FL", "ANL06FL")

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl))

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

## ADRECIST ----

adrecist <- haven::read_sas(read_path(a_in, "adrecist.sas7bdat")) |>
  filter(if_all(all_of(anlfl), ~ .x == "Y")) |>
  select(STUDYID, USUBJID, PARAMCD, PCHG, all_of(trtvar), all_of(anlfl))

adrecist <- adrecist |>
  select(USUBJID, PARAMCD, PCHG) |>
  filter(PARAMCD == resppar) |>
  filter(!is.na(PCHG)) |>
  select(-PARAMCD)

checkmate::assert_true(!any(duplicated(adrecist$USUBJID)))

## Analysis ----

categories <- c(
  "<= -30",
  "> -30 to <= -10",
  "> -10 to < 20",
  ">= 20"
)
ana <- adrecist |>
  left_join(adsl, by = "USUBJID") |>
  mutate(
    AVALC = case_when(
      PCHG <= -30 ~ categories[1],
      PCHG > -30 & PCHG <= -10 ~ categories[2],
      PCHG > -10 & PCHG < 20 ~ categories[3],
      PCHG >= 20 ~ categories[4],
    ),
    AVALC = factor(AVALC, levels = categories)
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
  analyze(
    vars = "AVALC",
    afun = s_proportion_factor,
    show_labels = "visible",
    var_labels = "Maximum % change in sum of diameters"
  ) |>
  append_topleft("Parameter")

# Output ----

result <- build_table(lyt, ana, alt_counts_df = adsl)

# Add titles and footnotes.

result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
