###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefcgis08.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefcgis08: Frequency Distribution Over
##                            Time[ - Double-blind Phase]; Full Analysis Set Analysis Set (Study mdd)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADCGIS
## Output:                    TEFCGIS08.rtf
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
tblid <- "TEFCGIS08"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Placebo"

# Define population flags used.
popfl <- "FAS1FL"

# Define analysis domain flags to be used (for ADCGIS).
domfl <- quote(ANL02FL == "Y" & APHASE == "Double Blind Phase")

# Define response parameter to be used (for ADMADRSI).
resppar <- "CGI0101"

# Mapping of response values to labels.
respvals <- c(
  "1" ~ "Normal (not at all ill)",
  "2" ~ "Borderline ill",
  "3" ~ "Mildly ill",
  "4" ~ "Moderately ill",
  "5" ~ "Markedly ill",
  "6" ~ "Severely ill",
  "7" ~ "Among the most extremely ill patients"
)

# Define time points to be used.
timepts <- c("Baseline (DB)", "Day 15", "Day 29", "Day 43", "Endpoint (DB)")

# Derived formats specifications.
formats <- list(
  n = "xx",
  percent = jjcsformat_xx("xx.x"),
  cum_percent = jjcsformat_xx("xx.x")
)

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  mutate(!!trtvar := forcats::fct_relevel(!!rlang::sym(trtvar), ctrlab)) |>
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

## ADCGIS ----

adcgis <- haven::read_sas(read_path(a_in, "adcgis.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(PARAMCD == resppar) |>
  filter(!!domfl) |>
  select(USUBJID, AVISIT, AVAL) |>
  mutate(
    USUBJID = factor(USUBJID),
    AVISIT = factor(AVISIT),
    RES = case_match(
      as.character(AVAL),
      !!!respvals,
      .ptype = factor(levels = sapply(respvals, rlang::f_rhs))
    )
  )

## Analysis ----

ana <- adcgis |>
  filter(AVISIT %in% timepts) |>
  inner_join(adsl, by = "USUBJID")

# Layout ----

lyt <- rtables::basic_table(top_level_section_div = " ") |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar, show_colcounts = TRUE, colcount_format = "N=xx") |>
  split_rows_by(
    "AVISIT",
    label_pos = "visible",
    split_label = "Severity of illness",
    section_div = " "
  ) |>
  split_cols_by(
    "STUDYID", # The particular choice of this is irrelevant.
    split_fun = prop_split_fun
  ) |>
  analyze("RES", afun = prop_table_afun, extra_args = list(formats = formats))

# Output ----

result <- build_table(lyt, df = ana, alt_counts_df = adsl)


# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
