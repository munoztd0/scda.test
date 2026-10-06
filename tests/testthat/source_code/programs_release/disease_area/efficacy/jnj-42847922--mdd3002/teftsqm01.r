###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              teftsqm01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create teftsqm01: Treatment Satisfaction Questionnaire of Medication
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADQS.
## Output:                    TEFTSQM01.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:                     1
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
## Description:
################################################################################

################################################################################
# Prep Environment
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

################################################################################
# - Define output ID and file location
# - Define treatment variable used (default=TRT01P)
# - Define population flag used (default=SDBFL)
# - Define control treatment arm
################################################################################

tblid <- "TEFTSQM01"
fileid <- write_path(opath, tblid)

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

trtvar <- "TRT01P"
popfl <- "FAS2FL"

ctrl_grp <- "Placebo"

################################################################################
# Read in ADSL and ADxxx datasets
################################################################################

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) %>%
  filter(!!rlang::sym(popfl) == "Y" & !!rlang::sym(trtvar) != "") %>%
  mutate(!!trtvar := as.factor(.data[[trtvar]])) %>%
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl))

adsl$colspan_trt <- factor(
  ifelse(adsl[[trtvar]] == ctrl_grp, " ", "Active Study Agent"),
  levels = c("Active Study Agent", " ")
)

################################################################################
# Use adphq9 for now until we find out where to get this data from - Remove later
################################################################################

ADQS <- haven::read_sas(read_path(a_in, "adphq9.sas7bdat")) %>%
  mutate(PARAMCD = as.factor(PARAMCD)) %>%
  group_by(USUBJID, PARAMCD) %>%
  slice(1) %>%
  ungroup() %>%
  select(USUBJID, PARAMCD, AVAL)

adqs <- ADQS %>%
  mutate(
    AVAL1 = case_when(PARAMCD == "PHQ9TOT" ~ AVAL),
    AVAL2 = case_when(PARAMCD == "PHQ0101" ~ AVAL),
    AVAL3 = case_when(PARAMCD == "PHQ0102" ~ AVAL)
  ) %>%
  select(USUBJID, starts_with("AVAL")) %>%
  select(-AVAL)

# join data together
qs <- adqs %>% inner_join(., adsl, by = c("USUBJID"))

colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = ctrl_grp,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

################################################################################
# Define layout and build table:
################################################################################

lyt <- rtables::basic_table(
  top_level_section_div = " ",
  show_colcounts = TRUE,
  colcount_format = "N=xx"
) %>%
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) %>%
  split_cols_by(trtvar) %>%
  analyze(
    "AVAL1",
    table_names = "AVAL1x",
    var_labels = "Global satisfaction score",
    show_labels = "visible",
    indent_mod = 0L,
    afun = function(x) {
      list(
        "N" = rcell(length(x), format = jjcsformat_xx("xx"))
      )
    }
  ) %>%
  analyze(
    "AVAL1",
    nested = TRUE,
    var_labels = "Global satisfaction score",
    show_labels = "hidden",
    indent_mod = 2L,
    afun = function(x) {
      list(
        "Mean (SD)" = rcell(
          c(mean(x), sd(x)),
          format = jjcsformat_xx("xx.x (xx.xx)")
        ),
        "Median" = rcell(median(x), format = jjcsformat_xx("xx.x")),
        "Min, max" = rcell(
          c(min(x), max(x)),
          format = jjcsformat_xx("xx., xx.")
        )
      )
    }
  ) %>%
  analyze(
    "AVAL2",
    nested = FALSE,
    table_names = "AVAL2x",
    var_labels = "Effectiveness score",
    show_labels = "visible",
    indent_mod = 0L,
    afun = function(x) {
      list(
        "N" = rcell(length(x), format = jjcsformat_xx("xx"))
      )
    }
  ) %>%
  analyze(
    "AVAL2",
    nested = TRUE,
    var_labels = "Effectiveness score",
    show_labels = "hidden",
    indent_mod = 2L,
    afun = function(x) {
      list(
        "Mean (SD)" = rcell(
          c(mean(x), sd(x)),
          format = jjcsformat_xx("xx.x (xx.xx)")
        ),
        "Median" = rcell(median(x), format = jjcsformat_xx("xx.x")),
        "Min, max" = rcell(
          c(min(x), max(x)),
          format = jjcsformat_xx("xx., xx.")
        )
      )
    }
  ) %>%
  analyze(
    "AVAL3",
    nested = FALSE,
    table_names = "AVAL3x",
    var_labels = "Convenience score",
    show_labels = "visible",
    indent_mod = 0L,
    afun = function(x) {
      list(
        "N" = rcell(length(x), format = jjcsformat_xx("xx"))
      )
    }
  ) %>%
  analyze(
    "AVAL3",
    nested = TRUE,
    var_labels = "Convenience score",
    show_labels = "hidden",
    indent_mod = 2L,
    afun = function(x) {
      list(
        "Mean (SD)" = rcell(
          c(mean(x), sd(x)),
          format = jjcsformat_xx("xx.x (xx.xx)")
        ),
        "Median" = rcell(median(x), format = jjcsformat_xx("xx.x")),
        "Min, max" = rcell(
          c(min(x), max(x)),
          format = jjcsformat_xx("xx., xx.")
        )
      )
    }
  )

result <- build_table(lyt, qs, alt_counts_df = adsl)

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

################################################################################
# Convert to tbl file and output table
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
