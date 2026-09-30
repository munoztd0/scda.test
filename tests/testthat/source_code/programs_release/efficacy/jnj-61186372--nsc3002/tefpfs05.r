###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefpfs05.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefpfs05: Progression-free Survival
##                            -Concordance Between Investigator and BICR
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     adsl.sas7bdat, adtteff.sas7bdat
## Output:                    tefpfs05.rtf
## Remarks:                   This program will need to be updated once we have
##                            JJCS-compliant data and analysis rules that are either
##                            applied in the data or documented in the DPS. Code is
##                            currently based on corresponding SAS code in lung cancer trial.
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

tblid <- "TEFPFS05"
fileid <- write_path(opath, tblid)
popfl <- "FASFL"
trtvar <- "TRT01P"

titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
titles$title <- titles$title[1]

################################################################################
# Process data:
################################################################################

# Read in SAS dataset and convert to R dataframe.
adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  mutate(
    !!popfl := factor(!!rlang::sym(popfl)),
    !!trtvar := factor(
      !!rlang::sym(trtvar),
      levels = c("ACP", "LACP/ACP-L", "CP")
    )
  ) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  create_colspan_var(
    non_active_grp = "CP",
    non_active_grp_span_lbl = " ",
    active_grp_span_lbl = "Active Study Agent",
    colspan_var = "colspan_trt",
    trt_var = trtvar
  ) |>
  select(USUBJID, !!rlang::sym(popfl), !!rlang::sym(trtvar), colspan_trt)

adtteef <- haven::read_sas(read_path(a_in, "adtteef.sas7bdat")) |>
  mutate(EVNTDESC = factor(EVNTDESC)) |>
  select(USUBJID, EVNTDESC, ADT, ADY, AVAL, PARAMCD, PARAM, CNSR)

adtteer <- adtteef |>
  filter(PARAMCD %in% c("RPFS"))

adtteev <- adtteef |>
  filter(PARAMCD %in% c("VPFS"))

adttee <- full_join(x = adtteer, y = adtteev, by = "USUBJID")

# Following rule based upon SAS code in Lung-cancer trial.
# If the difference in the progression-free survival (months/30.4375) according
# to IRC and the investigator is less than 7, then Status is "Complete agreement".

adttee <- inner_join(x = adsl, y = adttee, by = "USUBJID") |>
  mutate(
    pfsbodsys = factor(
      ifelse(
        CNSR.x == CNSR.y,
        "Agreement on event status~[super a]",
        "Disagreement on event status"
      ),
      levels = c(
        "Agreement on event status~[super a]",
        "Disagreement on event status"
      )
    ),
    pfsdecod1 = factor(
      case_when(
        CNSR.x == 0 & CNSR.y == 0 ~ "Event by both IRC and investigator",
        CNSR.x == 1 & CNSR.y == 1 ~ "No event by both IRC and investigator",
        CNSR.x == 1 & CNSR.y == 0 ~ "Event by investigator but not by IRC",
        CNSR.x == 0 & CNSR.y == 1 ~ "Event by IRC but not by investigator",
        TRUE ~ NA
      ),
      levels = c(
        "Event by both IRC and investigator",
        "No event by both IRC and investigator",
        "Event by investigator but not by IRC",
        "Event by IRC but not by investigator"
      )
    ),
    pfsdecod2 = factor(
      case_when(
        CNSR.x == 0 &
          CNSR.y == 0 &
          abs(AVAL.x * 30.4375 - AVAL.y * 30.4375) < 7 ~
          "Complete agreement",
        CNSR.x == 0 &
          CNSR.y == 0 &
          ADY.x > ADY.y ~
          "Agreement with later date by IRC",
        CNSR.x == 0 &
          CNSR.y == 0 &
          ADY.x < ADY.y ~
          "Agreement with earlier date by IRC",
        TRUE ~ NA
      ),
      levels = c(
        "Complete agreement",
        "Agreement with later date by IRC",
        "Agreement with earlier date by IRC"
      )
    )
  )

################################################################################
# Define layout and build table:
################################################################################

colspan_trt_map <- create_colspan_map(
  df = adsl,
  non_active_grp = "CP",
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

lyt <- basic_table(
  show_colcounts = TRUE,
  colcount_format = "N=xx",
  top_level_section_div = " "
) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  split_rows_by("pfsbodsys", section_div = " ") |>
  summarize_row_groups(
    "pfsbodsys",
    cfun = a_freq_j,
    extra_args = list(.stats = c("count_unique_fraction"))
  ) |>
  analyze(
    vars = "pfsdecod1",
    table_names = "event_by_both",
    show_labels = "hidden",
    afun = a_freq_j,
    extra_args = list(
      .stats = c("count_unique_fraction"),
      val = "Event by both IRC and investigator",
      label = "Event by both IRC and investigator"
    )
  ) |>
  analyze(
    vars = "pfsdecod2",
    show_labels = "hidden",
    afun = a_freq_j,
    extra_args = list(.stats = c("count_unique_fraction")),
    indent_mod = 1
  ) |>
  analyze(
    vars = "pfsdecod1",
    table_names = "no_event_by_both",
    show_labels = "hidden",
    afun = a_freq_j,
    extra_args = list(
      .stats = c("count_unique_fraction"),
      val = "No event by both IRC and investigator",
      label = "No event by both IRC and investigator"
    )
  ) |>
  analyze(
    vars = "pfsdecod1",
    table_names = "event_by_inv",
    show_labels = "hidden",
    afun = a_freq_j,
    extra_args = list(
      .stats = c("count_unique_fraction"),
      val = "Event by investigator but not by IRC",
      label = "Event by investigator but not by IRC"
    )
  ) |>
  analyze(
    vars = "pfsdecod1",
    table_names = "event_by_irc",
    show_labels = "hidden",
    afun = a_freq_j,
    extra_args = list(
      .stats = c("count_unique_fraction"),
      val = "Event by IRC but not by investigator",
      label = "Event by IRC but not by investigator"
    )
  )

result <- build_table(lyt, adttee, alt_counts_df = adsl)

################################################################################
# Post-Processing:
# - Prune any categories with all zeros.
################################################################################

result <- safe_prune_table(result, prune_func = count_pruner(cols = NULL))

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, titles)
################################################################################
# Convert to tbl file and output table:
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
