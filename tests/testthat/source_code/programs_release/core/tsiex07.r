###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsiex07.r
## R version:                 4.5.2
## junco version:             0.1.3
## Short Description:         Program to create tsiex07: Dose Levels Over Time
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adex
## Output:                    tsiex07.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
##  Rev #:
##  Modified By:
##  Reporting Effort:
##  Date:
##  Description:
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

################################################################################
# Define script level parameters:
################################################################################

################################################################################
# - Define output ID and file location
# - Define treatment variable used (default=TRT01A)
# - Define population flag used (default=SAFFL)
# - Define the column ordering of the dose levels
################################################################################

tblid <- "TSIEX07"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()


trtvar <- "TRT01A"
popfl <- "SAFFL"

dose_order <- c("7.5", "15", "Total")
################################################################################
# Process Data:
################################################################################

# Read in required data
adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) |>
  df_na() |>
  filter(!!rlang::sym(popfl) == "Y" & !!rlang::sym(trtvar) != "Placebo") |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl)) |>
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c(
        "Xanomeline Low Dose",
        "Xanomeline High Dose"
      )
    )
  )

adex <- haven::read_sas(envsetup::read_path(a_in, "adex.sas7bdat")) |>
  df_na() |>
  filter(!is.na(ADOSE) & !is.na(AVISIT) & !grepl("SCREENING", AVISIT, ignore.case = TRUE)) |>
  mutate(ADOSEC = as.factor(ADOSE)) |>
  select(STUDYID, USUBJID, ADOSEC, AVISIT, AVISITN) |>
  mutate(
    AVISIT := factor(
      stringr::str_to_sentence(as.character(AVISIT)),
      levels = stringr::str_to_sentence(unique(.data[['AVISIT']])[order(unique(.data[['AVISITN']]))])
    )
  )

# Create Total column manually as we cannot use the table layout option since patients can
# appear in multiple dose level columns

extot_ <- adex |>
  group_by(USUBJID, AVISIT) |>
  slice(1) |>
  mutate(ADOSEC = "Total") |>
  ungroup()

extot <- bind_rows(adex, extot_)

# join data together
ex <- extot |> inner_join(adsl, by = c("STUDYID", "USUBJID"))

ex$colspan_trt <- factor(
  ifelse(ex$ADOSEC == "Total", " ", "Active Study Agent"),
  levels = c("Active Study Agent", " ")
)


ex$ADOSEC <- factor(ex$ADOSEC, levels = dose_order)

# Create dataset to be used for column counts
adex_unique <- ex |>
  group_by(USUBJID, ADOSEC) |>
  slice(1) |>
  ungroup()

################################################################################
# Define layout and build table:
################################################################################

lyt <- rtables::basic_table(
  top_level_section_div = " ",
  show_colcounts = TRUE,
  colcount_format = "N=xx"
) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_in_group("ADOSEC")
  ) |>
  split_cols_by("ADOSEC") |>
  analyze(
    "AVISIT",
    afun = a_freq_j,
    show_labels = "hidden",
    extra_args = list(
      .stats = "count_unique_fraction",
      denom = "N_col"
    )
  ) |>
  append_topleft("Time Point, n (%)")

result <- build_table(lyt, ex, alt_counts_df = adex_unique, round_type = "sas")


################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

################################################################################
# Convert to tbl file and output table
################################################################################
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
