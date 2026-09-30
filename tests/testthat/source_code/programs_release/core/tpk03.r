###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tpk03.r
## R version:                 4.5.2
## junco version:             0.1.3
## Short Description:         Program to create tpk03:
##                            Subjects With [Matrix]
##                            [Active Study Agent/Analyte] Concentrations Below
##                            the Lowest Quantification Level ([LLOQ/LLOQxMRD
##                            and units]) Over Time
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adpc
## Output:                    tpk03.rtf
## Remarks:
##
## Modification History:
##  Rev #:
##  Modified By:
##  Reporting Effort:
##  Date:
##  Description:
################################################################################

################################################################################
# Prep environment:
################################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(dplyr)
library(forcats)
library(rtables)
library(tern)
library(junco)

################################################################################
# Define script level parameters:
################################################################################

tblid <- "TPK03"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
popfl <- "PKFL"
trtvar <- "TRT01A"
paramcd <- "XAN"

# flag to indicate if ATPT is present in the study
# if TRUE, time point is concatenation of AVISIT and ATPT
# if FALSE, time point is AVISIT only
use_atpt <- TRUE

# lloq_label: the threshold term to display in the label
# e.g. "LLOQxMRD" or "LLOQ"
lloq_label <- "LLOQxMRD"


################################################################################
# Process data:
################################################################################

adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) |>
  df_na() |>
  filter(.data[[popfl]] == "Y") |>
  select(STUDYID, USUBJID, all_of(c(trtvar, popfl))) |>
  # Drop the control group
  filter(.data[[trtvar]] != "Placebo") |>
  mutate({{ trtvar }} := fct_drop(.data[[trtvar]])) |>
  mutate(colspan_trt = "Active Study Agent") |>
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c(
        "Xanomeline Low Dose",
        "Xanomeline High Dose"
      )
    )
  )


# IF ATPT not present user can remove it from select statement

adpc <- haven::read_sas(envsetup::read_path(a_in, "adpc.sas7bdat")) |>
  df_na() |>
  filter(PARAMCD == paramcd) |> # If not required this filter can be removed
  select(STUDYID, USUBJID, AVISIT, AVISITN, any_of(c("ATPT", "ATPTN")), CRIT1FL) |>
  inner_join(adsl, by = c("STUDYID", "USUBJID"))


# derive selvisit from data, ordered by AVISITN then ATPTN
selvisit <- if (use_atpt) {
  adpc |>
    arrange(AVISITN, ATPTN) |>
    distinct(AVISIT, ATPT) |>
    mutate(AVISIT_ATPT = paste(AVISIT, ATPT, sep = ", ")) |>
    pull(AVISIT_ATPT)
} else {
  adpc |>
    arrange(AVISITN) |>
    distinct(AVISIT) |>
    pull(AVISIT)
}


adpc <- adpc |>
  # Concatenate or use AVISIT only depending on use_atpt
  mutate(
    AVISIT_ATPT = factor(
      if (use_atpt) paste(AVISIT, ATPT, sep = ", ") else as.character(AVISIT),
      levels = selvisit
    )
  )

################################################################################
# Define layout and build table:
################################################################################

lyt <- basic_table() |>
  split_cols_by(
    var = "colspan_trt",
    split_fun = drop_split_levels
  ) |>
  split_cols_by(
    var = trtvar,
    show_colcounts = TRUE,
    colcount_format = "N=xx",
    split_fun = add_overall_level("Combined", first = FALSE)
  ) |>
  split_rows_by(
    var = "AVISIT_ATPT",
    split_label = "Time Point",
    label_pos = "topleft",
    section_div = " "
  ) |>
  analyze(
    vars = "PKFL",
    show_labels = "hidden",
    afun = a_freq_j,
    extra_args = list(
      val = "Y",
      label = "N",
      .stats = "count_unique"
    )
  ) |>
  analyze(
    vars = "CRIT1FL",
    show_labels = "hidden",
    afun = a_freq_j,
    extra_args = list(
      val = "Y",
      label = paste0("Number of subjects with values below ", lloq_label, ", n (%)"),
      .stats = "count_unique_fraction",
      .indent_mods = 1
    )
  )

result <- build_table(lyt, df = adpc, alt_counts_df = adsl, round_type = "sas")

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

################################################################################
# Convert to tbl file and output table:
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
