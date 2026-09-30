###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tpk01b.r
## R version:                 4.5.2
## junco version:             0.1.3
## Short Description:         [Matrix] [Active Study Agent/Analyte] Concentrations ([units]) Over Time – [SAD/MAD][Part
##                            1]
## Disclaimer:                This script is a direct copy of the corresponding Core Standard output identifier. For
##                            SAD/MAD specific changes, refer to tsfvit02b.r, lsidm05.r, and gsfvit02.r for examples of
##                            STUDYPRT filtering, COHORT handling, treatment column structure modifications, pooled
##                            placebo derivations, combined treatment columns, dose-level updates, and other
##                            output-specific structural differences as applicable.
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adpc
## Output:                    tpk01b.rtf
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

tblid <- "tpk01b"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
popfl <- "PKFL"
trtvar <- "TRT01A"
add_geometric_mean <- TRUE
add_cv <- TRUE
paramcd <- "XAN"


# flag to indicate if ATPT is present in the study
# if TRUE, time point is concatenation of AVISIT and ATPT
# if FALSE, time point is AVISIT only
use_atpt <- TRUE


trt_grps <- c("Xanomeline Low Dose", "Xanomeline High Dose")

################################################################################
# Process data:
################################################################################

adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) |>
  df_na() |>
  filter(.data[[popfl]] == "Y") |>
  filter(.data[[trtvar]] == trt_grps) |>
  select(STUDYID, USUBJID, all_of(trtvar), PKFL) |>
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = trt_grps
    )
  )

adpc <- haven::read_sas(envsetup::read_path(a_in, "adpc.sas7bdat")) |>
  df_na() |>
  filter(PARAMCD == paramcd) |> # If not required this filter can be removed
  select(STUDYID, USUBJID, AVISIT, AVISITN, any_of(c("ATPT", "ATPTN")), AVAL) |>
  inner_join(adsl, by = c("USUBJID", "STUDYID"))

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
    var = trtvar,
    show_colcounts = TRUE,
    colcount_format = "N=xx",
    split_fun = drop_split_levels
  ) |>
  split_rows_by(
    var = "AVISIT_ATPT",
    split_label = "Time Point",
    label_pos = "topleft",
    child_labels = "hidden"
  ) |>
  analyze_vars_in_cols(
    vars = "AVAL",
    .stats = c(
      "n",
      "mean",
      "sd",
      "median",
      "min",
      "max",
      if (add_cv) "cv" else NULL,
      if (add_geometric_mean) "geom_mean" else NULL
    ),
    .labels = c(
      n = "N",
      mean = "Mean",
      sd = "SD",
      median = "Med",
      min = "Min",
      max = "Max",
      cv = "% CV",
      geom_mean = "Geometric Mean"
    ),
    .formats = c(
      n = jjcsformat_xx("xx"),
      mean = format_sigfig_j(3, format = "xx"),
      sd = format_sigfig_j(3, format = "xx"),
      median = format_sigfig_j(3, format = "xx"),
      min = format_sigfig_j(3, format = "xx"),
      max = format_sigfig_j(3, format = "xx"),
      cv = jjcsformat_xx("xx.x"),
      geom_mean = format_sigfig_j(3, format = "xx")
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

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape", combined_rtf = TRUE, nosplitin = list(cols = trtvar))
