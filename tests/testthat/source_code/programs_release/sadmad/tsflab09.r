###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsflab09.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Subjects in Each Quadrant for Cholestatic Drug-induced Liver Injury Screening Plot –
##                            [SAD/MAD][Part 1]
## Disclaimer:                This script is a direct copy of the corresponding Core Standard output identifier. For
##                            SAD/MAD specific changes, refer to tsfvit02b.r, lsidm05.r, and gsfvit02.r for examples of
##                            STUDYPRT filtering, COHORT handling, treatment column structure modifications, pooled
##                            placebo derivations, combined treatment columns, dose-level updates, and other
##                            output-specific structural differences as applicable.
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, addili
## Output:                    tsflab09.rtf
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
# Prep Environment
################################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)
library(dplyr)
library(tidyr)
library(lubridate)
library(rtables)
library(junco)

################################################################################
# Define script level parameters:
################################################################################

tblid <- "tsflab09"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Population flag variable (default=SAFFL).
popfl <- "SAFFL"
# Actual treatment variable (default=TRT01A).
trtvar <- "TRT01A"

ctrl_grp <- "Placebo"

combined_colspan_trt <- TRUE

if (combined_colspan_trt == TRUE) {
  add_combo <- add_combo_facet(
    "Combined",
    label = "Combined",
    levels = c("Xanomeline High Dose", "Xanomeline Low Dose")
  )

  rm_combo_from_placebo <- cond_rm_facets(
    facets = "Combined",
    ancestor_pos = NA,
    value = " ",
    split = "colspan_trt"
  )

  mysplit <- make_split_fun(post = list(add_combo, rm_combo_from_placebo))
}


################################################################################
# Process Data:
################################################################################

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  df_na() |>
  filter(.data[[popfl]] == "Y") |>
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c("Xanomeline Low Dose", "Xanomeline High Dose", "Placebo")
    )
  ) |>
  select(STUDYID, USUBJID, all_of(trtvar))

adsl$colspan_trt <- factor(
  ifelse(adsl[[trtvar]] == ctrl_grp, " ", "Active Study Agent"),
  levels = c("Active Study Agent", " ")
)

colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = ctrl_grp,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

addili <- haven::read_sas(read_path(a_in, "addili.sas7bdat")) |>
  df_na() |>
  filter(
    .data[[popfl]] == "Y",
    PARAMCD == "CDILI",
    !is.na(.data[[trtvar]])
  ) |>
  mutate(
    AVALCAT = reorder(factor(paste0(AVALC, " (", AVALCAT2, ")")), AVALCA2N)
  ) |>
  left_join(adsl |> select(USUBJID, colspan_trt), by = "USUBJID")


# Labels
quad_levels <- levels(addili$AVALCAT)
cdili_label <- paste0(unique(as.character(na.omit(addili$CRIT1))), "~[super a]")


################################################################################
# Define layout and build table:
################################################################################

extra_args <- list(
  denom = "n_parentdf",
  denom_by = "PARAMCD",
  .stats = "count_unique_fraction"
)

lyt <- basic_table(show_colcounts = TRUE, colcount_format = "N=xx") |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  )

if (combined_colspan_trt == TRUE) {
  lyt <- lyt |> split_cols_by(trtvar, split_fun = mysplit)
} else {
  lyt <- lyt |> split_cols_by(trtvar)
}

lyt <- lyt |>
  append_topleft(c(" ", "Quadrant, n (%)")) |>
  split_rows_by(
    "PARAMCD",
    split_fun = keep_split_levels("CDILI"),
    nested = FALSE,
    child_labels = "hidden"
  ) |>
  analyze(
    "PARAMCD",
    afun = a_freq_j,
    extra_args = list(
      val = "CDILI",
      .stats = "count",
      label = "Subjects with on-treatment ALP and TBILI values",
      extrablankline = TRUE
    ),
    show_labels = "hidden"
  ) |>
  split_rows_by(
    "AVALCAT",
    parent_name = "q1",
    split_fun = keep_split_levels(quad_levels[1]),
    nested = FALSE,
    child_labels = "hidden"
  ) |>
  analyze(
    "PARAMCD",
    afun = a_freq_j,
    extra_args = append(
      extra_args,
      list(val = "CDILI", label = quad_levels[1])
    ),
    show_labels = "hidden"
  ) |>
  analyze(
    "CRIT1FL",
    afun = a_freq_j,
    extra_args = append(
      extra_args,
      list(
        val = "Y",
        label = cdili_label
      )
    ),
    show_labels = "hidden",
    indent_mod = 1L
  ) |>
  split_rows_by(
    "AVALCAT",
    parent_name = "q2-q4",
    split_fun = keep_split_levels(quad_levels[2:4]),
    nested = FALSE,
    child_labels = "hidden"
  ) |>
  analyze(
    "AVALCAT",
    afun = a_freq_j,
    extra_args = append(extra_args, list(drop_levels = TRUE)),
    show_labels = "hidden"
  )


result <- build_table(
  lyt,
  addili,
  alt_counts_df = adsl,
  round_type = "sas"
)


################################################################################
# Titles and Footnotes
################################################################################

result <- set_titles(result, tab_titles)

################################################################################
# Output
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, label_width_ins = 2.4, orientation = "landscape")
