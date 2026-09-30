###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsfvit01p.r
## R version:                 4.5.2
## junco version:             0.1.3
## Short Description:         Program to create tsfvit01p: Mean and Mean Change
##                            From Baseline for Vital Sign Data Over Time.
## Author:                    C&SP Methodology
## Date:                      2026-09-302026-09-30
## Input:                     adaper, advs
## Output:                    tsfvit01p.rtf
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
library(rtables)
library(junco)
library(haven)


################################################################################
# Define script level parameters:
################################################################################

tblid <- "TSFVIT01P"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Safety population flag (default=SAFFL).
popfl <- "SAFFL"

# Actual treatment variable (default=TRTA).
trtvar <- "TRTA"
# Actual treatment levels (required to order the columns).
trtlev <- c("A", "B", "C", "D")

# Control arm name (NULL or one of the levels of the treatment variable).
ctrl_grp <- NULL

# Add Active Study Agent Combined column?
combined_colspan_trt <- FALSE

# Vital sign parameters of interest.
paramcd <- c("SYSBPU", "DIABPU")

################################################################################
# Process data:
################################################################################

adaper <- read_sas(read_path(a_in, "adaper.sas7bdat")) |>
  select(USUBJID, all_of(c(popfl, trtvar))) |>
  filter(
    toupper(.data[[popfl]]) == "Y",
    !is.na(.data[[trtvar]])
  ) |>
  mutate(
    !!trtvar := ordered(.data[[trtvar]], levels = trtlev),
  ) |>
  create_colspan_var(
    non_active_grp = ctrl_grp,
    non_active_grp_span_lbl = " ",
    active_grp_span_lbl = "Active Study Agent",
    colspan_var = "colspan_trt",
    trt_var = trtvar
  )

advs <- read_sas(read_path(a_in, "advs.sas7bdat")) |>
  select(
    USUBJID,
    all_of(c(popfl, trtvar)),
    AVISITN,
    AVISIT,
    ATPT,
    PARAMCD,
    PARAM,
    AVAL,
    BASE,
    CHG,
    starts_with("ANL"),
    ABLFL,
    APOBLFL
  ) |>
  filter(
    !is.na(USUBJID),
    (.data[[popfl]] == "Y"),
    !is.na(.data[[trtvar]]),
    PARAMCD %in% paramcd
  ) |>
  mutate(
    !!trtvar := ordered(.data[[trtvar]], levels = trtlev),
    ABLFL := factor(ifelse(ABLFL == "Y" & !is.na(ABLFL), "Y", "N")),
    APOBLFL := factor(ifelse(APOBLFL == "Y" & !is.na(APOBLFL), "Y", "N")),
    PARAMCD := ordered(PARAMCD, levels = paramcd),
    PARAM := factor(PARAM),
    AVISITTPT = ordered(
      case_when(
        ABLFL == "Y" ~ "Baseline",
        toupper(AVISIT) == "DAY 1" & toupper(ATPT) == "2H" ~ "2 hours postdose",
        toupper(AVISIT) == "DAY 1" & toupper(ATPT) == "4H" ~ "4 hours postdose",
        toupper(AVISIT) == "DAY 1" & toupper(ATPT) == "8H" ~ "8 hours postdose",
        toupper(AVISIT) == "DAY 2" & toupper(ATPT) == "24H" ~ "24 hours postdose",
        toupper(AVISIT) == "DAY 3" & toupper(ATPT) == "48H" ~ "48 hours postdose",
        toupper(AVISIT) == "EOS" ~ "EOS",
        TRUE ~ NA_character_
      ),
      levels = c("Baseline", paste(c("2", "4", "8", "24", "48"), "hours postdose"), "EOS")
    )
  )

advs <- inner_join(adaper, advs, by = c("USUBJID", popfl, trtvar))

# Safety checks on the data.

advs_checks_unique <- advs |>
  filter(ANL02FL == "Y" & (ABLFL == "Y" | APOBLFL == "Y")) |>
  group_by(USUBJID, PARAMCD, AVISIT) |>
  mutate(n_recsub = n()) |>
  filter(n_recsub > 1)

advs_checked <- advs |>
  filter(
    toupper(.data[[popfl]]) == "Y",
    toupper(ANL01FL) == "Y",
    toupper(ANL02FL) == "Y"
  ) |>
  mutate(
    flag = if_else(is.na(AVAL), 1L, 0L),
    flag1 = if_else(toupper(APOBLFL) == "Y" & is.na(BASE), 1L, 0L)
  )

if (nrow(filter(advs_checked, flag == 1)) > 0) {
  stop("Pay attention: There are tests with missing AVAL: ")
}

if (nrow(filter(advs_checked, flag1 == 1)) > 0) {
  stop("Pay attention: There are postbaseline tests without baseline: ")
}

# Keep only records with non-missing AVAL
advs <- filter(advs, !is.na(AVAL))

################################################################################
# Define layout and build table:
################################################################################

colspan_trt_map <- if (!is.null(ctrl_grp)) {
  create_colspan_map(
    advs,
    non_active_grp = ctrl_grp,
    non_active_grp_span_lbl = " ",
    active_grp_span_lbl = "Active Study Agent",
    colspan_var = "colspan_trt",
    trt_var = trtvar
  )
} else {
  tibble(
    colspan_trt = "Active Study Agent",
    !!trtvar := levels(advs[[trtvar]])
  )
}

split_combined <- if (combined_colspan_trt) {
  # Set up levels and label for the required combined columns.
  add_combo <- add_combo_facet(
    "Combined",
    label = "Combined",
    levels = setdiff(levels(advs[[trtvar]]), ctrl_grp)
  )

  # Choose if any facets need to be removed,
  # e.g remove the combined column for placebo.
  rm_combo_from_placebo <- cond_rm_facets(
    facets = "Combined",
    ancestor_pos = NA,
    value = " ",
    split = "colspan_trt"
  )

  make_split_fun(post = list(add_combo, rm_combo_from_placebo))
} else {
  NULL
}

# Common stats for time point / change from baseline to time point.
stats <- c("n", "mean_sd", "median_range")
formats <- junco_get_formats_from_stats(stats)
labels <- junco_get_labels_from_stats(stats, labels_in = c(n = "N"))
indent_mods <- c(n = 1L, mean_sd = 2L, median_range = 2L)

# Dynamic label for the "Change from baseline ..." sections.
a_chg_label <- function(x, .spl_context) {
  last_split <- length(.spl_context$split)
  label <- paste("Change from baseline to", .spl_context$value[last_split])
  rtables::rcell(NULL, label = label)
}

lyt <- basic_table(
  show_colcounts = TRUE,
  colcount_format = "N=xx",
  top_level_section_div = " "
) |>
  append_topleft("Parameter") |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar, split_fun = split_combined) |>
  add_overall_col("Total") |>
  split_rows_by(
    "PARAMCD",
    split_fun = keep_split_levels(paramcd),
    labels_var = "PARAM",
    child_labels = "visible"
  ) |>
  summarize_row_groups(
    "AVAL",
    cfun = c_summary_subset_label,
    extra_args = list(
      filter_expr = expression(ABLFL == "Y"),
      .stats = stats,
      .formats = formats,
      .labels = labels,
      .indent_mods = indent_mods,
      label = "Baseline"
    )
  ) |>
  split_rows_by(
    "AVISITTPT",
    split_fun = remove_split_levels("Baseline"),
    indent_mod = -1,
    section_div = " ",
  ) |>
  analyze(
    "AVAL",
    afun = tern::a_summary,
    extra_args = list(
      .stats = stats,
      .formats = formats,
      .labels = labels,
      .indent_mods = indent_mods - 1,
      na_rm = TRUE
    ),
    show_labels = "hidden",
    section_div = " "
  ) |>
  analyze(
    "CHG",
    afun = a_chg_label,
    table_names = "chg_bl_label",
    show_labels = "hidden",
    indent_mod = -1L
  ) |>
  analyze(
    "CHG",
    afun = tern::a_summary,
    extra_args = list(
      .stats = "n",
      .formats = formats,
      .labels = labels,
      na_rm = TRUE # Subjects with non-missing values at both baseline and the postbaseline time point.
    ),
    show_labels = "hidden"
  ) |>
  analyze(
    "BASE",
    afun = a_summary_subset,
    extra_args = list(
      filter_expr = expression(!is.na(CHG)),
      na_rm = TRUE,
      .stats = "mean_sd",
      .formats = formats,
      .labels = c(mean_sd = "Baseline Mean (SD)")
    ),
    table_names = "chg_bl_mean",
    show_labels = "hidden",
    indent_mod = 1L
  ) |>
  analyze(
    "CHG",
    afun = tern::a_summary,
    extra_args = list(
      .stats = c("mean_sd", "median_range"),
      .formats = formats,
      .labels = labels,
      na_rm = TRUE
    ),
    table_names = "chg_mm",
    show_labels = "hidden",
    indent_mod = 1L
  )

result <- build_table(lyt, advs, alt_counts_df = adaper, round_type = "sas")

# Add an empty line after each "Screening" section.
# Use `section_div_at_path()` as a workaround for the missing
# `section_div` argument in `summarize_row_groups()`.
# (Feature request: https://github.com/insightsengineering/rtables/issues/1083)
# The path below is determined from `rtables::row_paths_summary(result)`.
section_div_at_path(result, c("PARAMCD", "*", "@content", last(stats))) <- " "

################################################################################
# Post-Processing:
# - Adjust Combined (if displayed)  and Total columns Ns (When sequence of
#   treatments is used, one subject receives more than on treatment. Hence,
#   there are many rows for one unique subjects, while N should represent only
#   unique subjects).
################################################################################

if (combined_colspan_trt) {
  adaper_no_ctrl <- if (is.null(ctrl_grp)) {
    adaper
  } else {
    adaper[adaper[[trtvar]] != ctrl_grp, ]
  }
  N_combined <- length(unique(adaper_no_ctrl$USUBJID))
  colpath_combined <- c("colspan_trt", "Active Study Agent", trtvar, "Combined")
  facet_colcount(result, colpath_combined) <- N_combined
}

colpath_total <- c("Total", "Total")
N_total <- length(unique(adaper$USUBJID))
facet_colcount(result, colpath_total) <- N_total

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

################################################################################
# Convert to tbl file and output table:
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
