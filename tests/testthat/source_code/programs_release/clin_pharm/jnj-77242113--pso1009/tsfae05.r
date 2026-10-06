###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsfae05.r
## R version:                 4.5.0
## Short Description:         Program to create TSFAE05:
#                             Subjects With Treatment-emergent Adverse Events
#                             (Excluding Serious Adverse Events) With Frequency
#                             of ≥[x]% in Any Treatment Group by System Organ
#                             Class and Preferred Term; Safety Analysis Set
#                             (Study jjcs - clin pharm).
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     adaper.sas7bdat, adae.sas7bdat
## Output:                    tsfae05.rtf
## Remarks:                   Template R script version using rtables framework.
##
## Assumptions:               1. The script handles the sequence of treatments.
##                               The actual treatment variable (specified via
##                               the `trtvar` parameter in this R script) must
##                               be present in both the ADAPER and ADAE datasets.
##                            2. All statistics are computed based on the actual
##                               treatment variable.
##                            3. Rows with a missing value for the treatment
##                               variable are excluded from all calculations,
##                               including those for the Total column.
## R-functions:
## R-function Sample Call:
##
## Modification History:
##  Rev #:                    2
##  Reporting Effort:
##  Date:                     March 19, 2026
##  Description:              Update table name from TSFAE07 to TSFAE05.
################################################################################

################################################################################
# Prep environment:
################################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(dplyr)
library(rtables)
library(tern)
library(junco)
library(haven)


################################################################################
# Define script level parameters:
################################################################################

tblid <- "TSFAE05"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Actual treatment variable (default=TRTA).
trtvar <- "TRTA"

# Control arm name (NULL or one of the levels of the treatment variable).
ctrl_grp <- NULL

# Add Active Study Agent Combined column?
combined_colspan_trt <- FALSE

# Minimum frequency of adverse events (excluding serious adverse events) in any
# active treatment group, to be displayed in the table.
# The 'Combined' (if present) and 'Total' columns are ignored during this minimum
# frequency filtering.
min_freq <- 0.05

################################################################################
# Process data:
################################################################################

adaper <- read_sas(read_path(a_in, "adaper.sas7bdat")) |>
  select(USUBJID, SAFFL, TRTSEQA, all_of(trtvar)) |>
  filter(toupper(SAFFL) == "Y") |>
  filter(!is.na(.data[[trtvar]]), .data[[trtvar]] != "") |>
  mutate(!!trtvar := factor(.data[[trtvar]])) |>
  create_colspan_var(
    non_active_grp = ctrl_grp,
    non_active_grp_span_lbl = " ",
    active_grp_span_lbl = "Active Study Agent",
    colspan_var = "colspan_trt",
    trt_var = trtvar
  )

adae <- read_sas(read_path(a_in, "adae.sas7bdat")) |>
  select(
    USUBJID,
    SAFFL,
    TRTEMFL,
    all_of(trtvar),
    AESER,
    AEBODSYS,
    AEDECOD,
    AETERM
  ) |>
  filter(toupper(SAFFL) == "Y") |>
  filter(toupper(TRTEMFL) == "Y") |>
  filter(toupper(AESER) == "N") |>
  filter(!is.na(USUBJID)) |>
  filter(!is.na(.data[[trtvar]]), .data[[trtvar]] != "") |>
  mutate(across(c(AEBODSYS, AEDECOD), string_to_title)) |>
  # Handling missing coded term.
  mutate(
    AEBODSYS = ifelse(is.na(AEBODSYS) & !is.na(AETERM), "Uncoded", AEBODSYS),
    AEDECOD = ifelse(
      is.na(AEDECOD) & !is.na(AETERM),
      paste("Uncoded -", trimws(AETERM)),
      AEDECOD
    )
  ) |>
  mutate(across(c(!!trtvar, AEBODSYS, AEDECOD), factor))

adae <- inner_join(adaper, adae, by = c("USUBJID", "SAFFL", trtvar))

################################################################################
# Define layout and build table:
################################################################################

colspan_trt_map <- if (!is.null(ctrl_grp)) {
  create_colspan_map(
    adae,
    non_active_grp = ctrl_grp,
    non_active_grp_span_lbl = " ",
    active_grp_span_lbl = "Active Study Agent",
    colspan_var = "colspan_trt",
    trt_var = trtvar
  )
} else {
  tibble(
    colspan_trt = "Active Study Agent",
    !!trtvar := levels(adae[[trtvar]])
  )
}

split_combined <- if (combined_colspan_trt) {
  # Set up levels and label for the required combined columns.
  add_combo <- add_combo_facet(
    "Combined",
    label = "Combined",
    levels = setdiff(levels(adae[[trtvar]]), ctrl_grp)
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

acfun_extra_args <- list(
  .stats = "count_unique_fraction",
  denom = "n_altdf"
)

lyt <- basic_table(
  show_colcounts = TRUE,
  colcount_format = "N=xx",
  top_level_section_div = " "
) |>
  append_topleft(c("System Organ Class", "  Preferred Term, n (%)")) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar, split_fun = split_combined) |>
  add_overall_col("Total") |>
  analyze(
    vars = "TRTEMFL",
    afun = a_freq_j,
    extra_args = c(label = "Subjects with >=1 AE", val = "Y", acfun_extra_args),
    show_labels = "hidden",
  ) |>
  split_rows_by(
    var = "AEBODSYS",
    split_label = "",
    split_fun = trim_levels_in_group("AEDECOD"),
    nested = FALSE,
    label_pos = "topleft",
    section_div = c(" ")
  ) |>
  summarize_row_groups("AEBODSYS", cfun = a_freq_j, extra_args = acfun_extra_args) |>
  analyze("AEDECOD", afun = a_freq_j, extra_args = acfun_extra_args)

result <- build_table(lyt, adae, alt_counts_df = adaper)
result

################################################################################
# Post-Processing:
# - Adjust Total column Ns (when sequence of treatments is used, one subject
#   receives more than on treatment. Hence, there are many rows for one unique
#   subjects, while N should represent only unique subjects).
# - Sort by descending count on a chosen column: Combined (if displayed) or
#   Total (if Combined not displayed).
#   See function documentation for jj_complex_scorefun should your require a
#   different sorting behavior.
################################################################################

colpath_total <- c("Total", "Total")
sort_colpath <- if (combined_colspan_trt) {
  adaper_no_ctrl <- if (is.null(ctrl_grp)) {
    adaper
  } else {
    adaper[adaper[[trtvar]] != ctrl_grp, ]
  }
  n_combined <- length(unique(adaper_no_ctrl$USUBJID))
  colpath_combined <- c("colspan_trt", "Active Study Agent", trtvar, "Combined")
  facet_colcount(result, colpath_combined) <- n_combined
  colpath_combined
} else {
  colpath_total
}

n_total <- length(unique(adaper$USUBJID))
facet_colcount(result, colpath_total) <- n_total

if (nrow(adae) != 0) {
  result <- sort_at_path(
    result,
    path = c("AEBODSYS"),
    scorefun = jj_complex_scorefun(colpath = sort_colpath)
  ) |>
    sort_at_path(
      path = c("AEBODSYS", "*", "AEDECOD"),
      scorefun = jj_complex_scorefun(colpath = sort_colpath)
    )

  ##############################################################################
  # Prune table to only keep those that meet x% criteria for any active
  # treatment group (column)
  ##############################################################################
  acttrt <- colspan_trt_map[
    colspan_trt_map$colspan_trt == "Active Study Agent",
  ]
  more_than_x_percent <- has_fraction_in_any_col(
    atleast = min_freq,
    col_names = paste(acttrt$colspan_trt, acttrt[[trtvar]], sep = ".")
  )
  result <- safe_prune_table(result, keep_rows(more_than_x_percent))
}

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)
result

################################################################################
# Convert to tbl file and output table:
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
