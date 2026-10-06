###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsfae01a.r
## R version:                 4.5.0
## Short Description:         Program to create TSFAE01a:
##                            Overall Summary of Subjects With Treatment-emergent
##                            Adverse Events (Study jjcs - clin pharm)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     adaper.sas7bdat, adae.sas7bdat
## Output:                    tsfae01a.rtf
## Remarks:                   This variant includes summary of AEs by max severity.
##                            Template R script version using rtables framework.
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
##                            4. The variable `AESEV` may only contain the
##                               following values (case-insensitive): NA, "",
##                               "Mild", "Moderate", "Severe".
##                            5. In the "Worst severity" row group, for any given
##                               subject, missing severities (NA or "") are
##                               counted only once, and only if that subject has
##                               no other AEs with a non-missing severity
## R-functions:
## R-function Sample Call:
##
## Modification History:
##  Rev #:                    2
##  Reporting Effort:
##  Date:                     March 17, 2026
##  Description:              Update table row-structure.
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

tblid <- "TSFAE01a"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Actual treatment variable (default=TRTA).
trtvar <- "TRTA"

# Control arm name (NULL or one of the levels of the treatment variable).
ctrl_grp <- NULL

# Add Active Study Agent Combined column?
combined_colspan_trt <- FALSE

# Use test data instead of study data.
# This option is available to enhance the testing of this script.
use_test_data <- TRUE

################################################################################
# Process data:
################################################################################

if (use_test_data) {
  set.seed(1)
  adaper <- tibble(
    USUBJID = c(1:15, 1:14, 1:16),
    SAFFL := "Y",
    TRTEMFL := "Y",
    !!trtvar := rep(LETTERS[1:3], times = c(15, 14, 16))
  )
  ae_ids <- sample(nrow(adaper), size = 20, replace = TRUE)
  adae <- tibble(
    adaper[ae_ids, "USUBJID"],
    SAFFL := "Y",
    TRTEMFL := "Y",
    adaper[ae_ids, trtvar],
    AEREL = sample(
      c("RELATED", "NOT RELATED"),
      size = 20,
      replace = TRUE
    ),
    AESEV = sample(
      c("Mild", "Moderate", "Severe", ""),
      size = 20,
      replace = TRUE
    ),
    AESER = sample(c("N", "Y"), size = 20, replace = TRUE),
    AESDTH = sample(c("N", "Y"), size = 20, replace = TRUE),
    AESLIFE = sample(c("N", "Y"), size = 20, replace = TRUE),
    AESHOSP = sample(c("N", "Y"), size = 20, replace = TRUE),
    AESDISAB = sample(c("N", "Y"), size = 20, replace = TRUE),
    AESCONG = sample(c("N", "Y"), size = 20, replace = TRUE),
    AESMIE = sample(c("N", "Y"), size = 20, replace = TRUE),
    AEACN = sample(
      c("DOSE NOT CHANGED", "NOT APPLICABLE", "Drug Withdrawn"),
      size = 20,
      replace = TRUE
    )
  )
} else {
  adaper <- read_sas(read_path(a_in, "adaper.sas7bdat"))
  adae <- read_sas(read_path(a_in, "adae.sas7bdat"))
}

adaper <- adaper |>
  select(USUBJID, SAFFL, all_of(trtvar)) |>
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

adae <- adae |>
  select(
    USUBJID,
    SAFFL,
    TRTEMFL,
    all_of(trtvar),
    AEREL,
    AESER,
    AESDTH,
    AESLIFE,
    AESHOSP,
    AESDISAB,
    AESCONG,
    AESMIE,
    AEACN,
    AESEV
  ) |>
  filter(toupper(SAFFL) == "Y") |>
  filter(toupper(TRTEMFL) == "Y") |>
  filter(!is.na(USUBJID)) |>
  filter(!is.na(.data[[trtvar]]), .data[[trtvar]] != "") |>
  mutate(!!trtvar := factor(.data[[trtvar]])) |>
  mutate(AEREL := factor(AEREL, levels = c("RELATED", "NOT RELATED"))) |>
  mutate(
    across(
      c(AESER, AESDTH, AESLIFE, AESHOSP, AESDISAB, AESCONG, AESMIE),
      ~ factor(toupper(.x), levels = c("N", "Y"))
    )
  ) |>
  mutate(
    REL := factor(ifelse(AEREL == "RELATED", "Y", "N"), levels = c("N", "Y")),
    RELSER := factor(ifelse(REL == "Y" & AESER == "Y", "Y", "N"), levels = c("N", "Y")),
    DST := factor(ifelse(toupper(AEACN) == "DRUG WITHDRAWN", "Y", "N"), levels = c("N", "Y")),
    RELDTH := factor(ifelse(REL == "Y" & AESDTH == "Y", "Y", "N"), levels = c("N", "Y"))
  ) |>
  mutate(
    AESEV := ifelse(is.na(AESEV) | AESEV == "", "Missing", AESEV),
    AESEV := string_to_title(AESEV),
    AESEV := ordered(AESEV, levels = c("Missing", "Mild", "Moderate", "Severe"))
  )

# Variables for the "AEs" group.
# To add custom events (rows) under "AEs", append the variable names and their
# matching labels below.
# Variables must be dichotomous and encoded as 'Y'/'N' (NA allowed).
aes_vars <- c("AESER", "REL", "RELSER", "DST", "RELDTH")
var_labels(adae[, aes_vars]) <- c(
  "SAEs",
  "Related AEs",
  "Related SAEs",
  "AE leading to permanent discontinuation of study treatment",
  "Related AEs leading to death"
)

# Variables for SEA classification.
sae_classification_vars <- c("AESDTH", "AESLIFE", "AESHOSP", "AESDISAB", "AESCONG", "AESMIE")
var_labels(adae[, sae_classification_vars]) <- c(
  "Death",
  "Life-threatening",
  "Requires or prolongs hospitalization",
  "Persistent or significant disability/incapacity",
  "Congenital anomaly or birth defect",
  "Other medically important event"
)

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

aesevall_spf <- make_combo_splitfun(
  nm = "AESEV_ALL",
  label = "Worst severity",
  levels = NULL,
)

aeserall_spf <- make_combo_splitfun(
  nm = "AESER_ALL",
  label = "SAE classification~[super a]",
  levels = "Y",
)

acfun_extra_args <- list(
  .stats = "count_unique_fraction",
  denom = "n_altdf"
)

ae_label_map <- data.frame(
  var = aes_vars,
  value = "Y",
  label = var_labels(adae[, aes_vars])
)

sae_label_map <- data.frame(
  var = sae_classification_vars,
  value = "Y",
  label = var_labels(adae[, sae_classification_vars])
)

lyt <- basic_table(
  show_colcounts = TRUE,
  colcount_format = "N=xx",
  top_level_section_div = " "
) |>
  append_topleft(c(" ", " ", "Event, n (%)")) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar, split_fun = split_combined) |>
  add_overall_col("Total") |>
  analyze(
    "TRTEMFL",
    afun = a_freq_j,
    extra_args = c(val = "Y", label = "AEs", acfun_extra_args),
    show_labels = "hidden"
  ) |>
  analyze(
    unique(ae_label_map$var),
    a_freq_j,
    extra_args = c(val = "Y", label_map = list(ae_label_map), acfun_extra_args),
    show_labels = "hidden",
    indent_mod = 1L
  ) |>
  split_rows_by("AESEV", split_fun = aesevall_spf) |>
  analyze("AESEV", afun = a_maxlev) |>
  split_rows_by("AESER", split_fun = aeserall_spf) |>
  analyze(
    unique(sae_label_map$var),
    a_freq_j,
    extra_args = c(val = "Y", label_map = list(sae_label_map), acfun_extra_args),
    show_labels = "hidden"
  )

result <- build_table(lyt, adae, alt_counts_df = adaper)
result

################################################################################
# Post-Processing:
# - Adjust Combined (if displayed) and Total columns Ns (When sequence of
#   treatments is used, one subject receives more than on treatment. Hence,
#   there are many rows for one unique subjects, while N should represent only
#   unique subjects).
# - Prune "Missing AEs" category with all zeros.
################################################################################

if (combined_colspan_trt) {
  adaper_no_ctrl <- if (is.null(ctrl_grp)) {
    adaper
  } else {
    adaper[adaper[[trtvar]] != ctrl_grp, ]
  }
  n_combined <- length(unique(adaper_no_ctrl$USUBJID))
  colpath_combined <- c("colspan_trt", "Active Study Agent", trtvar, "Combined")
  facet_colcount(result, colpath_combined) <- n_combined
}

n_total <- length(unique(adaper$USUBJID))
colpath_total <- c("Total", "Total")
facet_colcount(result, colpath_total) <- n_total

result <- safe_prune_table(
  result,
  prune_func = count_pruner(cols = trtvar, cat_include = "Missing")
)

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)
result

################################################################################
# Convert to tbl file and output table:
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
