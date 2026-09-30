###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsfdth01.r
## R version:                 4.2.1
## junco version:             0.1.3
## Short Description:         Program to create tsfdth01: Table of deaths
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl
## Output:                    tsfdth01.rtf
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

# Section --- sentinel row removal utility ---------------------------------
# Recursively removes DataRows matching a label from the table tree.
# Used to strip placeholder rows ("NA_VALUES") inserted for categories
# that have no preferred term underneath.
remove_datarows_by_label <- function(tt, label) {
  if (inherits(tt, "TableTree") || inherits(tt, "ElementaryTable")) {
    kids <- tree_children(tt)
    kids <- Filter(
      function(k) {
        !(inherits(k, "DataRow") && obj_label(k) == label)
      },
      kids
    )
    tree_children(tt) <- lapply(kids, remove_datarows_by_label, label = label)
  }
  tt
}

################################################################################
# Define script level parameters:
################################################################################

tblid <- "TSFDTH01"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

trtvar <- "TRT01A"
popfl <- "SAFFL"

days <- 30

# Either Combined column or highest JnJ Dose level based on study
sort_colpath <- "Xanomeline High Dose"

combined_colspan_trt <- TRUE
risk_diff <- TRUE
rr_method <- "wald"
ctrl_grp <- "Placebo"
include_dth_first_dose <- TRUE

if (combined_colspan_trt) {
  # Set up levels and label for the required combined columns
  add_combo <- add_combo_facet(
    "Combined",
    label = "Combined",
    levels = c("Xanomeline High Dose", "Xanomeline Low Dose")
  )

  # choose if any facets need to be removed - e.g remove the combined column for placebo
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

adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) %>%
  df_na() %>%
  filter(!!rlang::sym(popfl) == "Y") %>%
  select(
    STUDYID,
    USUBJID,
    all_of(trtvar),
    all_of(popfl),
    DTHFL,
    DTHTRTFL,
    DTHAFTFL,
    DTHB60FL,
    DDPCDTHC,
    DTHCAUS
  ) %>%
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c(
        "Xanomeline Low Dose",
        "Xanomeline High Dose",
        "Placebo"
      )
    )
  ) %>%
  create_colspan_var(
    non_active_grp = "Placebo",
    non_active_grp_span_lbl = " ",
    active_grp_span_lbl = "Active Study Agent",
    colspan_var = "colspan_trt",
    trt_var = trtvar
  ) %>%
  mutate(
    rrisk_header = "Risk Difference (%) (95% CI)",
    rrisk_label = paste(!!rlang::sym(trtvar), "vs Placebo")
  )

# Sentinel "NA_VALUES" inserted where DDPCDTHC exists but DTHCAUS is missing;
# these rows are removed post-build so the parent category row remains.
adsl <- adsl %>%
  mutate(
    DTHCAUS = as.character(DTHCAUS),
    DTHCAUS = case_when(
      !is.na(DTHCAUS) ~ stringr::str_to_sentence(tolower(DTHCAUS)),
      !is.na(DDPCDTHC) ~ "NA_VALUES",
      TRUE ~ NA_character_
    )
  )

# Set DDPCDTHC factor levels in spec-defined order
adsl$DDPCDTHC <- factor(
  stringr::str_to_sentence(as.character(adsl$DDPCDTHC)),
  levels = intersect(
    c("Adverse event", "Disease progression of trial indication", "Treatment failure/relapse", "Other"),
    unique(stringr::str_to_sentence(as.character(adsl$DDPCDTHC[!is.na(adsl$DDPCDTHC)])))
  )
)
################################################################################
# Define layout and build table:
################################################################################
colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = ctrl_grp,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

ref_path <- c("colspan_trt", " ", trtvar, ctrl_grp)

extra_args_rr <- list(
  method = rr_method,
  .stats = c("count_unique_fraction"),
  ref_path = ref_path,
  riskdiff = TRUE,
  denom = "n_altdf"
)

has_deaths <- nrow(filter(adsl, DTHFL == "Y")) > 0

lyt <- basic_table(
  top_level_section_div = " ",
  show_colcounts = TRUE,
  colcount_format = "N=xx"
) %>%
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  )

if (combined_colspan_trt) {
  lyt <- lyt %>%
    split_cols_by(trtvar, split_fun = mysplit)
} else {
  lyt <- lyt %>%
    split_cols_by(trtvar)
}

if (risk_diff) {
  lyt <- lyt %>%
    split_cols_by("rrisk_header", nested = FALSE) %>%
    split_cols_by(
      trtvar,
      labels_var = "rrisk_label",
      split_fun = remove_split_levels(ctrl_grp)
    )
}

# Section --- death block definitions --------------------
death_blocks <- list(
  list(var = "DTHFL", label = "Total deaths"),
  list(var = "DTHTRTFL", label = paste0("Deaths within ", days, " days of last dose")),
  list(var = "DTHAFTFL", label = paste0("Deaths beyond ", days, " days of last dose"))
)

if (include_dth_first_dose) {
  death_blocks <- death_blocks |>
    append(
      list(
        list(var = "DTHB60FL", label = "Deaths within 60 days of first dose")
      )
    )
}

# Section --- build row layout via loop --------------------
if (has_deaths) {
  for (i in seq_along(death_blocks)) {
    block <- death_blocks[[i]]

    if (i == 1) {
      lyt <- lyt %>%
        split_rows_by(
          block$var,
          split_fun = keep_split_levels("Y"),
          split_label = "Deaths",
          label_pos = "topleft",
          section_div = " "
        )
    } else {
      lyt <- lyt %>%
        split_rows_by(
          block$var,
          split_fun = keep_split_levels("Y"),
          section_div = " "
        )
    }

    lyt <- lyt %>%
      summarize_row_groups(
        block$var,
        cfun = a_freq_j,
        extra_args = append(extra_args_rr, list(label = block$label))
      ) %>%
      split_rows_by(
        "DDPCDTHC",
        split_fun = drop_split_levels
      ) %>%
      summarize_row_groups(
        "DDPCDTHC",
        cfun = a_freq_j,
        extra_args = extra_args_rr
      ) %>%
      analyze(
        "DTHCAUS",
        afun = a_freq_j,
        extra_args = extra_args_rr,
        indent_mod = 0,
        show_labels = "hidden"
      )
  }
} else {
  lyt <- lyt %>%
    analyze(
      "DTHFL",
      a_freq_j,
      show_labels = "hidden",
      extra_args = append(extra_args_rr, list(label = "Total deaths"))
    )
}

lyt <- lyt %>%
  append_topleft("  Cause of Death, n (%)")

result <- build_table(lyt, adsl, alt_counts_df = adsl, round_type = "sas")

# Section --- post-process table -------------------------------------------
result <- remove_datarows_by_label(result, "NA_VALUES")
result <- suppressWarnings(remove_col_count(result))

# Section --- sort preferred terms by descending incidence -----------------
# DDPCDTHC categories: fixed spec order (factor levels).
# DTHCAUS within each category: sorted by decreasing incidence on sort_colpath.
if (has_deaths) {
  for (block in death_blocks) {
    for (cat in levels(adsl$DDPCDTHC)) {
      tryCatch(
        result <- suppressWarnings(
          sort_at_path(
            result,
            c(block$var, "Y", "DDPCDTHC", cat, "DTHCAUS"),
            scorefun = jj_complex_scorefun(colpath = sort_colpath)
          )
        ),
        error = function(e) NULL
      )
    }
  }
}

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

################################################################################
# Convert to tbl file and output table
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
