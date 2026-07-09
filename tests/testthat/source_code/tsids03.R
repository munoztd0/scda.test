################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsids03.r
## R version:                 4.5.3
## junco Version:             0.1.3
## Short Description:         Program to create tsids03: Analysis Sets
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl
## Output:                    tsids03.rtf
## Remarks:
##
## Modification History:
##  Rev #:
##  Modified By:
##  Reporting Effort:
##  Date:
##  Description:
################################################################################

###############################################################################
#Prep environment
###############################################################################
library(envsetup)
library(tern)
library(dplyr)
library(rtables)
library(junco)

###############################################################################
#Define script level parameters
###############################################################################

tblid <- "TSIDS03"
fileid <- write_path(opath, tblid)
popfl <- "RANDFL"
trtvar <- "TRT01P"
combined_colspan_trt <- TRUE
total_col <- TRUE
ctrl_grp <- "Placebo"
tab_titles <- list(title = "Dummy Title",
                     subtitles = NULL,
                     main_footer = "Dummy Note: On-treatment is defined as ~{optional treatment-emergent}")
analysis_set_vars <- c("FASFL", "SAFFL", "PPROTFL")
analysis_set_vars_lbls <- c(FASFL = "Full", SAFFL = "Safety", PPROTFL = "Per protocol")
anl_vars_label_map <- data.frame(
  var = analysis_set_vars,
  value = "Y",
  label = c("Full", "Safety", "Per protocol")
)
if (combined_colspan_trt == TRUE) {
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

###############################################################################
#Process data
###############################################################################

adsl <- adsl_jnj %>%
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
  filter(!!rlang::sym(popfl) == "Y")
adsl$colspan_trt <- factor(
  ifelse(adsl[[trtvar]] == "Placebo", " ", "Active Study Agent"),
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

###############################################################################
#Define layout and build table
###############################################################################

lyt <- basic_table(
  show_colcounts = TRUE,
  colcount_format = "N=xx",
  top_level_section_div = " "
) %>%
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  )
if (combined_colspan_trt == TRUE) {
  lyt <- lyt %>%
    split_cols_by(trtvar, split_fun = mysplit)
} else {
  lyt <- lyt %>%
    split_cols_by(trtvar)
}
if (total_col == TRUE) {
  lyt <- lyt %>%
    add_overall_col(label = "Total")
}
lyt <- lyt %>%
  analyze(
    vars = analysis_set_vars,
    afun = a_freq_j,
    extra_args = append(
      list(
        .stats = c("count_unique_fraction")
      ),
      list(
        val = "Y",
        label_map = anl_vars_label_map,
        riskdiff = FALSE,
        extrablankline = FALSE,
        NULL
      )
    ),
    show_labels = "hidden",
    na_str = " "
  ) %>%
  append_topleft("Analysis Set")

result <- build_table(lyt, df = adsl, round_type = "sas")

###############################################################################
#Post-processing
###############################################################################

if (nrow(adsl) == 0) {
  #Post-processing step to remove table rows with all 0 or NA values
  result <- safe_prune_table(result, prune_func = prune_empty_level)
}
###############################################################################
#Retrieve titles and footnotes
###############################################################################

result <- set_titles(result, tab_titles)

###############################################################################
#Convert to tbl file and output table
###############################################################################


colwidth <- c(20, 23, 23, 25, 23, 25)

tt_to_tlgrtf(colwidths = colwidth, result, file = fileid)
