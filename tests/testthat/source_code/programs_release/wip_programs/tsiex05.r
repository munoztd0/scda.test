###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsiex05.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Relative Dose Intensity
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADEXSUM.
## Output:                    TSIEX05.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:                     1
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
## Description:
################################################################################

################################################################################
# Prep Environment
################################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)

library(tern)
library(stringr)
library(tidytlg)
library(grid)

################################################################################
# - Define output ID and file location
# - Choose whether or not you want to present a combined active treatment column (default=Y)
# - Define how to create combined treatment columns (if required)
################################################################################

tblid <- "TSIEX05"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()


combined_active_trt <- "Y"

if (combined_active_trt == "Y") {
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
    split = "active_trt"
  )

  mysplit <- make_split_fun(post = list(add_combo, rm_combo_from_placebo))
}

################################################################################
# Process Data:
################################################################################

# Read in required data
sas_vs_rds_check("adsl", a_in)
adsl <- readRDS(read_path(a_in, "adsl.rds")) %>%
  filter(SAFFL == "Y" & TRT01A != "Placebo") %>%
  select(STUDYID, USUBJID, TRT01A, SAFFL)

# If AVISIT is not present in ADEXSUM than create 'Overall' as the Visit
sas_vs_rds_check("adexsum", a_in)
adexsum <- readRDS(read_path(a_in, "adexsum.rds")) %>%
  mutate(VISIT = if (exists("AVISIT")) AVISIT else "Overall") %>%
  filter(stringr::str_detect(PARAMCD, "RDINTE")) %>%
  select(USUBJID, VISIT, PARAMCD, AVAL)

adsl$active_trt <- factor(
  ifelse(adsl$TRT01A == "Placebo", " ", "Active Study Agent"),
  levels = c("Active Study Agent", " ")
)

# join data together
# User to add to PARAMLBL if any RDINTESx PARAMCDs exist due to different study agents
ex <- adexsum %>%
  inner_join(., adsl, by = c("USUBJID")) %>%
  mutate(
    PARAMLBL = case_when(
      PARAMCD == "RDINTE" ~ "Relative dose intensity~[super a] (%)"
    )
  )

################################################################################
# Define layout and build table:
################################################################################

lyt <- rtables::basic_table(
  top_level_section_div = " ",
  show_colcounts = TRUE,
  colcount_format = "N=xx"
) %>%
  split_cols_by("active_trt", split_fun = trim_levels_in_group("TRT01A"))

if (combined_active_trt == "Y") {
  lyt <- lyt %>%
    split_cols_by("TRT01A", split_fun = mysplit)
} else {
  lyt <- lyt %>%
    split_cols_by("TRT01A")
}

lyt <- lyt %>%
  split_rows_by(
    "VISIT",
    split_label = " ",
    split_fun = drop_split_levels,
    label_pos = "topleft",
    section_div = c(" ")
  ) %>%
  split_rows_by(
    "PARAMLBL",
    split_label = " ",
    split_fun = drop_split_levels,
    label_pos = "topleft",
    section_div = c(" "),
    nested = TRUE
  ) %>%
  analyze(
    "AVAL",
    show_labels = "hidden",
    nested = TRUE,
    indent_mod = 1L,
    afun = function(x) {
      list(
        "N" = rcell(length(x), format = jjcsformat_xx("xx"))
      )
    }
  ) %>%
  analyze(
    "AVAL",
    show_labels = "hidden",
    nested = TRUE,
    indent_mod = 3L,
    afun = function(x) {
      list(
        "Mean (SD)" = rcell(
          c(mean(x), sd(x)),
          format = jjcsformat_xx("xx.x (xx.xx)")
        ),
        "Median" = rcell(median(x), format = jjcsformat_xx("xx.x")),
        "Min, max" = rcell(
          c(min(x), max(x)),
          format = jjcsformat_xx("xx., xx.")
        ),
        "Interquartile range" = rcell(
          c(quantile(x, c(0.25, 0.75), type = 2)),
          format = jjcsformat_xx("xx.x, xx.x")
        )
      )
    }
  )

result <- build_table(lyt, ex, alt_counts_df = adsl)

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)
# Viewer(result)

################################################################################
# Convert to tbl file and output table
################################################################################
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
