###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefsubtx01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefsubtx01: Best Response to First
##                            Subsequent Antimyeloma Therapy
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADEFF
## Output:                    TEFSUBTX01.rtf
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


library(dplyr)
library(rtables)
library(junco)
library(haven)

################################################################################
# Define script level parameters:
################################################################################

################################################################################
# - Define output ID and file location
# - Define treatment variable used (default=TRT01P)
# - Define population flag used (default=FASFL)
# - Define control treatment arm
################################################################################

tblid <- "TEFSUBTX01"
fileid <- write_path(opath, tblid)

titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
titles$title <- titles$title[1]

trtvar <- "TRT01P"
popfl <- "RANDFL"
ctrl_grp <- "Dummy B - DPd"

################################################################################
# Read in ADSL and ADEFF datasets
################################################################################

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) %>%
  filter(!!rlang::sym(popfl) == "Y") %>%
  mutate(!!trtvar := as.factor(.data[[trtvar]])) %>%
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl))

adsl$colspan_trt <- factor(
  ifelse(adsl[[trtvar]] == ctrl_grp, " ", "Active Study Agent"),
  levels = c("Active Study Agent", " ")
)

adrs <- haven::read_sas(read_path(a_in, "adrs.sas7bdat")) %>%
  filter(PARAMCD == "IRCCONSR") %>%
  mutate(
    PARAMCD = as.factor(PARAMCD),
    AVALC = case_when(
      AVALC == "SCR" ~ "Stringent complete response (sCR)",
      AVALC == "CR" ~ "Complete response (CR)",
      AVALC == "VGPR" ~ "Very good partial response (VGPR)",
      AVALC == "PR" ~ "Partial response (PR)",
      AVALC == "MR" ~ "Minimal response (MR)",
      AVALC == "SD" ~ "Stable disease (SD)",
      AVALC == "PD" ~ "Progressive disease (PD)",
      AVALC == "NE" ~ "Not applicable (NA)",
      AVALC == " " ~ "Unknown"
    ),
    AVALC = factor(
      AVALC,
      levels = c(
        "Stringent complete response (sCR)",
        "Complete response (CR)",
        "Very good partial response (VGPR)",
        "Partial response (PR)",
        "Minimal response (MR)",
        "Stable disease (SD)",
        "Progressive disease (PD)",
        "Not applicable (NA)",
        "Unknown"
      )
    )
  ) %>%
  mutate(RESPONSE = case_when(PARAMCD == "IRCCONSR" ~ "Y")) %>%
  select(USUBJID, RESPONSE, AVALC)

# Create best response if needed - for now just take first until we find out if this should already
# be derived in the ADaM dataset
adrsbr <- adrs %>%
  mutate(
    resporder = case_when(
      AVALC == "Stringent complete response (sCR)" ~ 1,
      AVALC == "Complete response (CR)" ~ 2,
      AVALC == "Very good partial response (VGPR)" ~ 3,
      AVALC == "Partial response (PR)" ~ 4,
      AVALC == "Minimal response (MR)" ~ 5,
      AVALC == "Stable disease (SD)" ~ 6,
      AVALC == "Progressive disease (PD)" ~ 7,
      AVALC == "Not applicable (NA)" ~ 8,
      AVALC == "Unknown" ~ 9
    )
  )

adrsbr <- adrsbr %>%
  arrange(USUBJID, resporder) %>%
  group_by(USUBJID) %>%
  slice(1) %>%
  ungroup()

# join data together
eff <- adrsbr %>% inner_join(., adsl, by = c("USUBJID"))

################################################################################
# Define layout and build table:
################################################################################

lyt <- rtables::basic_table(
  top_level_section_div = " ",
  show_colcounts = TRUE,
  colcount_format = "N=xx"
) %>%
  split_cols_by("colspan_trt", split_fun = trim_levels_in_group(trtvar)) %>%
  split_cols_by(trtvar) %>%
  add_overall_col("Total") %>%
  analyze(
    "AVALC",
    var_labels = "Best response to first subsequent antimyeloma therapy~[super a]",
    show_labels = "visible",
    afun = a_freq_j,
    extra_args = list(
      denom = "n_df",
      .stats = c("n_df", "count_unique_denom_fraction")
    ),
    indent_mod = 0L
  ) %>%
  append_topleft("Response")

result <- build_table(lyt, eff, alt_counts_df = adsl)

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, titles)

################################################################################
# Convert to tbl file and output table
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
