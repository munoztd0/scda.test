###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lirada04.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Listing of Subjects With Samples Positive for Antibodies to [Active Study Agent] at
##                            Baseline – [SAD/MAD][Part 1]
## Disclaimer:                This script is a direct copy of the corresponding Core Standard output identifier. For
##                            SAD/MAD specific changes, refer to tsfvit02b.r, lsidm05.r, and gsfvit02.r for examples of
##                            STUDYPRT filtering, COHORT handling, treatment column structure modifications, pooled
##                            placebo derivations, combined treatment columns, dose-level updates, and other
##                            output-specific structural differences as applicable.
## Author:                    C&SP Methodology
## Date:                      2026-06-18
## Input:                     adishum
## Output:                    lirada04.rtf
## Remarks:
## R-functions:
## R-function Sample Call:
##
## Modification History:
##  Rev #:
##  Modified By:
##  Reporting Effort:
##  Date:
##  Description:
###############################################################################

###############################################################################
# Prep environment
###############################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)
library(dplyr)
library(rtables)
library(rlistings)
library(junco)

###############################################################################
# Define script level parameters
###############################################################################

tblid <- "lirada04"
fileid <- write_path(opath, tblid)
popfl <- "IMFL" # Immunogenicity Analysis Set flag
trtvar <- "TRT01A"
trtvarnum <- "TRT01AN"

parqual_value <- "XANOMELINE"
key_cols <- c("COL0", "COL1")
sort_cols <- c(trtvarnum, "COL1", "ADT")
disp_cols <- "COL2"

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

###############################################################################
# Process data
###############################################################################

adishum <- haven::read_sas(envsetup::read_path(a_in, "adishum.sas7bdat")) |>
  df_na()

# top level filter
main_subject_list <- adishum |>
  filter(
    !!rlang::sym(popfl) == "Y",
    PARQUAL == parqual_value,
    PARAMCD == "ADABL",
    AVALC == "Y"
  ) |>
  pull(
    USUBJID
  ) |>
  unique()

# adishum_filtered on ADABLT
adishum_filtered_adablt <- adishum |>
  filter(
    stringi::stri_trans_toupper(PARAMCD) == "ADABLT",
    USUBJID %in% main_subject_list
  ) |>
  select(
    !!trtvar,
    !!trtvarnum,
    ADT,
    USUBJID,
    ADABLT_AVAL = AVAL
  )

# Section --- Build Dataset Base From Adishum ------------
# collect data after all joins
lsting <- adishum_filtered_adablt

# create necessary columns for Listing
lsting <- lsting |>
  mutate(
    COL0 = explicit_na(.data[[trtvar]]),
    COL1 = explicit_na(USUBJID, ""),
    COL2 = if_else(
      is.na(ADABLT_AVAL),
      "",
      formatC(ADABLT_AVAL, format = "fg")
    )
  )

# define column names for listing
lsting <- var_relabel(
  lsting,
  COL0 = "Treatment Group",
  COL1 = "Subject ID",
  COL2 = "Baseline Titer"
)

###############################################################################
# Build listing
###############################################################################

result <- rlistings::as_listing(
  df = lsting,
  key_cols = key_cols,
  disp_cols = disp_cols,
  sort_cols = sort_cols,
  round_type = "sas",
  col_formatting = list(
    COL2 = fmt_config(
      na_str = "",
      align = "right"
    )
  )
)

###############################################################################
# Add titles and footnotes
###############################################################################

result <- set_titles(result, tab_titles)

###############################################################################
# Output listing
###############################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
