###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lirada05.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Listing of Humoral Immunogenicity Through [Time Point] – [SAD/MAD][Part 1]
## Disclaimer:                This script is a direct copy of the corresponding Core Standard output identifier. For
##                            SAD/MAD specific changes, refer to tsfvit02b.r, lsidm05.r, and gsfvit02.r for examples of
##                            STUDYPRT filtering, COHORT handling, treatment column structure modifications, pooled
##                            placebo derivations, combined treatment columns, dose-level updates, and other
##                            output-specific structural differences as applicable.
## Author:                    C&SP Methodology
## Date:                      2026-06-18
## Input:                     adishum
## Output:                    lirada05.rtf
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

tblid <- "lirada05"
fileid <- write_path(opath, tblid)
popfl <- "IMFL" # Immunogenicity Analysis Set flag
trtvar <- "TRT01A"
trtvarnum <- "TRT01AN"

key_cols <- c("COL0", "COL1", "COL2", "COL3")
sort_cols <- c(trtvarnum, "COL1", "AVISITN", "ATPTN")
disp_cols <- paste0("COL", 4:7)

study_agent <- "active study agent"
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
    !!rlang::sym(popfl) == "Y"
  ) |>
  pull(
    USUBJID
  ) |>
  unique()

# base: visit-level, one row per subject-visit-timepoint
adishum_filtered_main <- adishum |>
  filter(
    USUBJID %in% main_subject_list,
    !is.na(AVISIT)
  ) |>
  select(
    !!trtvar,
    !!trtvarnum,
    AVISIT,
    AVISITN,
    USUBJID,
    ATPT,
    ATPTN,
    ADA_STATUS = ADATRES,
    NAB_STATUS = NABSTAT
  ) |>
  distinct()

# adishum_filtered on titer
adishum_filtered_titer <- adishum |>
  filter(
    stringi::stri_trans_toupper(PARAMCD) == "TITER",
    USUBJID %in% main_subject_list
  ) |>
  select(
    USUBJID,
    AVISIT,
    TITER_AVAL = AVAL
  )

# adishum_filtered on NABCNR
adishum_filtered_nabcnr <- adishum |>
  filter(
    stringi::stri_trans_toupper(PARAMCD) == "NABCNR",
    USUBJID %in% main_subject_list
  ) |>
  select(
    USUBJID,
    AVISIT,
    NABCNR_AVALC = AVALC
  )

# Section --- Build Dataset Base From Adishum ------------
# collect data after all joins
lsting <- adishum_filtered_main |>
  left_join(adishum_filtered_titer, by = c("USUBJID", "AVISIT")) |>
  left_join(adishum_filtered_nabcnr, by = c("USUBJID", "AVISIT"))

# create necessary columns for Listing
lsting <- lsting |>
  mutate(
    COL0 = explicit_na(.data[[trtvar]]),
    COL1 = explicit_na(USUBJID, ""),
    COL2 = explicit_na(ADA_STATUS, ""),
    # Optional Column: COL3/NAB_STATUS
    COL3 = explicit_na(NAB_STATUS, ""),
    COL4 = explicit_na(AVISIT, ""),
    # Optional Column: COL5/ATPT
    COL5 = explicit_na(ATPT, ""),
    COL6 = if_else(
      is.na(TITER_AVAL),
      "",
      formatC(TITER_AVAL, format = "fg")
    ),
    # Optional Column: COL7/NABCNR_AVALC
    COL7 = explicit_na(NABCNR_AVALC, "")
  )

# define column names for listing
lsting <- var_relabel(
  lsting,
  COL0 = "Treatment Group",
  COL1 = "Subject ID",
  COL2 = "Subject Status~[super a] for Treatment-Emergent ADA",
  # Optional Column: COL3/NAB_STATUS
  COL3 = "Subject Status for NAb",
  COL4 = "Visit",
  # Optional Column: COL5/ATPT
  COL5 = "Time Point",
  COL6 = sprintf(
    "Antibody to %s Status/Titer",
    study_agent
  ),
  # Optional Column: COL7/NABCNR_AVALC
  COL7 = "Sample Status for NAb"
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
    COL6 = fmt_config(
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
