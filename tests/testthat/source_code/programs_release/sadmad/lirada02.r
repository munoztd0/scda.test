###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lirada02.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Listing of Subjects Who Discontinued the Study by Treatment-emergent Antibodies to [Active
##                            Study Agent] Status – [SAD/MAD][Part 1]
## Disclaimer:                This script is a direct copy of the corresponding Core Standard output identifier. For
##                            SAD/MAD specific changes, refer to tsfvit02b.r, lsidm05.r, and gsfvit02.r for examples of
##                            STUDYPRT filtering, COHORT handling, treatment column structure modifications, pooled
##                            placebo derivations, combined treatment columns, dose-level updates, and other
##                            output-specific structural differences as applicable.
## Author:                    C&SP Methodology
## Date:                      2026-06-18
## Input:                     adishum
## Output:                    lirada02.rtf
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

tblid <- "lirada02"
fileid <- write_path(opath, tblid)
popfl <- "IMFL" # Immunogenicity Analysis Set flag
trtvar <- "TRT01A"
trtvarnum <- "TRT01AN"
adpc_paramcd_criteria <- "DOSE"

key_cols <- c("COL0", "COL1", "COL2")
sort_cols <- c(trtvarnum, "COL1", "AVISITN")
disp_cols <- paste0("COL", 0:8)
study_specific_col <- c(
  "COL7",
  "COL8"
)

study_agent <- "active study agent"
antibody_study_agent <- sprintf("antibodies to %s", study_agent)
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

# dataset adpc
adpc <- haven::read_sas(envsetup::read_path(a_in, "adpc.sas7bdat")) |>
  df_na() |>
  filter(
    toupper(PARAMCD) == toupper(adpc_paramcd_criteria),
    USUBJID %in% main_subject_list
  ) |>
  select(
    USUBJID,
    adpc_AVISIT = AVISIT,
    adpc_AVAL = AVAL
  )

# dataset adae
adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) |>
  df_na() |>
  filter(
    !is.na(DCSREAS),
    USUBJID %in% main_subject_list
  ) |>
  select(
    USUBJID,
    DCSREAS
  )

adishum_filtered_titer <- adishum |>
  filter(
    toupper(PARAMCD) == "TITER",
    USUBJID %in% main_subject_list
  ) |>
  select(
    !!trtvar,
    !!trtvarnum,
    AVISIT,
    AVISITN,
    USUBJID,
    TRTA,
    TITER_AVAL = AVAL
  )

# Section --- Build Dataset Base From Adishum ------------
# collect data after all joins
lsting <- adsl |>
  left_join(
    adishum_filtered_titer,
    by = "USUBJID"
  ) |>
  left_join(
    adpc,
    by = c(
      USUBJID = "USUBJID",
      AVISIT = "adpc_AVISIT"
    )
  )

# create necessary columns for Listing
lsting <- lsting |>
  mutate(
    COL0 = explicit_na(.data[[trtvar]]),
    COL1 = explicit_na(USUBJID, ""),
    COL2 = explicit_na(AVISIT, ""),
    COL3 = explicit_na(DCSREAS, ""),
    COL4 = explicit_na(TRTA, ""),
    COL5 = if_else(
      is.na(adpc_AVAL),
      "",
      formatC(adpc_AVAL, format = "fg")
    ),
    COL6 = if_else(
      is.na(TITER_AVAL),
      "",
      formatC(TITER_AVAL, format = "fg")
    ) #,
    # STUDY-SPECIFIC: replace later
    # COL7 = "STUDY-SPECIFIC",
    # COL8 = "STUDY-SPECIFIC"
  )

# define column names for listing
lsting <- var_relabel(
  lsting,
  COL0 = "Treatment Group",
  COL1 = "Subject ID",
  COL2 = "Visit",
  COL3 = "Reason for Discontinuation",
  COL4 = "Last Scheduled Study Treatment",
  COL5 = sprintf(
    "Matrix %s Conc. (\u03bcg/mL)",
    stringi::stri_trans_totitle(study_agent)
  ),
  COL6 = sprintf(
    "%s Status/Titer",
    stringi::stri_trans_totitle(antibody_study_agent)
  )
  # STUDY-SPECIFIC: update label as appropriate
  # COL7 = "[Success Criteria / Continuous Response]",
  # COL8 = "[Concomitant Medication]"
)

###############################################################################
# Build listing
###############################################################################

result <- rlistings::as_listing(
  df = lsting,
  key_cols = key_cols,
  disp_cols = disp_cols[!disp_cols %in% study_specific_col],
  sort_cols = sort_cols,
  round_type = "sas",
  col_formatting = list(
    COL5 = fmt_config(
      na_str = "",
      align = "right"
    ),
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
