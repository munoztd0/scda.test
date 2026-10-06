###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lsidem05.r
## R version:                 4.5.2
## junco version:             0.1.6
## Short Description:         Listing of Comments – [SAD/MAD] [Part 1]
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, co
## Output:                    lsidem05.rtf
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
library(tidyr)
library(rtables)
library(rlistings)
library(junco)

###############################################################################
# Define script level parameters
###############################################################################

tblid <- "LSIDEM05"
fileid <- write_path(opath, tblid)
popfl <- "FASFL"
trtvar <- "TRT01P"
key_cols <- paste0("COL", 0:2)
disp_cols <- paste0("COL", 0:4)
concat_sep <- " / "
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
studypart_var <- "PARTC" # need to be changed to STUDYPRT variable
studyprt_val <- "PART 1"
cohort_var <- "COHORTC"
actarm_str <- "SAD"


###############################################################################
# Process data
###############################################################################


adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) %>%
  df_na() %>%
  dplyr::filter(
    !!rlang::sym(studypart_var) == studyprt_val,
    grepl(actarm_str, ACTARM, ignore.case = TRUE),
    !!rlang::sym(popfl) == "Y"
  ) %>%
  dplyr::mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = unique(sort(.data[[trtvar]]))
    ),
    !!rlang::sym(cohort_var) := factor(
      .data[[cohort_var]],
      levels = unique(.data[[cohort_var]])
    )
  )

co <- haven::read_sas(envsetup::read_path(d_in, "co.sas7bdat")) %>%
  df_na()

co <- co %>%
  dplyr::mutate(
    COVAL_COMBINED = apply(
      dplyr::pick(starts_with("COVAL")),
      1,
      function(r) {
        vals <- as.character(r)
        vals <- vals[!is.na(vals) & vals != ""]
        paste(vals, collapse = "")
      }
    )
  )

adsl_co <- co %>%
  dplyr::inner_join(adsl, by = c("STUDYID", "USUBJID"))

lsting <- adsl_co %>%
  dplyr::mutate(
    COL0 = explicit_na(.data[[trtvar]], ""), # or cohort_var
    COL1 = explicit_na(USUBJID, ""),
    COL2 = explicit_na(RDOMAIN, ""),
    COL3 = explicit_na(COREF, ""),
    COL4 = stringr::str_to_sentence(explicit_na(COVAL_COMBINED, ""))
  )

lsting <- var_relabel(
  lsting,
  COL0 = "Treatment Group", # or Cohort
  COL1 = "Subject ID",
  COL2 = "Related Domain",
  COL3 = "Comment Reference",
  COL4 = "Comment"
)

###############################################################################
# Build listing
###############################################################################

result <- rlistings::as_listing(
  df = lsting,
  key_cols = key_cols,
  disp_cols = disp_cols,
  round_type = "sas"
)

###############################################################################
# Add titles and footnotes
###############################################################################

result <- set_titles(result, tab_titles)

###############################################################################
# Output listing
###############################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
