###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lsidem03.r
## R version:                 4.5.2
## junco version:             0.1.6
## Short Description:         Listing of Meal Intake – [SAD/MAD] [Part 1]
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, ml
## Output:                    lsidem03.rtf
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

tblid <- "LSIDEM03"
fileid <- write_path(opath, tblid)
popfl <- "FASFL"
trtvar <- "TRT01P"
key_cols <- paste0("COL", 0:4)
disp_cols <- paste0("COL", 0:7)
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

ml <- haven::read_sas(envsetup::read_path(d_in, "ml.sas7bdat")) %>%
  df_na()

adsl_ml <- ml %>%
  dplyr::inner_join(adsl, by = c("STUDYID", "USUBJID"))

lsting <- adsl_ml %>%
  dplyr::mutate(
    MLSTDT = as.Date(MLSTDTC, format = "%Y-%m-%d"),
    MLSTTM = substr(as.character(MLSTDTC), 12, 16),
    MLSTYR = dplyr::if_else(
      stringr::str_length(MLSTDTC) >= 4,
      substr(MLSTDTC, 1, 4),
      NA
    ),
    MLSTMO = month.abb[
      as.numeric(dplyr::if_else(
        stringr::str_length(sub("T.*", "", MLSTDTC)) >= 7 &
          substr(sub("T.*", "", MLSTDTC), 6, 7) != "--",
        substr(sub("T.*", "", MLSTDTC), 6, 7),
        NA
      ))
    ],
    MLSTDAY = dplyr::if_else(
      stringr::str_length(MLSTDTC) >= 10,
      substr(MLSTDTC, 9, 10),
      NA
    ),
    MLSTDY = dplyr::if_else(MLSTDT > TRTSDT, MLSTDT - TRTSDT + 1, MLSTDT - TRTSDT),
    MLSTDYL = dplyr::case_when(
      !is.na(MLSTDY) ~ as.character(MLSTDY),
      TRUE ~ "-"
    ),
    MLENDT = as.Date(MLENDTC, format = "%Y-%m-%d"),
    MLENTM = substr(as.character(MLENDTC), 12, 16),
    MLENYR = dplyr::if_else(
      stringr::str_length(MLENDTC) >= 4,
      substr(MLENDTC, 1, 4),
      NA
    ),
    MLENMO = month.abb[
      as.numeric(dplyr::if_else(
        stringr::str_length(sub("T.*", "", MLENDTC)) >= 7 &
          substr(sub("T.*", "", MLENDTC), 6, 7) != "--",
        substr(sub("T.*", "", MLENDTC), 6, 7),
        NA
      ))
    ],
    MLENDAY = dplyr::if_else(
      stringr::str_length(MLENDTC) >= 10,
      substr(MLENDTC, 9, 10),
      NA
    ),
    MLENDY = dplyr::if_else(MLENDT > TRTSDT, MLENDT - TRTSDT + 1, MLENDT - TRTSDT),
    MLENDYL = dplyr::case_when(
      !is.na(MLENDY) ~ as.character(MLENDY),
      TRUE ~ "-"
    ),
    MLDOSE = as.factor(MLDOSE)
  ) %>%
  unite(
    "MLSTDTL",
    MLSTDAY,
    MLSTMO,
    MLSTYR,
    sep = "",
    na.rm = TRUE,
    remove = FALSE
  ) %>%
  unite(
    "MLENDTL",
    MLENDAY,
    MLENMO,
    MLENYR,
    sep = "",
    na.rm = TRUE,
    remove = FALSE
  ) %>%
  mutate(
    COL0 = explicit_na(.data[[trtvar]], ""), # or cohort_var
    COL1 = explicit_na(USUBJID, ""),
    COL2 = explicit_na(MLTPT, ""),
    COL3 = case_when(
      (MLSTDTL != "" | !is.na(MLSTDTL)) & MLSTTM != "" & !is.na(MLSTDYL) ~
        paste0(toupper(MLSTDTL), concat_sep, MLSTTM, " (", MLSTDYL, ")"),
      MLSTDTL != "" & MLSTTM != "" & is.na(MLSTDYL) ~
        paste0(toupper(MLSTDTL), concat_sep, MLSTTM, ""),
      (MLSTDTL != "" | !is.na(MLSTDTL)) & MLSTTM == "" & !is.na(MLSTDYL) ~
        paste0(toupper(MLSTDTL), concat_sep, "--:--", " (", MLSTDYL, ")"),
      TRUE ~ ""
    ),
    COL4 = case_when(
      (MLENDTL != "" | !is.na(MLENDTL)) & MLENTM != "" & !is.na(MLENDYL) ~
        paste0(toupper(MLENDTL), concat_sep, MLENTM, " (", MLENDYL, ")"),
      MLENDTL != "" & MLENTM != "" & is.na(MLENDYL) ~
        paste0(toupper(MLENDTL), concat_sep, MLENTM, ""),
      (MLENDTL != "" | !is.na(MLENDTL)) & MLENTM == "" & !is.na(MLENDYL) ~
        paste0(toupper(MLENDTL), concat_sep, "--:--", " (", MLENDYL, ")"),
      TRUE ~ ""
    ),
    COL5 = stringr::str_to_sentence(explicit_na(MLTRT, "")),
    COL6 = explicit_na(MLDOSE, ""),
    COL7 = explicit_na(MLDOSU, "")
  ) %>%
  dplyr::arrange(COL0, COL1, COL3) %>%
  dplyr::select(dplyr::all_of(disp_cols))

lsting <- var_relabel(
  lsting,
  COL0 = "Treatment Group", # or Cohort
  COL1 = "Subject ID",
  COL2 = "Planned Time Point",
  COL3 = "Start Date/Time of Meal (Study Day~[super a])",
  COL4 = "End Date/Time of Meal (Study Day~[super a])",
  COL5 = "Type of Meal",
  COL6 = "Amount Consumed",
  COL7 = "Amount Consumed Unit"
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
