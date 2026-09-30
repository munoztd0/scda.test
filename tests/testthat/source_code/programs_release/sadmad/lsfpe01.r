###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lspe01.r
## R version:                 4.5.2
## junco version:             0.1.6
## Short Description:         Listing of Subjects With Abnormal Physical
##                            Examination Results – [SAD/MAD] [Part 1]
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, pe, supppe
## Output:                    lspe01.rtf
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

tblid <- "LSFPE01"
fileid <- write_path(opath, tblid)
popfl <- "SAFFL"
trtvar <- "TRT01A"
key_cols <- paste0("COL", 0:2)
disp_cols <- paste0("COL", 0:7)
concat_sep <- " / "
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
studypart_var <- "PARTC"
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
    ),
    SEX = factor(
      case_when(SEX == "M" ~ "Male", SEX == "F" ~ 'Female', TRUE ~ SEX),
      levels = c("Male", "Female", "Intersex", "Unknown")
    ),
    RACE = factor(
      case_when(
        RACE == "AMERICAN INDIAN OR ALASKA NATIVE" ~ "American Indian or Alaska Native",
        RACE == "ASIAN" ~ "Asian",
        RACE == "BLACK OR AFRICAN AMERICAN" ~ "Black or African American",
        RACE == "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" ~ "Native Hawaiian or other Pacific Islander",
        RACE == "WHITE" ~ "White",
        RACE == "MULTIPLE" ~ "Multiple",
        RACE == "NOT REPORTED" ~ "Not reported",
        RACE == "UNKNOWN" ~ "Unknown",
        RACE == "OTHER" ~ "Other"
      ),
      levels = c(
        "American Indian or Alaska Native",
        "Asian",
        "Black or African American",
        "Native Hawaiian or other Pacific Islander",
        "White",
        "Multiple",
        "Not reported",
        "Unknown",
        "Other"
      )
    )
  ) %>%
  select(
    USUBJID,
    !!rlang::sym(trtvar),
    !!rlang::sym(cohort_var),
    AGE,
    AGEU,
    SEX,
    RACE
  )

pe <- haven::read_sas(envsetup::read_path(d_in, "pe.sas7bdat")) %>%
  df_na() %>%
  filter(PESTRESC == "ABNORMAL") %>%
  select(USUBJID, VISIT, PEDTC, PETEST, PEORRES)

supppe <- haven::read_sas(envsetup::read_path(d_in, "supppe.sas7bdat")) %>%
  df_na() %>%
  tidyr::pivot_wider(
    names_from = "QNAM",
    values_from = "QVAL"
  ) %>%
  select(USUBJID, PECLSIG) %>%
  distinct()

lsting <- pe %>%
  inner_join(adsl, by = "USUBJID") %>%
  inner_join(supppe, by = "USUBJID") %>%
  mutate(
    PEDTC = as.Date(PEDTC, format = "%Y-%m-%d"),
    PEDTM = substr(as.character(PEDTC), 12, 16),
    PEDTYR = dplyr::if_else(
      stringr::str_length(PEDTC) >= 4,
      substr(PEDTC, 1, 4),
      NA
    ),
    PEDTMO = month.abb[
      as.numeric(dplyr::if_else(
        stringr::str_length(sub("T.*", "", PEDTC)) >= 7 &
          substr(sub("T.*", "", PEDTC), 6, 7) != "--",
        substr(sub("T.*", "", PEDTC), 6, 7),
        NA
      ))
    ],
    PEDTDAY = dplyr::if_else(
      stringr::str_length(PEDTC) >= 10,
      substr(PEDTC, 9, 10),
      NA
    )
  ) %>%
  unite(
    "PEDT",
    PEDTDAY,
    PEDTMO,
    PEDTYR,
    sep = "",
    na.rm = TRUE,
    remove = FALSE
  ) %>%
  mutate(
    AGE = explicit_na(as.character(AGE), ""),
    AGEU = explicit_na(AGEU, ""),
    SEX = explicit_na(SEX, ""),
    RACE = explicit_na(RACE, ""),
    COL0 = explicit_na(.data[[trtvar]], ""), # or cohort_var
    COL1 = explicit_na(USUBJID, ""),
    COL2 = paste(AGE, SEX, RACE, sep = concat_sep),
    COL3 = stringr::str_to_sentence(
      explicit_na(as.character(VISIT), "")
    ),
    COL4 = toupper(PEDT),
    COL5 = explicit_na(as.character(PETEST), ""),
    COL6 = stringr::str_to_sentence(
      explicit_na(as.character(PEORRES), "")
    ),
    COL7 = explicit_na(as.character(PECLSIG), "")
  ) %>%
  arrange(COL0, COL1, PEDTC)

lsting <- var_relabel(
  lsting,
  COL0 = "Treatment Group", # or Cohort
  COL1 = "Subject ID",
  COL2 = paste("Age (years)", "Sex", "Race", sep = concat_sep),
  COL3 = "Visit",
  COL4 = "Assessment Date",
  COL5 = "Body System",
  COL6 = "Abnormal Findings",
  COL7 = "Clinically Significant"
)

###############################################################################
# Build listing
###############################################################################

result <- rlistings::as_listing(
  df = lsting,
  key_cols = key_cols,
  disp_cols = disp_cols
)

###############################################################################
# Add titles and footnotes:
###############################################################################

result <- set_titles(result, tab_titles)

###############################################################################
# Output listing
###############################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
