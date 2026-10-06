###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lsfvit01.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Listing of Subjects With Clinically Important Vital Signs – [SAD/MAD][Part 1]
## Disclaimer:                This script is a direct copy of the corresponding Core Standard output identifier. For
##                            SAD/MAD specific changes, refer to tsfvit02b.r, lsidm05.r, and gsfvit02.r for examples of
##                            STUDYPRT filtering, COHORT handling, treatment column structure modifications, pooled
##                            placebo derivations, combined treatment columns, dose-level updates, and other
##                            output-specific structural differences as applicable.
##
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     advs
## Output:                    lsfvit01.rtf
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
library(tidytlg)


###############################################################################
# Define script level parameters
###############################################################################

tblid <- "lsfvit01"
fileid <- write_path(opath, tblid)
popfl <- "SAFFL"
trtvar <- "TRT01A"
key_cols <- c("COL0", "COL1", "COL2", "COL3")
sort_cols <- c("COL0", "COL1", "COL2", "COL3")
disp_cols <- paste0("COL", 0:11)
concat_sep <- " / "
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Parameters to be included by default
param_cds <- c("SYSBP", "DIABP", "PULSE", "RESP", "TEMP", "WEIGHT", "HEIGHT")

###############################################################################
# Process data
###############################################################################

advs <- haven::read_sas(envsetup::read_path(a_in, "advs.sas7bdat")) |>
  df_na() |>
  filter(!!rlang::sym(popfl) == "Y" & PARAMCD %in% param_cds) |>
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c("Xanomeline Low Dose", "Xanomeline High Dose", "Placebo")
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
    ),
    PARAM = factor(
      .data$PARAM,
      levels = unique(.data[['PARAM']])[order(unique(.data[['PARAMN']]))]
    ),
    AVISIT = factor(
      .data[['AVISIT']],
      levels = unique(.data[['AVISIT']])[order(unique(.data[['AVISITN']]))]
    )
  )

advs_crit <- advs |>
  filter(CRIT7FL == 'Y' | CRIT8FL == 'Y') |>
  select(STUDYID, USUBJID, PARAMCD) |>
  distinct()

advs_dig <- tidytlg:::make_precision_data(
  df = advs,
  decimal = 4,
  precisionby = "PARAMCD",
  precisionon = "AVAL"
) |>
  rename(c(VALDIGMAX = "decimal"))

advs_list <- advs_crit |>
  inner_join(
    advs,
    by = c(
      "STUDYID" = "STUDYID",
      "USUBJID" = "USUBJID",
      "PARAMCD" = "PARAMCD"
    )
  ) |>
  inner_join(advs_dig, by = c("PARAMCD" = "PARAMCD"))

lsting <- advs_list |>
  mutate(
    AGE = explicit_na(as.character(AGE), ""),
    SEX = explicit_na(SEX, ""),
    RACE = explicit_na(RACE, ""),
    ADT = ifelse(
      nchar(as.character(ADT)) == 10,
      toupper(format(ADT, "%d%b%Y")),
      ""
    ),
    ATM = ifelse(!is.na(ADTM), substr(as.character(ADTM), 12, 16), ""),
    ADYN = ifelse(!is.na(ADY), ADY, NA),
    ADY = ifelse(!is.na(ADY), ADY, "--"),
    VAL_RES = case_when(
      is.na(VALDIGMAX) & !is.na(AVAL) ~ tidytlg::roundSAS(AVAL, digits = 0, as_char = TRUE, na_char = NULL),
      VALDIGMAX == 0 & !is.na(AVAL) ~ tidytlg::roundSAS(AVAL, digits = 0, as_char = TRUE, na_char = NULL),
      VALDIGMAX == 1 & !is.na(AVAL) ~ tidytlg::roundSAS(AVAL, digits = 1, as_char = TRUE, na_char = NULL),
      VALDIGMAX == 2 & !is.na(AVAL) ~ tidytlg::roundSAS(AVAL, digits = 2, as_char = TRUE, na_char = NULL),
      VALDIGMAX == 3 & !is.na(AVAL) ~ tidytlg::roundSAS(AVAL, digits = 3, as_char = TRUE, na_char = NULL),
      VALDIGMAX >= 4 & !is.na(AVAL) ~ tidytlg::roundSAS(AVAL, digits = 4, as_char = TRUE, na_char = NULL),
      !is.na(AVALC) ~ AVALC
    ),
    VAL_HL = case_when(
      !is.na(ANRIND) ~ substr(ANRIND, 1, 1),
      .default = NA
    ),
    VAL = case_when(
      !is.na(VAL_RES) & !is.na(VAL_HL) ~ paste(VAL_RES, VAL_HL, sep = " "),
      !is.na(VAL_RES) & is.na(VAL_HL) ~ VAL_RES,
      .default = NA
    ),
    CRIT7L = ifelse(CRIT7FL == "Y", as.character(CRIT7), NA),
    CRIT8L = ifelse(CRIT8FL == "Y", as.character(CRIT8), NA)
  ) |>
  unite(
    "CRITL",
    CRIT7L,
    CRIT8L,
    sep = ", ",
    na.rm = TRUE,
    remove = FALSE
  ) |>
  mutate(
    COL0 = explicit_na(.data[[trtvar]], ""),
    COL1 = explicit_na(USUBJID, ""),
    COL2 = paste(AGE, SEX, RACE, sep = concat_sep),
    COL3 = explicit_na(PARAM, ""),
    # Optional Variable: ATM
    COL4 = case_when(
      ADT == "" ~ "",
      ADT != "" & ATM != "" & ADY != "" ~ paste0(ADT, concat_sep, ATM, " (", ADY, ")"),
      ADT != "" & ATM == "" & ADY != "" ~ paste0(ADT, concat_sep, "--:--", " (", ADY, ")"),
      ADT != "" & ATM != "" & ADY == "" ~ paste0(ADT, concat_sep, ATM, " (-)"),
      ADT != "" & ATM == "" & ADY == "" ~ paste0(ADT, concat_sep, "--:--", " (-)"),
    ),
    COL5 = explicit_na(AVISIT, ""),
    # Optional Column: COL6/ATPT
    COL6 = explicit_na(stringr::str_to_sentence(ATPT), ""),
    COL7 = VAL,
    COL8 = explicit_na(VSCLSIG, ""),
    COL9 = case_when(
      is.na(CHG) ~ "",
      is.na(VALDIGMAX) & !is.na(AVAL) & !is.na(CHG) ~
        tidytlg::roundSAS(CHG, digits = 0, as_char = TRUE, na_char = NULL),
      VALDIGMAX == 0 & !is.na(CHG) ~ tidytlg::roundSAS(CHG, digits = 0, as_char = TRUE, na_char = NULL),
      VALDIGMAX == 1 & !is.na(CHG) ~ tidytlg::roundSAS(CHG, digits = 1, as_char = TRUE, na_char = NULL),
      VALDIGMAX == 2 & !is.na(CHG) ~ tidytlg::roundSAS(CHG, digits = 2, as_char = TRUE, na_char = NULL),
      VALDIGMAX == 3 & !is.na(CHG) ~ tidytlg::roundSAS(CHG, digits = 3, as_char = TRUE, na_char = NULL),
      VALDIGMAX >= 4 & !is.na(CHG) ~ tidytlg::roundSAS(CHG, digits = 4, as_char = TRUE, na_char = NULL)
    ),
    # Optional Column: COL9/CRITy/ATOXGR
    COL10 = explicit_na(CRITL, ""),
    # COL10 = explicit_na(ATOXGR, ""),
    COL11 = explicit_na(TRTEMFL, "")
  ) |>
  arrange(COL0, COL1, COL2, COL3, !is.na(ADYN), ADYN, ADTM)

lsting <- var_relabel(
  lsting,
  COL0 = "Treatment Group",
  COL1 = "Subject ID",
  COL2 = paste("Age (years)", "Sex", "Race", sep = concat_sep),
  COL3 = "Vital Sign (unit)",
  # Optional Variable: ATM
  COL4 = paste(
    "Assessment Date",
    "Time (Study Day~[super a])",
    sep = concat_sep
  ),
  COL5 = "Visit",
  # Optional Column: COL6/ATPT
  COL6 = "Time Point",
  COL7 = "Result",
  COL8 = "Clinically Significant?",
  COL9 = "Change From Baseline",
  # Optional Column: COL10/CRITy/ATOXGR
  # COL10 = "Grade"
  COL10 = "Criteria",
  COL11 = "Treatment-emergent?"
)

###############################################################################
# Build listing
###############################################################################

result <- rlistings::as_listing(
  df = lsting,
  key_cols = key_cols,
  sort_cols = sort_cols,
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
