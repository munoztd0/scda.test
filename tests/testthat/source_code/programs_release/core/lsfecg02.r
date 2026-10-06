###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lsfecg02.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create lsfecg02: Listing of ECG Values
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adeg
## Output:                    lsfecg02.rtf
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
library(tidytlg)

###############################################################################
# Define script level parameters
###############################################################################

tblid <- "LSFECG02"
fileid <- write_path(opath, tblid)
popfl <- "SAFFL"
trtvar <- "TRT01A"
key_cols <- c("COL0", "COL1", "COL2", "COL3")
sort_cols <- c("COL0", "COL1", "COL2", "COL3")
disp_cols <- paste0("COL", 0:10)
concat_sep <- " / "
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()


###############################################################################
# Process data
###############################################################################

adeg <- haven::read_sas(envsetup::read_path(a_in, "adeg.sas7bdat")) |>
  df_na() |>
  filter(!!rlang::sym(popfl) == "Y" & PARAMCD != "EGALL") |>
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

adeg_dig <- tidytlg:::make_precision_data(
  df = adeg,
  decimal = 4,
  precisionby = "PARAMCD",
  precisionon = "AVAL"
) |>
  rename(c(VALDIGMAX = "decimal"))

adeg_list <- adeg |>
  inner_join(adeg_dig, by = c("PARAMCD" = "PARAMCD"))

lsting <- adeg_list |>
  mutate(
    AGE = explicit_na(as.character(AGE), ""),
    SEX = explicit_na(SEX, ""),
    RACE = explicit_na(RACE, ""),
    ADT = ifelse(
      nchar(as.character(ADT)) == 10,
      toupper(format(ADT, "%d%b%Y")),
      "---------"
    ),
    ATM = ifelse(!is.na(ADTM), substr(as.character(ADTM), 12, 16), "--:--"),
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
    VAL_CS = case_when(
      EGCLSIG == "Y" & PARAMCD == "INTP" ~ "CS",
      EGCLSIG == "N" & PARAMCD == "INTP" ~ "NCS",
      .default = NA
    ),
    VAL = case_when(
      !is.na(VAL_RES) & !is.na(VAL_CS) ~
        paste(
          VAL_RES,
          VAL_CS,
          sep = " "
        ),
      !is.na(VAL_RES) & is.na(VAL_CS) ~ VAL_RES,
      .default = NA
    ),
    CRIT = case_when(
      (is.na(CRIT1) | CRIT1FL == "N") & (is.na(CRIT2) | CRIT2FL == "N") ~ "",
      !is.na(CRIT1) & CRIT1FL == "Y" ~ CRIT1,
      !is.na(CRIT2) & CRIT2FL == "Y" ~ CRIT2,
    ),
    TREM_FL = case_when(
      TRTEMFL == "Y" ~ "Yes",
      .default = NA
    ),
    COL0 = explicit_na(.data[[trtvar]], ""),
    COL1 = explicit_na(USUBJID, ""),
    COL2 = paste(AGE, SEX, RACE, sep = concat_sep),
    COL3 = explicit_na(PARAM, ""),
    # Optional Variable: ATM
    COL4 = paste(ADT, concat_sep, ATM, " (", ADY, ")", sep = ""),
    COL5 = explicit_na(AVISIT, ""),
    # Optional Column: COL6/ATPT
    COL6 = explicit_na(ATPT, ""),
    COL7 = explicit_na(VAL, ""),
    # Optional Column: COL8/CHG
    COL8 = case_when(
      is.na(CHG) ~ "",
      is.na(VALDIGMAX) & !is.na(AVAL) & !is.na(CHG) ~
        tidytlg::roundSAS(CHG, digits = 0, as_char = TRUE, na_char = NULL),
      VALDIGMAX == 0 & !is.na(CHG) ~ tidytlg::roundSAS(CHG, digits = 0, as_char = TRUE, na_char = NULL),
      VALDIGMAX == 1 & !is.na(CHG) ~ tidytlg::roundSAS(CHG, digits = 1, as_char = TRUE, na_char = NULL),
      VALDIGMAX == 2 & !is.na(CHG) ~ tidytlg::roundSAS(CHG, digits = 2, as_char = TRUE, na_char = NULL),
      VALDIGMAX == 3 & !is.na(CHG) ~ tidytlg::roundSAS(CHG, digits = 3, as_char = TRUE, na_char = NULL),
      VALDIGMAX >= 4 & !is.na(CHG) ~ tidytlg::roundSAS(CHG, digits = 4, as_char = TRUE, na_char = NULL)
    ),
    # Optional Column: COL9/CRITy
    COL9 = explicit_na(CRIT, ""),
    COL10 = explicit_na(TREM_FL, "")
  ) |>
  arrange(COL0, COL1, COL2, COL3, !is.na(ADYN), ADYN, ADTM)

lsting <- var_relabel(
  lsting,
  COL0 = "Treatment Group",
  COL1 = "Subject ID",
  COL2 = paste("Age (years)", "Sex", "Race", sep = concat_sep),
  COL3 = "ECG Parameter (unit)",
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
  # Optional Column: COL8/CHG
  COL8 = "Change From Baseline",
  # Optional Column: COL9/CRITy
  COL9 = "Criteria",
  COL10 = "Treatment-emergent?"
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
# If resulting Listing output is too large of a file (>20MB) then the listing
# should be split into multiple parts
# The split below is based on the treatment groups
# Update as-needed for your study
result1 <- result |>
  filter(toupper(.data[[trtvar]]) == "XANOMELINE LOW DOSE")

tt_to_tlgrtf(string_map = string_map, tt = 
  result1,
  file = paste0(fileid, "PART1OF3"),
  orientation = "landscape"
)

result2 <- result |>
  filter(toupper(.data[[trtvar]]) == "XANOMELINE HIGH DOSE")

tt_to_tlgrtf(string_map = string_map, tt = 
  result2,
  file = paste0(fileid, "PART2OF3"),
  orientation = "landscape"
)

result3 <- result |>
  filter(toupper(.data[[trtvar]]) == "PLACEBO")

tt_to_tlgrtf(string_map = string_map, tt = 
  result3,
  file = paste0(fileid, "PART3OF3"),
  orientation = "landscape"
)
