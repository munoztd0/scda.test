###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

###############################################################################
## Original Reporting Effort: Standards
## Program Name:               lsflab02.r
## R version:                  4.5.2
## junco version:              0.1.3
## Short Description:          Program to create lsflab02: Listing of Liver Function Results
##                             for Subjects With ≥1 On-treatment [AST or ALT]
##                             Values ≥3x ULN, TBILI Values ≥2x ULN or ALP Values ≥2x ULN
## Author:                     C&SP Methodology
## Date:                       2026-09-30
## Input:                      addili
## Output:                     lsflab02.rtf
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
library(rlistings)
library(junco)

###############################################################################
# Define script level parameters
###############################################################################

tblid <- "LSFLAB02"
fileid <- write_path(opath, tblid)
popfl <- "SAFFL"
trtvar <- "TRT01A"
key_cols <- c("COL0", "COL1", "COL2")
sort_cols <- c("COL0", "COL1", "COL2", "ADT", "ADY")
disp_cols <- paste0("COL", 0:8)
concat_sep <- " / "
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

central_lab <- c("CENTRAL")
###############################################################################
# Process data
###############################################################################

addili_raw <- haven::read_sas(envsetup::read_path(a_in, "addili.sas7bdat")) |>
  df_na()

addili <- addili_raw |>
  filter(
    !!rlang::sym(popfl) == "Y",
    PARAMCD %in% c("ALT", "AST", "ALP", "BILI"),
    !is.na(.data[[trtvar]]),
    is.na(DTYPE)
  ) |>
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c("Xanomeline Low Dose", "Xanomeline High Dose", "Placebo")
    ),
    SEX = factor(
      case_when(SEX == "M" ~ "Male", SEX == "F" ~ "Female", TRUE ~ SEX),
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
    AVISIT = factor(
      .data[["AVISIT"]],
      levels = unique(.data[["AVISIT"]])[order(unique(.data[["AVISITN"]]))]
    )
  )

addili_crit <- addili |>
  filter(ONTRTFL == "Y" & CRIT1FL == "Y") |>
  select(STUDYID, USUBJID) |>
  distinct()

addili_joined <- addili_crit |>
  inner_join(
    addili,
    by = c(
      "STUDYID" = "STUDYID",
      "USUBJID" = "USUBJID"
    )
  )

addili_list <- addili_joined |>
  mutate(
    R2ANRHI_CHAR = sprintf("%.2f", R2ANRHI),
    # For BILI records of Hy's Law subjects where the visit-level flag is set,
    # concatenate R2ANRHI value with the flag "^"
    # Non-central lab values are additionally suffixed with "*"
    VAL = paste0(
      R2ANRHI_CHAR,
      if_else(!is.na(LBNAM) & toupper(LBNAM) != toupper(central_lab), "*", ""),
      if_else(
        PARAMCD == "BILI" & !is.na(ANL05FL) & ANL05FL == 'Y' & !is.na(CRIT1FL) & CRIT1FL == 'Y',
        "^",
        ""
      ),
      if_else(
        PARAMCD == "BILI" & !is.na(ANL07FL) & ANL07FL == 'Y' & !is.na(CRIT1FL) & CRIT1FL == 'Y',
        "#",
        ""
      )
    )
  ) |>
  group_by(
    STUDYID,
    USUBJID,
    TRT01A,
    AGE,
    SEX,
    RACE,
    AVISIT,
    ADT,
    ADY,
    PARAMCD
  ) |>
  summarise(
    VAL = paste(unique(na.omit(VAL)), collapse = ", "),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = PARAMCD,
    values_from = VAL
  )

lsting <- addili_list |>
  mutate(
    AGE = as.character(AGE),
    SEX = as.character(SEX),
    RACE = as.character(RACE),
    ADT_CHAR = ifelse(!is.na(ADT), toupper(format(ADT, "%d%b%Y")), ""),
    ADY_CHAR = ifelse(!is.na(ADY), as.character(ADY), "--"),
    COL0 = .data[[trtvar]],
    COL1 = explicit_na(USUBJID, ""),
    COL2 = paste(AGE, SEX, RACE, sep = concat_sep),
    COL3 = case_when(
      ADT_CHAR == "" ~ "",
      ADT_CHAR != "" & ADY_CHAR != "" ~ paste0(ADT_CHAR, " (", ADY_CHAR, ")"),
      ADT_CHAR != "" & ADY_CHAR == "" ~ paste0(ADT_CHAR, " (-)")
    ),
    COL4 = AVISIT,
    COL5 = if ("ALT" %in% names(addili_list)) ALT else "",
    COL6 = if ("AST" %in% names(addili_list)) AST else "",
    COL7 = if ("ALP" %in% names(addili_list)) ALP else "",
    COL8 = if ("BILI" %in% names(addili_list)) BILI else ""
  ) |>
  mutate(across(c(COL5, COL6, COL7, COL8), ~ replace_na(., ""))) |>
  arrange(COL0, COL1, ADT, ADY)

lsting <- var_relabel(
  lsting,
  COL0 = "Treatment Group",
  COL1 = "Subject ID",
  COL2 = paste("Age (years)", "Sex", "Race", sep = concat_sep),
  COL3 = "Assessment Date (Study Day~[super a])",
  COL4 = "Visit",
  COL5 = "ALT\\line(x ULN)",
  COL6 = "AST\\line(x ULN)",
  COL7 = "ALP\\line(x ULN)",
  COL8 = "TBILI\\line(x ULN)"
)


###############################################################################
# Build listing object
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

tt_to_tlgrtf(string_map = string_map, tt = 
  result,
  file = fileid,
  orientation = "landscape"
)
