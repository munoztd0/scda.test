###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lsids02.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create lsids02: Listing of Subjects Who
##                            Discontinued Study Participation Prematurely
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl
## Output:                    lsids02.rtf
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
library(tern)
library(dplyr)
library(rtables)
library(rlistings)
library(junco)

###############################################################################
# Define script level parameters
###############################################################################

tblid <- "LSIDS02"
fileid <- write_path(opath, tblid)
popfl <- "FASFL"
trtvar <- "TRT01P"
key_cols <- c("COL0", "COL1")
disp_cols <- paste0("COL", 0:9)
concat_sep <- " / "
tab_titles <- list(title = "Dummy Title",
                     subtitles = NULL,
                     main_footer = "Dummy Note: On-treatment is defined as ~{optional treatment-emergent}")

###############################################################################
# Process data
###############################################################################

adsl <- adsl_jnj %>%
  filter(!!rlang::sym(popfl) == "Y" & EOSSTT %in% c("DISCONTINUED")) %>%
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c(
        "Xanomeline Low Dose",
        "Xanomeline High Dose",
        "Placebo"
      )
    ),
    SEX = factor(
      case_when(
        SEX == "F" ~ "Female",
        SEX == "M" ~ "Male"
      ),
      levels = c("Female", "Male")
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
  )

adexsum <- adexsum_jnj %>%
  filter(PARAMCD == "CUMDOSE") %>%
  select(STUDYID, USUBJID, PARAMCD, PARAM, AVAL)

adsl_adexsum <- left_join(
  adsl,
  adexsum,
  by = c(
    "STUDYID" = "STUDYID",
    "USUBJID" = "USUBJID"
  )
)

lsting <- adsl_adexsum %>%
  mutate(
    AGE = explicit_na(as.character(AGE), ""),
    SEX = explicit_na(SEX, ""),
    RACE = explicit_na(RACE, ""),
    AVAL = explicit_na(as.character(AVAL), ""),
    AVALU = case_when(
      !is.na(AVAL) ~ stringr::str_extract(PARAM, "(?<=\\()([^()]*?)(?=\\)[^()]*$)"),
      is.na(AVAL) ~ ""
    ),
    DCSREAS = explicit_na(DCSREAS, ""),
    DCSREASP = explicit_na(DCSREASP, ""),
    COL0 = explicit_na(.data[[trtvar]], ""),
    COL1 = explicit_na(USUBJID, ""),
    COL2 = paste(AGE, SEX, RACE, sep = concat_sep),
    # Optional Column: COL3/LSVISIT
    COL3 = explicit_na(LSVISIT, ""),
    #COL4 = explicit_na(as.character(EOSDY), ""),
    # Optional Column: COL5/LTVISIT
    COL4 = explicit_na(LTVISIT, ""),
    # Optional Column: COL6/TRTEDY
    #COL5 = explicit_na(as.character(TRTEDY), ""),
    COL5 = ifelse(
      is.na(TRTEDT),
      "",
      toupper(format(as.Date(TRTEDT), format = "%d%b%Y"))
    ),
    # Optional Column: COL7/CUMDOSE/CUMDOSU
    COL6 = paste0(AVAL, " ", AVALU),
    COL7 = ifelse(
      is.na(EOSDT),
      "",
      toupper(format(as.Date(EOSDT), format = "%d%b%Y"))
    ),
    COL8 = case_when(
      DCSREAS == "OTHER" ~ paste0(DCSREAS, " (", stringr::str_to_sentence(DCSREASP), ")"),
      DCSREAS != "OTHER" ~ DCSREAS
    ),
    # Optional Column: COL10/UNBLNDFL
    COL9 = ifelse(is.na(UNBLNDFL), "No", "Yes")
  ) %>%
  arrange(COL0, COL1)

lsting <- lsting |>
  mutate(
    COL5 = ifelse(is.na(TRTEDY), COL5, sprintf("%s (%s)", COL5, TRTEDY)),
    COL7 = ifelse(is.na(EOSDY), COL7, sprintf("%s (%s)", COL7, EOSDY))
  )

lsting <- var_relabel(
  lsting,
  COL0 = "Treatment Group",
  COL1 = "Subject ID",
  COL2 = paste("Age (years)", "Sex", "Race", sep = concat_sep),
  # Optional Column: COL3/LSVISIT
  COL3 = "Last Study Visit~[super a]",
  #COL4 = "Study Day~[super b] of Discontinuation",
  # Optional Column: COL5/LTVISIT
  COL4 = "Last Treatment Visit~[super b]",
  # Optional Column: COL6/TRTEDY
  COL5 = "Date of Last Study Agent Administered (Study Day~[super c])",
  # Optional Column: COL7/CUMDOSE/CUMDOSU
  COL6 = "Cumulative Dose (unit)",
  COL7 = "Date of Discontinuation (Study Day~[super c])",
  COL8 = "Primary Reason for Discontinuation",
  # Optional Column: COL10/UNBLNDFL
  COL9 = "Was Blind Broken?"
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


colwidth <- c(21, 13, 67, 36, 25, 33, 20, 27, 27, 15)

tt_to_tlgrtf(colwidths = colwidth, head(result, 100), file = fileid, orientation = "landscape")
