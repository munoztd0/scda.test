###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lsfae07a.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create lsfae07a: Selected [Narrow / Broad] OCMQs
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adaeocmq
## Output:                    lsfae07a.rtf
## Remarks:                   This variant is for Severity
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
library(rlistings)

###############################################################################
# Define script level parameters
###############################################################################

tblid <- "LSFAE07a"
fileid <- write_path(opath, tblid)
tab_titles <- list(title = "Dummy Title",
                     subtitles = NULL,
                     main_footer = "Dummy Note: On-treatment is defined as ~{optional treatment-emergent}")
popfl <- "SAFFL"
trtvar <- "TRT01A"
ocmqclass <- "Narrow"
key_cols <- c("COL0", "COL1", "COL2", "COL3")
disp_cols <- paste0("COL", 0:8)
concat_sep <- " / "
# Parameter to control whether time should be displayed
include_time <- TRUE
# Parameter to control whether imputed date to be displayed
include_imputed_dates <- FALSE

#Combination treatment parameters

# Provide whether this is a combination treatments
combination_trt <- FALSE

if (combination_trt) {
  comb_trtvars <- c("AEDRGS1", "AEDRGS2") # Provide the variables containing combination treatment information in sequence

  n_comb_trt <- length(comb_trtvars)

  comb_acnvars <- paste0("AEACNS", seq_len(n_comb_trt)) # This will create a variable list like AEACNS1, AEACNS2 etc
}

###############################################################################
# Process data
###############################################################################

adaeocmq <- adaeocmq_jnj %>%
  mutate(
    AEDECOD = case_when(
      AEDECOD == "" ~ paste0("Uncoded: ", AETERM),
      .default = AEDECOD
    )
  ) %>%
  filter(
    !!sym(popfl) == "Y" &
      !is.na(OCMQNAM) &
      OCMQCLSS == ocmqclass &
      TRTEMFL == "Y"
  ) %>%
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
      dplyr::case_when(
        SEX == "F" ~ "Female",
        SEX == "M" ~ "Male",
        SEX == "U" ~ "Unknown",
        SEX == "INTERSEX" ~ "Intersex"
      ),
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
    AEACN = factor(
      case_when(
        AEACN == "DOSE NOT CHANGED" ~ "Dose Not Changed",
        AEACN == "NOT APPLICABLE" ~ "Not Applicable",
        AEACN == "DRUG WITHDRAWN" ~ "Drug Withdrawn",
        AEACN == "DOSE REDUCED" ~ "Dose Reduced",
        AEACN == "DOSE RATE REDUCED" ~ "Dose Rate Reduced",
        AEACN == "DRUG INTERRUPTED" ~ "Drug Interrupted",
        AEACN == "DOSE INCREASED" ~ "Dose Increased",
        AEACN == "UNKNOWN" ~ "Unknown",
        AEACN == "NOT APPLICABLE" ~ "Not Applicable"
      )
    ),
    AESER = case_when(
      AESER == "Y" ~ "Yes",
      AESER == "N" ~ "No",
      .default = NA
    ),
    AEOUT = as.factor(case_when(
      AEOUT == "NOT RECOVERED/NOT RESOLVED" ~ "Not Recovered/Not Resolved",
      AEOUT == "RECOVERED/RESOLVED" ~ "Recovered/Resolved",
      AEOUT == "FATAL" ~ "Fatal",
      TRUE ~ "Other"
    ))
  )

###################################
# Section for combination treatment
###################################

if (combination_trt) {
  adaeocmq <- adaeocmq %>%
    mutate(
      across(
        all_of(comb_acnvars),
        ~ factor(
          case_when(
            . == "DOSE NOT CHANGED" ~ "Dose Not Changed",
            . == "NOT APPLICABLE" ~ "Not Applicable",
            . == "DRUG WITHDRAWN" ~ "Drug Withdrawn",
            . == "DOSE REDUCED" ~ "Dose Reduced",
            . == "DOSE RATE REDUCED" ~ "Dose Rate Reduced",
            . == "DRUG INTERRUPTED" ~ "Drug Interrupted",
            . == "DOSE INCREASED" ~ "Dose Increased",
            . == "UNKNOWN" ~ "Unknown"
          ),
          levels = c(
            "Dose Not Changed",
            "Not Applicable",
            "Drug Withdrawn",
            "Dose Reduced",
            "Dose Rate Reduced",
            "Drug Interrupted",
            "Dose Increased",
            "Unknown"
          )
        )
      )
    )
}

###################################

lsting <- adaeocmq %>%
  mutate(
    AGE = explicit_na(as.character(AGE), ""),
    SEX = explicit_na(SEX, ""),
    RACE = explicit_na(RACE, ""),
    AEDECOD = explicit_na(AEDECOD, ""),
    AETERM = ifelse(
      is.na(AETERM),
      "",
      string_to_title(gsub("\\$", "", AETERM))
    ),
    ASTDT = ifelse(
      nchar(as.character(ASTDT)) == 10,
      toupper(format(ASTDT, "%d%b%Y")),
      ""
    ),
    ASTTM = ifelse(
      include_time & !is.na(ASTDTM),
      substr(as.character(ASTDTM), 12, 16),
      ""
    ),
    ASTDYN = ifelse(!is.na(ASTDY), ASTDY, NA),
    ASTDY = ifelse(!is.na(ASTDY), ASTDY, ""),
    ASTDTFS = ifelse(!is.na(ASTDTF), "*", ""),
    ASTDTFSC = ifelse(!is.na(ASTDTF), paste0(ASTDTF, "*"), ""),
    AENDT = explicit_na(as.character(AENDT), ""),
    AENDT = ifelse(
      nchar(AENDT) == 10,
      toupper(format(as.Date(AENDT), "%d%b%Y")),
      ""
    ),
    AENTM = ifelse(
      include_time & !is.na(AENDTM),
      substr(as.character(AENDTM), 12, 16),
      ""
    ),
    AENDY = ifelse(!is.na(AENDY), AENDY, ""),
    AENDTFS = ifelse(!is.na(AENDTF), "*", ""),
    AENDTFSC = ifelse(!is.na(AENDTF), paste0(AENDTF, "*"), ""),
    AEOUTC = explicit_na(AEOUT, ""),
    AESEV = explicit_na(AESEV, ""),
    AESER = explicit_na(AESER, "")
  )

lsting <- lsting %>%
  {
    if (!combination_trt) {
      mutate(
        .,
        COL7 = explicit_na(AEACN, "")
      )
    } else if (combination_trt) {
      rowwise(.) %>%
        mutate(
          COL7 = paste(
            paste0(c_across(all_of(comb_trtvars)), ": ", c_across(all_of(comb_acnvars))),
            collapse = "; "
          )
        )
    }
  } %>%
  mutate(
    COL0 = OCMQNAM,
    COL1 = explicit_na(.data[[trtvar]], ""),
    COL2 = explicit_na(USUBJID, ""),
    COL3 = paste(AGE, SEX, RACE, sep = concat_sep),
    COL4 = paste(AEDECOD, AETERM, sep = concat_sep),
    COL5 = case_when(
      ASTDT == "" ~ "",
      include_time & include_imputed_dates & ASTDT != "" & ASTTM != "" & ASTDY != "" & ASTDTFSC != "" ~
        paste0(ASTDT, concat_sep, ASTTM, " (", ASTDY, ")", concat_sep, ASTDTFSC),
      include_time & include_imputed_dates & ASTDT != "" & ASTTM != "" & ASTDY != "" & ASTDTFSC == "" ~
        paste0(ASTDT, concat_sep, ASTTM, " (", ASTDY, ")"),
      include_time & include_imputed_dates & ASTDT != "" & ASTTM == "" & ASTDY != "" & ASTDTFSC != "" ~
        paste0(ASTDT, concat_sep, "--:--", " (", ASTDY, ")", concat_sep, ASTDTFSC),
      include_time & include_imputed_dates & ASTDT != "" & ASTTM == "" & ASTDY != "" & ASTDTFSC == "" ~
        paste0(ASTDT, concat_sep, "--:--", " (", ASTDY, ")"),
      include_time & ASTDT != "" & ASTTM != "" & ASTDY != "" ~ paste0(ASTDT, concat_sep, ASTTM, " (", ASTDY, ")"),
      include_time & ASTDT != "" & ASTTM == "" & ASTDY != "" ~ paste0(ASTDT, concat_sep, "--:--", " (", ASTDY, ")"),
      include_imputed_dates & ASTDT != "" & ASTDY != "" & ASTDTFSC != "" ~
        paste0(ASTDT, " (", ASTDY, ")", concat_sep, ASTDTFSC),
      include_imputed_dates & ASTDT != "" & ASTDY != "" & ASTDTFSC == "" ~ paste0(ASTDT, " (", ASTDY, ")"),
      ASTDT != "" & ASTDY != "" ~ paste0(ASTDT, " (", ASTDY, ")"),
    ),
    COL6 = case_when(
      AENDT == "" ~ "",
      include_time & include_imputed_dates & AENDT != "" & AENTM != "" & AENDY != "" & AENDTFSC != "" ~
        paste0(AENDT, concat_sep, AENTM, " (", AENDY, ")", concat_sep, AENDTFSC),
      include_time & include_imputed_dates & AENDT != "" & AENTM != "" & AENDY != "" & AENDTFSC == "" ~
        paste0(AENDT, concat_sep, AENTM, " (", AENDY, ")"),
      include_time & include_imputed_dates & AENDT != "" & AENTM == "" & AENDY != "" & AENDTFSC != "" ~
        paste0(AENDT, concat_sep, "--:--", " (", AENDY, ")", concat_sep, AENDTFSC),
      include_time & include_imputed_dates & AENDT != "" & AENTM == "" & AENDY != "" & AENDTFSC == "" ~
        paste0(AENDT, concat_sep, "--:--", " (", AENDY, ")"),
      include_time & AENDT != "" & AENTM != "" & AENDY != "" ~ paste0(AENDT, concat_sep, AENTM, " (", AENDY, ")"),
      include_time & AENDT != "" & AENTM == "" & AENDY != "" ~ paste0(AENDT, concat_sep, "--:--", " (", AENDY, ")"),
      include_imputed_dates & AENDT != "" & AENDY != "" & AENDTFSC != "" ~
        paste0(AENDT, " (", AENDY, ")", concat_sep, AENDTFSC),
      include_imputed_dates & AENDT != "" & AENDY != "" & AENDTFSC == "" ~ paste0(AENDT, " (", AENDY, ")"),
      AENDT != "" & AENDY != "" ~ paste0(AENDT, " (", AENDY, ")"),
    ),
    COL8 = paste(AEOUTC, AESEV, AESER, sep = concat_sep)
  ) %>%
  arrange(
    COL0,
    COL1,
    COL2,
    !is.na(ASTDYN),
    ASTDYN,
    if (include_time) ASTDTM else ASTDT,
    AEDECOD,
    AETERM
  )

lsting <- var_relabel(
  lsting,
  COL0 = "OCMQ",
  COL1 = "Treatment Group",
  COL2 = "Subject ID",
  COL3 = paste("Age (years)", "Sex", "Race", sep = concat_sep),
  COL4 = paste("Preferred Term", "Reported Term", sep = concat_sep),
  COL5 = if (include_time) {
    paste("Start Date", "Time (Study Day~[super a])", sep = concat_sep)
  } else {
    "Start Date (Study Day~[super a])"
  },
  COL6 = if (include_time) {
    paste("End Date", "Time (Study Day~[super a])", sep = concat_sep)
  } else {
    "End Date (Study Day~[super a])"
  },
  COL7 = "Action Taken With Study Treatment",
  COL8 = paste("Outcome", "Severity", "Serious", sep = concat_sep)
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
# Add titles and footnotes:
###############################################################################

result <- set_titles(result, tab_titles)

###############################################################################
# Output listing
###############################################################################


colwidth <- c(22, 21, 13, 42, 78, 25, 25, 22, 37)

tt_to_tlgrtf(colwidths = colwidth, head(result, 100), file = fileid, orientation = "landscape")
