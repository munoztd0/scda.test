###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lsimh02.r
## R version:                 4.2.1
## junco version:             0.1.3
## Short Description:         Create lsimh02: Listing of Procedures
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, pr, supppr
## Output:                    lsimh02.rtf
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

tblid <- "LSIMH02"
fileid <- write_path(opath, tblid)
popfl <- "FASFL"
trtvar <- "TRT01A"
key_cols <- c("COL0", "COL1", "COL2")
sort_cols <- c("COL0", "COL1", "COL2", "PRSTDTC", "PRDECOD")
disp_cols <- paste0("COL", 0:11)
concat_sep <- " / "
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

###############################################################################
# Process data
###############################################################################

adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) |>
  df_na() |>
  filter(!!rlang::sym(popfl) == "Y") |>
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
      levels = c("Male", "Female")
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

pr <- haven::read_sas(envsetup::read_path(d_in, "pr.sas7bdat")) |>
  df_na()

supppr <- haven::read_sas(envsetup::read_path(d_in, "supppr.sas7bdat")) |>
  select(STUDYID, USUBJID, IDVAR, IDVARVAL, QNAM, QVAL) |>
  pivot_wider(
    id_cols = c(STUDYID, USUBJID, IDVAR, IDVARVAL),
    names_from = QNAM,
    values_from = QVAL
  ) |>
  df_na()


adsl_pr <- adsl |>
  inner_join(pr |> mutate(PRSEQ_CHR = as.character(PRSEQ)), by = c("STUDYID", "USUBJID")) |>
  left_join(
    supppr,
    by = c("STUDYID", "USUBJID", "PRSEQ_CHR" = "IDVARVAL")
  )

###############################################################################
# ISO 8601 duration parsing helper
# P1D -> duration="1", unit="(Days)"
# PT1H -> duration="1", unit="(Hours)"
# PT15M -> duration="15", unit="(Minutes)"
# PT1H30M -> duration="1, 30", unit="(Hours, Minutes)"
###############################################################################

parse_iso_duration <- function(x) {
  if (is.na(x) || x == "") {
    return(list(duration = "", unit = ""))
  }
  parts <- list(
    Days = stringr::str_match(x, "P(\\d+)D")[, 2],
    Hours = stringr::str_match(x, "T(\\d+)H")[, 2],
    Minutes = stringr::str_match(x, "(\\d+)M")[, 2]
  )
  parts <- Filter(Negate(is.na), parts)
  if (length(parts) == 0) {
    return(list(duration = "", unit = ""))
  }
  list(
    duration = paste(unname(parts), collapse = ", "),
    unit = paste0("(", paste(names(parts), collapse = ", "), ")")
  )
}

lsting <- adsl_pr |>
  mutate(
    AGE = explicit_na(as.character(AGE), ""),
    SEX = explicit_na(SEX, ""),
    RACE = explicit_na(RACE, ""),
    # Parse ISO 8601 PRDUR
    PRDUR_PARSED = lapply(PRDUR, parse_iso_duration),
    PRDUR_VAL = sapply(PRDUR_PARSED, `[[`, "duration"),
    PRDUR_UNIT = sapply(PRDUR_PARSED, `[[`, "unit"),
    # Start date components
    PRSTDT = as.Date(PRSTDTC, format = "%Y-%m-%d"),
    PRSTYR = ifelse(stringr::str_length(PRSTDTC) >= 4, substr(PRSTDTC, 1, 4), NA),
    PRSTMO = month.abb[as.numeric(ifelse(
      stringr::str_length(sub("T.*", "", PRSTDTC)) >= 7 &
        substr(sub("T.*", "", PRSTDTC), 6, 7) != "--",
      substr(sub("T.*", "", PRSTDTC), 6, 7),
      NA
    ))],
    PRSTDAY = ifelse(stringr::str_length(PRSTDTC) >= 10, substr(PRSTDTC, 9, 10), NA),
    PRSTDY = ifelse(PRSTDT > TRTSDT, PRSTDT - TRTSDT + 1, PRSTDT - TRTSDT),
    PRSTDYL = ifelse(!is.na(PRSTDT) & is.na(PRSTDY), "-", as.character(PRSTDY)),
    # End date components
    PRENDT = as.Date(PRENDTC, format = "%Y-%m-%d"),
    PRENYR = ifelse(stringr::str_length(PRENDTC) >= 4, substr(PRENDTC, 1, 4), NA),
    PRENMO = month.abb[as.numeric(ifelse(
      stringr::str_length(sub("T.*", "", PRENDTC)) >= 7 &
        substr(sub("T.*", "", PRENDTC), 6, 7) != "--",
      substr(sub("T.*", "", PRENDTC), 6, 7),
      NA
    ))],
    PRENDAY = ifelse(stringr::str_length(PRENDTC) >= 10, substr(PRENDTC, 9, 10), NA),
    PRENDY = ifelse(PRENDT > TRTSDT, PRENDT - TRTSDT + 1, PRENDT - TRTSDT),
    PRENDYL = case_when(
      !is.na(PRENDY) ~ as.character(PRENDY),
      is.na(PRENDY) & !is.na(PRENDT) ~ "-",
      TRUE ~ NA_character_
    )
  ) |>
  unite("PRSTDTL", PRSTDAY, PRSTMO, PRSTYR, sep = "", na.rm = TRUE, remove = FALSE) |>
  unite("PRENDTL", PRENDAY, PRENMO, PRENYR, sep = "", na.rm = TRUE, remove = FALSE) |>
  mutate(
    COL0 = explicit_na(.data[[trtvar]], ""),
    COL1 = explicit_na(USUBJID, ""),
    COL2 = paste(AGE, SEX, RACE, sep = concat_sep),
    COL3 = ifelse(toupper(PREVINTX) == "BEFORE INTO THE STUDY", "Yes", "No"),
    COL4 = paste(
      explicit_na(stringr::str_to_sentence(PRDECOD), ""),
      explicit_na(stringr::str_to_sentence(PRTRT), ""),
      sep = concat_sep
    ),
    COL5 = explicit_na(stringr::str_to_sentence(PRINDC), ""),
    COL6 = case_when(
      PRSTDTL != "" & !is.na(PRSTDYL) ~ paste0(toupper(PRSTDTL), " (", PRSTDYL, ")"),
      PRSTDTL != "" & is.na(PRSTDYL) ~ toupper(PRSTDTL),
      TRUE ~ ""
    ),
    COL7 = case_when(
      PRENDTL != "" & !is.na(PRENDYL) ~ paste0(toupper(PRENDTL), " (", PRENDYL, ")"),
      PRENDTL != "" & is.na(PRENDYL) ~ toupper(PRENDTL),
      TRUE ~ ""
    ),
    COL8 = ifelse(
      PRDUR_VAL != "" & PRDUR_UNIT != "",
      paste0(PRDUR_VAL, " ", PRDUR_UNIT),
      explicit_na(PRDUR_VAL, "")
    ),
    COL9 = ifelse(toupper(PRPLN) == "Y", "Yes", "No"),
    COL10 = explicit_na(stringr::str_to_sentence(PRFIND), ""),
    COL11 = ifelse(toupper(PRAEFIND) == "Y", "Yes", "No")
  ) |>
  arrange(COL0, COL1, COL2, PRSTDTC, PRDECOD)

lsting <- var_relabel(
  lsting,
  COL0 = "Treatment Group",
  COL1 = "Subject ID",
  COL2 = paste("Age (years)", "Sex", "Race", sep = concat_sep),
  COL3 = "Procedure/Surgery Planned Before Study Entry?",
  COL4 = paste("Preferred Term", "Reported Term", sep = concat_sep),
  COL5 = "Indication",
  COL6 = "Start Date/Time (Study Day~[super a])",
  COL7 = "End Date/Time (Study Day~[super a])",
  COL8 = "Procedure Duration (Unit)",
  COL9 = "Procedure Elective?",
  COL10 = "Diagnostic Findings",
  COL11 = "Findings Adverse Event?"
)

###############################################################################
# Build listing
###############################################################################

result <- rlistings::as_listing(
  df = lsting,
  key_cols = key_cols,
  disp_cols = disp_cols,
  sort_cols = sort_cols,
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
