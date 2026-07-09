###############################################################################
# Prep environment
###############################################################################

library(envsetup)
library(tern)
library(dplyr)
library(tidyr)
library(rtables)
library(rlistings)
library(junco)

###############################################################################
# Define script level parameters
###############################################################################

tblid <- "LSIDS06"
fileid <- write_path(opath, tblid)
popfl <- "RANDFL"
trtvar <- "TRT01P"
key_cols <- c("COL0", "COL1", "COL2")
disp_cols <- paste0("COL", 0:4)
tab_titles <- list(title = "Dummy Title",
                     subtitles = NULL,
                     main_footer = "Dummy Note: On-treatment is defined as ~{optional treatment-emergent}")
flag_reason_map <- c(SAFFL = "SAFEXRS", FASFL = "FASEXRS", PPROTFL = "PPREXRS", PKFL = "PKEXRES", IMFL = "IMEXRES")
flag_vars <- names(flag_reason_map)
reason_vars <- unname(flag_reason_map)
concat_sep <- " / "

###############################################################################
# Process data
###############################################################################

adsl <- adsl_jnj %>%
  filter(!!rlang::sym(popfl) == "Y") %>%
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

labels <- sapply(flag_vars, function(v) {
  lbl <- attr(adsl[[v]], "label")
  lbl <- gsub("\\bflag\\b", "", lbl, ignore.case = TRUE) # remove word flag
  lbl <- gsub("\\s{2,}", " ", lbl) # collapse extra spaces
  trimws(lbl) # trim ends
})

adsl_excl_ana <- adsl %>%
  select(USUBJID, .data[[trtvar]], AGE, SEX, RACE, all_of(flag_vars), all_of(reason_vars)) %>%
  pivot_longer(cols = all_of(flag_vars), names_to = "FLAG_VAR", values_to = "FLAG_VAL") %>%
  filter(toupper(FLAG_VAL) != "Y") %>%
  mutate(FLAG_DESC = labels[FLAG_VAR]) %>%
  rowwise() %>%
  mutate(EXCL_REASON = adsl[[flag_reason_map[FLAG_VAR]]][adsl$USUBJID == USUBJID]) %>%
  ungroup() %>%
  select(USUBJID, .data[[trtvar]], AGE, SEX, RACE, FLAG_VAR, FLAG_VAL, FLAG_DESC, EXCL_REASON)


lsting <- adsl_excl_ana %>%
  mutate(
    AGE = explicit_na(as.character(AGE), ""),
    SEX = explicit_na(SEX, ""),
    RACE = explicit_na(RACE, ""),
    COL0 = explicit_na(.data[[trtvar]], ""),
    COL1 = explicit_na(USUBJID, ""),
    COL2 = paste(AGE, SEX, RACE, sep = concat_sep),
    COL3 = explicit_na(FLAG_DESC, ""),
    COL4 = explicit_na(EXCL_REASON, ""),
  )

lsting <- var_relabel(
  lsting,
  COL0 = "Treatment Group",
  COL1 = "Subject ID",
  COL2 = paste("Age (years)", "Sex", "Race", sep = concat_sep),
  COL3 = "Analysis Set",
  COL4 = "Reason(s) for Exclusion"
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


colwidth <- c(18, 63, 67, 53, 80)

tt_to_tlgrtf(colwidths = colwidth, head(result, 100), file = fileid, orientation = "landscape")
