###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefmadsg01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefmadsg01:
##                            [Montgomery-Asberg Depression Rating Scale (MADRS) Total Score]
##                            [(Analysis Under [Estimand X])] - Means and Mean Changes From Baseline[(DB)]
##                            Over Time by [Subgroup][ - Double-blind Phase]; Full Analysis Set Analysis Set
##                            (Study mdd)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADMADRSI
## Output:                    TEFMADSG01.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:                     1
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
## Description:
################################################################################

# Environment ----

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)


library(dplyr)
library(rtables)
library(junco)
library(haven)

# Parameters ----

# Define output ID and file location.
tblid <- "TEFMADSG01"
fileid <- write_path(opath, tblid)

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Placebo"

# Define population flags used.
popfl <- "FAS1FL"

# Define subgroup variable used, alongside the label to be used and the order of the
# categories.
subgroup <- "AGEGR1"
subgrlbl <- "Age Group"
subgrorder <- c("18-34 years", "35-54 years", "55-64 years", ">=65 years")

# Define analysis domain flags to be used (for ADMADRSI).
domfl <- quote(ANL03FL == "Y" & APHASEN == 1)

# Define visit flags to be used.
visitfl <- quote(!(AVISIT %in% c("Endpoint (DBFU)")))

# Define response parameter to be used (for ADMADRSI).
resppar <- "MADES1"

# Double blind labeling in the output or not.
doubleblind <- TRUE

# Choose decimal places:
dp <- 0 # Number of decimal places used for min/max.
# Mean will use dp+1 decimal places.
# SD will use dp+2 decimal places.

# Derived formats specifications.
dblabel <- if (doubleblind) " (DB)" else ""
set_default_na_str("NE")

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  mutate(
    !!trtvar := forcats::fct_relevel(!!rlang::sym(trtvar), ctrlab, after = Inf)
  ) |>
  mutate(!!subgroup := factor(.data[[subgroup]], levels = subgrorder)) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(subgroup))

## ADMADRSI ----

admadrsi <- haven::read_sas(read_path(a_in, "admadrsi.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(!!domfl, !!visitfl) |>
  filter(PARAMCD == resppar) |>
  select(USUBJID, AVISIT, AVAL, BASE, CHG) |>
  mutate(
    # Remove the possible (DB) suffix here, such that it won't be duplicate below.
    AVISIT = gsub(" (DB)", "", x = AVISIT, fixed = TRUE),
    USUBJID = factor(USUBJID)
  )

## Analysis ----

ana <- admadrsi |>
  inner_join(adsl, by = "USUBJID") |>
  filter(!!visitfl) |>
  mutate(AVISIT = factor(paste0(AVISIT, dblabel)))

# Layout ----

ana$dp <- dp
mysplitfun <- make_split_fun(post = list(junco:::postfun_eq5d))

lyt <- basic_table() |>
  split_cols_by_multivar(
    c("AVAL", "BASE", "CHG"),
    varlabels = c("Visit Values", " ", "Change from Baseline (DB)")
  ) |>
  split_cols_by("STUDYID", split_fun = mysplitfun) |>
  split_rows_by(
    subgroup,
    page_by = TRUE,
    split_fun = drop_split_levels,
    split_label = subgrlbl,
    section_div = " "
  ) |>
  # Note: This line is currently needed to keep the above page split labels in the rtf,
  # as a workaround for an issue in rtables.
  summarize_row_groups(cfun = junco:::ac_blank_line) |>
  split_rows_by("TRT01P", section_div = " ") |>
  analyze("AVISIT", afun = column_stats(exclude_visits = c("Baseline (DB)")))

# Output ----

result <- build_table(lyt, ana)

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
