###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefmad07.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefmad07: MADRS Total Score [(Analysis under [Estimand X])]: Means
##                            and Mean Changes From Baseline[(DB)] Over Time
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADMADRS
## Output:                    TEFMAD07.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:                     1
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
## Description:
################################################################################

################################################################################
# Prep Environment
################################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)


library(dplyr)
library(rtables)
library(junco)
library(haven)

################################################################################
# Define script level parameters:
################################################################################

################################################################################
# - Define output ID and file location
# - Define treatment variable used (default=TRT01P)
# - Define population flag used (default=FAS1FL)
# - Define control treatment arm
# - Define parameter code PARAMCD needed for output
################################################################################

tblid <- "TEFMAD07"
fileid <- write_path(opath, tblid)

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

trtvar <- "TRT01P"
popfl <- "FAS1FL"

ctrl_grp <- "Treatment B"

param_code <- "MADES1"

################################################################################
# Read in ADSL and ADMADRS datasets
################################################################################
adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y" & !!rlang::sym(trtvar) != "") |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl))

admadrsi <- haven::read_sas(read_path(a_in, "admadrsi.sas7bdat")) |>
  filter(ANL03FL == "Y" & PARAMCD == param_code & APHASEN == 1) |>
  mutate(AVISIT = as.factor(AVISIT)) |>
  # adapt below as needed
  mutate(
    AVISIT = factor(
      AVISIT,
      levels = c(
        "Baseline (DB)",
        "Day 15",
        "Day 29",
        "Day 43",
        "Endpoint (DBFU)"
      )
    )
  ) |>
  select(USUBJID, AVISIT, APHASEN, AVAL, BASE, CHG)

# adapt below as needed
admadrsi$AVISIT <- forcats::fct_recode(
  admadrsi$AVISIT,
  "Baseline (DB)" = "Baseline (DB)",
  "Day 15" = "Day 15",
  "Day 29" = "Day 29",
  "Day 43" = "Day 43",
  "Endpoint" = "Endpoint (DBFU)"
)

# join data together
madrs <- admadrsi |> inner_join(adsl, by = c("USUBJID"))

# Choose a dp
# dp=number of dps used for min/max
# mean will be dp+1
# sd will be dp+2
madrs$dp <- 0

################################################################################
# Define layout and build table:
################################################################################

mysplitfun <- make_split_fun(post = list(junco:::postfun_eq5d))

lyt <- basic_table() |>
  split_cols_by_multivar(
    c("AVAL", "BASE", "CHG"),
    varlabels = c("Visit Values", " ", "Change from Baseline (DB)")
  ) |>
  split_cols_by("STUDYID", split_fun = mysplitfun) |>
  split_rows_by("TRT01P", section_div = " ") |>
  analyze("AVISIT", afun = column_stats(exclude_visits = c("Baseline (DB)")))

result <- build_table(lyt, madrs)

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

# adjust the column-width for first column and Base Mean (SD) columns to fit all on landscape
fontspec <- font_spec("Times", 9L, 1.2)
col_gap <- 7L
label_width_ins <- 2

colwidths <- def_colwidths(
  result,
  fontspec,
  col_gap = col_gap,
  label_width_ins = label_width_ins
)

acolwidths <- colwidths
acolwidths[1] <- 55
acolwidths[8] <- 18

################################################################################
# Convert to tbl file and output table
################################################################################
tt_to_tlgrtf(string_map = string_map, tt = 
  result,
  file = fileid,
  orientation = "landscape",
  colwidths = acolwidths
)
