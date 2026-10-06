###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefeq5d01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefeq5d01: EQ-5D-5L Individual Items, Health Status Index,
##                            and EQ-VAS: Means and Mean Changes From Baseline Over Time
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADEQ5D
## Output:                    TEFEQ5D01.rtf
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
# - Define population flag used (default=FAS2FL)
# - Define control treatment arm
################################################################################

tblid <- "TEFEQ5D01"
fileid <- write_path(opath, tblid)

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

trtvar <- "TRT01P"
popfl <- "FAS2FL"

ctrl_grp <- "Treatment B"

paramcdlist <- c("EQ5D0206", "EQ5DHSI", "EQ5DSUM")

################################################################################
# Read in ADSL and ADMADRS datasets
################################################################################
adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y" & !!rlang::sym(trtvar) != "") |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl))

adeq5d <- haven::read_sas(read_path(a_in, "adeq5d.sas7bdat")) |>
  filter(APHASEN == 1 & ANL02FL == "Y" & PARAMCD %in% paramcdlist) |>
  mutate(
    PARAM = as.factor(PARAM),
    AVISIT = as.factor(AVISIT)
  ) |>
  # adapt below as needed
  mutate(
    AVISIT = factor(
      AVISIT,
      levels = c("Baseline (DB)", "Day 43", "Endpoint (DB)")
    )
  ) |>
  select(USUBJID, PARAMCD, AVISIT, APHASEN, PARAM, AVAL, BASE, CHG)

# join data together
eq5d <- adeq5d |> inner_join(adsl, by = c("USUBJID"))

# Set up dataset to use for precision for each parameter
# dp=number of dps used for min/max
# mean will be dp+1
# sd will be dp+2
dp <- data.frame(
  PARAM = c(
    "EQ-5D Health Status Index",
    "EQ5D02-EQ VAS Score",
    "EQ-5D Sum Score"
  ),
  dp = c(2, 0, 0)
)

eq5d <- eq5d |> left_join(dp, by = c("PARAM"))

# adapt below as needed
eq5d$PARAM <- factor(
  eq5d$PARAM,
  levels = c(
    "EQ-5D Health Status Index",
    "EQ5D02-EQ VAS Score",
    "EQ-5D Sum Score"
  )
)
eq5d$PARAM <- forcats::fct_recode(
  eq5d$PARAM,
  "Health Status Index" = "EQ-5D Health Status Index",
  "EQ-VAS Score" = "EQ5D02-EQ VAS Score",
  "Sum Score" = "EQ-5D Sum Score"
)

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
  split_rows_by("PARAM", section_div = " ") |>
  split_rows_by("TRT01P", section_div = " ") |>
  analyze("AVISIT", afun = column_stats(exclude_visits = c("Baseline (DB)")))

result <- build_table(lyt, eq5d)

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

# adjust the column-width for first column and Base Mean (SD) columns to fit all on landscape
# fontspec <- font_spec("Times", 9L, 1.2)
# col_gap <- 7L
# label_width_ins <- 2
#
# colwidths = def_colwidths(result, fontspec,
#                           col_gap = col_gap,
#                           label_width_ins = label_width_ins
# )
#
# acolwidths <- colwidths
# acolwidths[1] <- 55
# acolwidths[8] <- 25

# adjust the column-widths if you want to try and fit onto 1 page landscape
colwidths <- c(40, 8, 14, 14, 14, 12, 12, 24, 8, 14, 14, 14, 14, 12, 12)

################################################################################
# Convert to tbl file and output table
################################################################################
tt_to_tlgrtf(string_map = string_map, tt = 
  result,
  file = fileid,
  orientation = "landscape",
  colwidths = colwidths
)
