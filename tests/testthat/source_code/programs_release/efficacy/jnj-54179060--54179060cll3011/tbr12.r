###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tbr12.r
## R Version:                 4.4.2
## junco Version:             0.1.7
## Short Description:         Program to create tbr12: Cumulative Time Spent [Months] at
##                            Each Health State; Analysis Set (Study jjcs - benefit risk)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADAE, ADTTEEFF
## Output:                    TBR12.rtf
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
library(survival)
library(junco)
library(haven)



source(read_path(cl, "junco_benefit_risk.r"))


# Parameters ----

# Define output ID and file location.
tblid <- "TBR12"
fileid <- write_path(opath, tblid)

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define treatment arms to use, the first one will be the active and the
# difference will be first minus second.
trtlvls <- c("Ibrutinib + Venetoclax", "Chlorambucil + Obinutuzumab")

# Define population flag used.
popfl <- "ITTFL"

# Define AE filtering flags to be used.
aefl <- quote(
  TRTEMFL == "Y" &
    AETOXGRN %in% c(3, 4) &
    !(is.na(ASTDT) | is.na(AENDT)) &
    ASTDT <= (TRTEDT + 30)
)

# Define PFS parameters to use.
pfsparamcd <- "PFSINV"
pfsevent <- "PD"

# Define time scale.
timeunit <- "Months" # Alternatives: Days, Weeks, Months, Years.

# Specify the length of the time intervals (on the time scale, e.g. months).
timeinterval <- 3

# Define the cutoff date.
cutoff_date <- as.Date("2024-10-31")

# Define the category row labels.
labels <- c(
  "noaepd" = "No AEs/ No PD",
  "aenopd" = "AEs/ No PD",
  "pd" = "Progressive Disease",
  "discfu" = "Discontinued",
  "dead" = "Death"
)

# Confidence level for the CIs.
conflvl <- 0.95

# Formats to use.
formats <- list(
  mean = jjcsformat_xx("xx.xx"),
  mean_ci = jjcsformat_xx("(xx.xx, xx.xx)")
)

# Automatically derived time factor:

timefactor <- switch(
  timeunit,
  Days = 1,
  Weeks = 7,
  Months = 365.25 / 12,
  Years = 365.25,
  stop(paste("time unit", timeunit, "is not implemented"))
)

# Data ----

## ADSL ----

adsl <- read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y", !!rlang::sym(trtvar) %in% trtlvls) |>
  mutate(!!trtvar := factor(.data[[trtvar]], levels = trtlvls)) |>
  select(
    STUDYID,
    USUBJID,
    all_of(trtvar),
    all_of(popfl),
    RANDDT,
    DTHDT,
    EOSDT,
    EOSSTT
  )

## ADAE ----

adae <- read_sas(read_path(a_in, "adae.sas7bdat")) |>
  filter(!!aefl) |>
  select(USUBJID, ASTDT, AENDT, AEDECOD, TRTEDT) |>
  arrange(USUBJID, ASTDT)

## ADTTEEFF ----

adtteeff <- read_sas(read_path(a_in, "adtteeff.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  select(USUBJID, AVAL, PARAMCD, CNSR, EVNTDESC, ADT) |>
  filter(PARAMCD == pfsparamcd)

## Analysis ----

### Step 1: Time grid ----

ana_dur <- adsl |>
  mutate(
    TRTEND = pmin(DTHDT, EOSDT, cutoff_date, na.rm = TRUE),
    TRTDUR = as.numeric(TRTEND - RANDDT) / timefactor
  )
max_trtdur <- max(ana_dur$TRTDUR)
endtime <- ceiling(max_trtdur / timeinterval) * timeinterval
timegrid <- seq(from = 0, to = endtime, by = timeinterval)
intervals <- seq_len(length(timegrid) - 1)

### Step 2: Randomization, Death, PD and Discontinued events ----

ana_rand_death <- adsl |>
  select(USUBJID, RANDDT, DTHDT, all_of(trtvar))

ana_pd <- adtteeff |>
  filter(CNSR == 0, EVNTDESC == pfsevent) |>
  select(USUBJID, ADT) |>
  rename(PDDT = ADT)

ana_disc <- adsl |>
  filter(!is.na(EOSDT) & EOSSTT == "DISCONTINUED") |>
  select(USUBJID, EOSDT) |>
  rename(DISCDT = EOSDT)

### Step 3: AE days ----

ana_ae <- adsl |>
  left_join(adae, by = "USUBJID") |>
  as_tibble() |>
  select(STUDYID, USUBJID, ASTDT, AENDT) |>
  group_by(STUDYID, USUBJID) |>
  summarize(AEDAYS = list(all_days(ASTDT, AENDT))) |>
  ungroup()

### Step 4: Classify status per interval ----

ana_status <- ana_ae |>
  left_join(ana_rand_death, by = "USUBJID") |>
  left_join(ana_pd, by = "USUBJID") |>
  left_join(ana_disc, by = "USUBJID") |>
  rowwise(everything()) |>
  reframe(
    create_health_status(
      RANDDT,
      DTHDT,
      DISCDT,
      PDDT,
      AEDAYS,
      intervals = intervals,
      timegrid = timegrid,
      timefactor = timefactor
    )
  )

### Step 5: Count time per individual per health state ----

ana <- ana_status |>
  select(STUDYID, USUBJID, all_of(trtvar), STATUS) |>
  tidyr::unnest(cols = STATUS) |>
  mutate(STATUS = factor(STATUS, levels = names(labels))) |>
  group_by(STUDYID, USUBJID, across(all_of(trtvar))) |>
  reframe(COUNT = as.integer(table(STATUS)), CATEGORY = labels) |>
  mutate(TIME = COUNT * timeinterval)

# Layout ----

lyt <- basic_table() |>
  split_cols_by(
    trtvar,
    show_colcounts = TRUE,
    split_fun = add_overall_level(
      "Overall", # Note: It needs to stay "Overall" here (this won't be shown anyway).
      label = paste("Difference", paste(trtlvls, collapse = " - ")),
      first = FALSE
    )
  ) |>
  split_cols_by(
    "STUDYID", # The particular choice of this is irrelevant.
    split_fun = tbr12_split_fun_fct(conf_level = conflvl)
  ) |>
  split_rows_by(
    "CATEGORY",
    split_label = "Health State",
    label_pos = "topleft"
  ) |>
  summarize_row_groups(
    "TIME",
    cfun = tbr12_time_acfun,
    extra_args = list(
      arm = trtvar,
      conf_level = conflvl,
      formats = formats
    )
  )

# Output ----

result <- build_table(lyt, ana, alt_counts_df = adsl)

# Post-hoc we can suppress the column count for the Overall column.
colcount_visible(result, c(trtvar, "Overall")) <- FALSE

# Please select titles and footnotes as appropriate.
tab_titles <- get_titles_internal(tblid)

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
