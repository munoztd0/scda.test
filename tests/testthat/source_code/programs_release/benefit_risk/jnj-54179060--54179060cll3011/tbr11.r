###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tbr11.r
## R Version:                 4.4.2
## junco Version:             0.1.7
## Short Description:         Program to create tbr11: Q-TWiST Analysis Results
##                            Restricted at [X Years];  Analysis Set (Study jjcs - benefit risk)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADAE, ADTTEEFF
## Output:                    TBR11.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:                     1
## Modified By:
## Reporting effort:          Code Refactoring
## Date:                      2026-09-302025-08-13
## Description:               Removed source(read_path(cl, ...)) lines and added appropriate library imports
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



# Benefit risk ---
source(read_path(cl, "junco_benefit_risk.r"))

# Parameters ----

# Define output ID and file location.
tblid <- "TBR11"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define treatment arms to use, the first one will be the active and the
# difference will be first minus second.
trtlvls <- c("Ibrutinib + Venetoclax", "Chlorambucil + Obinutuzumab")

# Define population flag used.
popfl <- "ITTFL"

# Define AE flags to be used.
aefl <- quote(
  TRTEMFL == "Y" & AETOXGRN %in% c(3, 4) & !(is.na(ASTDT) | is.na(AENDT))
)

# Define PFS parameter to use.
pfsparamcd <- "PFSINV"

# Define OS parameter to use.
osparamcd <- "OS"

# Define time scale.
timeunit <- "Months" # Alternatives: Days, Weeks, Months, Years.

# Define the cutoff date.
cutoff_date <- as.Date("2021-02-26")

# Number of bootstrap samples.
nsamples <- 1500

# Confidence level to use.
conflvl <- 0.95

# Utility scores for health states (between 0 and 1).
# These are used to calculate Q-TWiST as follows:
# Q-TWiST = UTWiST × TWiST + UTOX × TOX + UREL × REL
# The base case is:
# c(utwist = 1, utox = 1, urel = 1)
utility <- c(utwist = 1, utox = 0.2, urel = 0.5)

# Column label to use for RMST estimates.
rmstlabel <- "RMST Value"

# Row labels for the QTWiST portion of the table.
qtwist_labels <- c(
  twist = "TWiST~[super a,b]",
  rel = "REL~[super c]",
  qtwist = "Q-TWiST~[super d]"
)

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Derived formats and methods specifications.
formats <- list(
  est = jjcsformat_xx("xx.xx"),
  se = jjcsformat_xx("xx.xx"),
  ci = jjcsformat_xx("(xx.xx, xx.xx)"),
  est_se = jjcsformat_xx("xx.xx (xx.xx)"),
  pval = jjcsformat_pval_fct(alpha)
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
  droplevels() |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl), RANDDT, DTHDT)

## ADAE ----

adae <- read_sas(read_path(a_in, "adae.sas7bdat")) |>
  filter(!!aefl) |>
  select(USUBJID, ASTDT, AENDT, AEDECOD, TRTEDT) |>
  arrange(USUBJID, ASTDT)

## ADTTEEFF ----

adtteeff <- read_sas(read_path(a_in, "adtteeff.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  select(USUBJID, AVAL, PARAMCD, CNSR, EVNTDESC, ADT) |>
  filter(PARAMCD %in% c(pfsparamcd, osparamcd))

## Analysis ----

# Subset to PFS and OS.
anapfs <- adtteeff |>
  filter(PARAMCD == pfsparamcd) |>
  select(-PARAMCD) |>
  rename(
    PFS_AVAL = AVAL,
    PFS_ADT = ADT,
    PFS_CNSR = CNSR,
    PFS_EVNTDESC = EVNTDESC
  )

anaos <- adtteeff |>
  filter(PARAMCD == osparamcd) |>
  select(-PARAMCD) |>
  rename(
    OS_AVAL = AVAL,
    OS_ADT = ADT,
    OS_CNSR = CNSR,
    OS_EVNTDESC = EVNTDESC
  )

# Merge datasets.
anatox01 <- adsl |>
  left_join(adae, by = "USUBJID") |>
  left_join(anapfs, by = "USUBJID") |>
  left_join(anaos, by = "USUBJID") |>
  filter(is.na(ASTDT) | (ASTDT >= RANDDT & ASTDT <= PFS_ADT)) |>
  # The duration of an AE (AEDURN) is defined as the duration of time between the
  # start and end date (resolution date) of the AE if the AE had resolved by the time
  # of progression or death. If the AE had not been resolved by the time of PD or death,
  # then the date of progression is used in the place of the AE end date.
  mutate(
    # The end date for each AE was the resolution date, disease progression date,
    # death date, or end of follow-up (26Feb 2021), whichever occurred first.
    minenddt = pmin(AENDT, PFS_ADT, cutoff_date, na.rm = TRUE),
    # Set last date of AE and duration to pfs_adt or os_date;
    AENDT = if_else(
      !is.na(ASTDT) & !is.na(AENDT) & (AENDT > minenddt),
      minenddt,
      AENDT
    ),
  ) |>
  # Only consider AEs starting before last treatment dose + 30 days.
  filter(is.na(ASTDT) | ASTDT <= (TRTEDT + 30)) |>
  arrange(USUBJID, ASTDT)

# For each patient we sum all days spent with the AEs of interest.
ana <- anatox01 |>
  group_by(USUBJID) |>
  summarize(
    toxdays = calc_tox_days(ASTDT, AENDT)
  ) |>
  left_join(adsl, by = "USUBJID") |>
  left_join(anaos, by = "USUBJID") |>
  left_join(anapfs, by = "USUBJID") |>
  mutate(
    USUBJID = factor(USUBJID),
    mincnsrdt = pmin(OS_ADT, PFS_ADT, cutoff_date, na.rm = TRUE),
    tox_cnsr = if_else(is.na(toxdays) | toxdays == 0, 1, 0), # Censor subjects without any toxicity days
    tox_dur = if_else(
      tox_cnsr == 0,
      toxdays / timefactor,
      0 # Patients who did not experience the AEs of interest
      # before disease progression were assigned a duration of zero for the TOX state.
    )
  ) |>
  select(
    STUDYID,
    USUBJID,
    all_of(trtvar),
    PFS_AVAL,
    PFS_CNSR,
    OS_AVAL,
    OS_CNSR,
    tox_dur,
    tox_cnsr
  )

# Slightly reformat: event indicator instead of censoring indicator, long format.
ana_long <- ana |>
  mutate(
    TOX_AVAL = tox_dur,
    PFS_EVENT = as.logical(1 - PFS_CNSR),
    OS_EVENT = as.logical(1 - OS_CNSR),
    TOX_EVENT = as.logical(1 - tox_cnsr)
  ) |>
  select(
    STUDYID,
    USUBJID,
    all_of(trtvar),
    PFS_AVAL,
    PFS_EVENT,
    OS_AVAL,
    OS_EVENT,
    TOX_AVAL,
    TOX_EVENT
  ) |>
  tidyr::pivot_longer(
    cols = starts_with(c("PFS_", "OS_", "TOX_")),
    names_to = c("PARAM", ".value"),
    names_pattern = "(.+)_(.+)"
  ) |>
  mutate(
    PARAM = factor(PARAM)
  )

# Layout ----

# Generate once all IDs of the stratified bootstrap samples.
set.seed(346)
all_boot_ids <- replicate(
  n = nsamples,
  boot_ids(ana, arm = trtvar),
  simplify = FALSE
)

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
    split_fun = tbr11_split_fun_fct(
      rmst_label = rmstlabel,
      conf_level = conflvl
    )
  ) |>
  split_rows_by(
    "PARAM",
    split_label = paste0("Restricted mean (", timeunit, ")"),
    label_pos = "visible"
  ) |>
  summarize_row_groups(
    "AVAL",
    cfun = tbr11_rmst_acfun,
    extra_args = list(
      arm = trtvar,
      is_event = "EVENT",
      id = "USUBJID",
      boot_ids = all_boot_ids,
      conf_level = conflvl,
      formats = formats
    )
  ) |>
  analyze(
    "AVAL",
    afun = tbr11_qtwist_acfun,
    var_labels = paste0(
      "Restricted mean health state duration (",
      timeunit,
      ")"
    ),
    show_labels = "visible",
    nested = FALSE,
    extra_args = list(
      arm = trtvar,
      param = "PARAM",
      is_event = "EVENT",
      scores = utility,
      id = "USUBJID",
      labels = qtwist_labels,
      boot_ids = all_boot_ids,
      conf_level = conflvl,
      formats = formats
    )
  )

# Output ----

result <- build_table(lyt, ana_long, alt_counts_df = adsl)

# Post-hoc we can suppress the column count for the Overall column.
colcount_visible(result, c(trtvar, "Overall")) <- FALSE

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
