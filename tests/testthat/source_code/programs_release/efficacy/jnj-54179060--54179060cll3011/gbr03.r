###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              gbr03.r
## R Version:                 4.4.2
## junco Version:             0.1.7
## Short Description:         Program to create gbr03: Mean Duration (Unweighted) of Health States;
##                            [Analysis Set] Analysis Set (Study jjcs - benefit risk)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADAE, ADTTEEFF
## Output:                    GBR03.rtf
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
library(ggplot2)
library(patchwork)
library(scales)
library(tidyr)
library(tidytlg)
library(junco)
library(haven)
library(survival)



# Benefit risk ---
source(read_path(cl, "junco_benefit_risk.r"))

# Parameters ----

# Define output ID and file location.
tblid <- "GBR03"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

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

# Define colors for the events, in the order "REL", "TWIST", "TOX".
event_colors <- c("blue", "orange", "darkgray")

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
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
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
    USUBJID,
    all_of(trtvar),
    PFS_AVAL,
    PFS_CNSR,
    OS_AVAL,
    OS_CNSR,
    tox_dur,
    tox_cnsr
  )

# Functions ----

# Function to plot RMST for health states with one bar per treatment arm.
plot_rmst_bars <- function(data) {
  data_by_arm <- split(data, data[[trtvar]])
  rmst_by_arm <- lapply(
    data_by_arm,
    function(df) {
      km_pfs <- survfit(Surv(PFS_AVAL, 1 - PFS_CNSR) ~ 1, data = df)
      km_os <- survfit(Surv(OS_AVAL, 1 - OS_CNSR) ~ 1, data = df)
      km_ae <- survfit(Surv(tox_dur, 1 - tox_cnsr) ~ 1, data = df)

      rmst_pfs <- summary(km_pfs)$table["rmean"]
      rmst_os <- summary(km_os)$table["rmean"]
      rmst_ae <- summary(km_ae)$table["rmean"]

      c(PFS = unname(rmst_pfs), OS = unname(rmst_os), TOX = unname(rmst_ae))
    }
  )
  rmst <- do.call(rbind, rmst_by_arm) |>
    as.data.frame() |>
    dplyr::mutate(
      TWIST = PFS - TOX,
      REL = OS - PFS
    ) |>
    select(TOX, TWIST, REL) |>
    as_tibble(rownames = "Treatment") |>
    pivot_longer(
      cols = c(TOX, TWIST, REL),
      names_to = "HealthState",
      values_to = "Value"
    ) |>
    mutate(HealthState = factor(HealthState, levels = c("REL", "TWIST", "TOX")))

  # Create the stacked bar plot
  p <- ggplot(rmst, aes(x = Treatment, y = Value, fill = HealthState)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(
      # Add value labels
      aes(label = sprintf("%.2f", Value)),
      position = position_stack(vjust = 0.5),
      color = "white"
    ) +
    labs(
      x = "Treatment",
      y = paste0("Mean Duration (", timeunit, ")")
    ) +
    scale_fill_manual(values = event_colors) +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Output ----

title_footer <- tab_titles

# Produce the plot, save it as png, export to rtf.
png_name <- paste0(tolower(tblid), ".png")
barplot <- plot_rmst_bars(ana)
png_file <- write_path(opath, png_name)
ggsave(
  png_file,
  barplot,
  height = 5,
  device = png,
  type = "cairo",
  dpi = 300
)
print(png_file)

# Export to rtf.
rtf_name <- tblid
tidytlg::gentlg(
  tlf = "g",
  plotnames = png_file,
  plotheight = 5,
  orientation = "landscape",
  opath = write_path(opath),
  file = rtf_name,
  title = title_footer$title,
  footers = title_footer$main_footer
)
