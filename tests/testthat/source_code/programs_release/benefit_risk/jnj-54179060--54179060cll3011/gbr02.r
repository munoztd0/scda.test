###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              gbr02.r
## R Version:                 4.4.2
## junco Version:             0.1.7
## Short Description:         Program to create gbr02: Partitioned Survival Curves for
##                            [Treatment Group A] Based on the Q-TWIST Analysis;
##                            [Analysis Set] Analysis Set (Study jjcs - benefit risk)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADAE, ADTTEEFF
## Output:                    GBR02.rtf
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
tblid <- "GBR02"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define population flag used.
popfl <- "ITTFL"

# Define AE filtering flags to be used.
aefl <- quote(
  TRTEMFL == "Y" & AETOXGRN %in% c(3, 4) & !(is.na(ASTDT) | is.na(AENDT))
)

# Define PFS parameter to use.
pfsparamcd <- "PFSINV"

# Define OS parameter to use.
osparamcd <- "OS"

# Define time scale.
timeunit <- "Months" # Alternatives: Days, Weeks, Months, Years.
starttime <- "Date of Randomization"

# Specify the desired time points for subjects at risk table.
timegrid <- seq(0, 78, by = 6)

# Define the cutoff date.
cutoff_date <- as.Date("2021-02-26")

# Colors to use, in the order "OS", "PFS", and "TOX".
event_colors <- c("blue", "orange", "darkgray")

# Line types to use, in the same order as colors.
event_linetypes <- c("solid", "dashed", "dotted")

# Whether to fill the areas with color between the survival curves.
fill_areas <- TRUE

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
    # death date, or end of follow-up (cut off date), whichever occurred first.
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

# Function to plot survival curves for a treatment arm.
plot_survival_curves <- function(data) {
  km_pfs <- survfit(Surv(PFS_AVAL, 1 - PFS_CNSR) ~ 1, data = data)
  km_os <- survfit(Surv(OS_AVAL, 1 - OS_CNSR) ~ 1, data = data)
  km_ae <- survfit(Surv(tox_dur, 1 - tox_cnsr) ~ 1, data = data)

  km_data <- bind_rows(
    data.frame(time = km_ae$time, surv = km_ae$surv * 100, event = "TOX"), # Multiply by 100
    data.frame(time = km_pfs$time, surv = km_pfs$surv * 100, event = "PFS"), # Multiply by 100
    data.frame(time = km_os$time, surv = km_os$surv * 100, event = "OS") # Multiply by 100
  )

  # Create a data frame for subjects at risk after each event comparison.
  generate_risk_data <- function(km_fit, event_name) {
    this_stepfun <- stepfun(
      x = km_fit$time,
      y = c(km_fit$n.risk, km_fit$n.risk[length(km_fit$n.risk)])
    )
    risk_df <- data.frame(
      Time = timegrid,
      SubjectsAtRisk = this_stepfun(timegrid)
    )
    risk_df$Event <- event_name
    risk_df
  }

  risk_data <- bind_rows(
    generate_risk_data(km_ae, "TOX"),
    generate_risk_data(km_pfs, "PFS"),
    generate_risk_data(km_os, "OS")
  )

  # Creating a summary table for subjects at risk for specific time points
  risk_summary <- risk_data %>%
    group_by(Time) %>%
    summarize(
      TOX = first(SubjectsAtRisk[Event == "TOX"]),
      PFS = first(SubjectsAtRisk[Event == "PFS"]),
      OS = first(SubjectsAtRisk[Event == "OS"])
    ) %>%
    replace_na(list(TOX = 0, PFS = 0, OS = 0)) %>% # Replace NA with 0
    pivot_longer(
      cols = c(TOX, PFS, OS),
      names_to = "Event",
      values_to = "SubjectsAtRisk"
    )

  # Create the survival plot
  p1 <- ggplot(km_data) +
    geom_step(
      aes(x = time, y = surv, color = event, linetype = event),
      linewidth = 1
    ) +
    labs(
      x = paste(timeunit, "from", starttime),
      y = "% of Subjects Without Event"
    ) +
    scale_color_manual(values = event_colors) +
    scale_linetype_manual(values = event_linetypes) +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
    ) +
    scale_y_continuous(
      breaks = seq(from = 0, to = 100, by = 10),
      limits = c(0, 100)
    ) +
    scale_x_continuous(
      breaks = timegrid,
      limits = c(min(timegrid) - 1, max(timegrid) + 1)
    )

  if (fill_areas) {
    # Create the data to fill the area in between the curves.
    km_ribbon_data <- gen_step_ribbon_data(
      df = km_data,
      curve_type = "event",
      x = "time",
      y = "surv",
      baseline_y = 0,
      start_x = 0,
      curve_diff = list(
        OS = c("PFS", "OS"),
        PFS = c("TOX", "PFS"),
        TOX = c("baseline", "TOX")
      )
    )

    # Add this to the plot.
    p1 <- p1 +
      geom_ribbon(
        data = km_ribbon_data,
        aes(
          x = time,
          ymin = y_min,
          ymax = y_max,
          group = curve_diff,
          fill = curve_diff
        ),
        alpha = 0.5
      ) +
      scale_fill_manual(values = event_colors, guide = "none")
  }

  p2 <- ggplot(
    risk_summary,
    aes(x = Time, y = Event, label = SubjectsAtRisk, color = Event)
  ) +
    scale_color_manual(values = event_colors, guide = "none") +
    geom_text(size = 3, na.rm = TRUE) +
    theme_minimal() +
    labs(x = NULL, y = "Subjects at Risk") +
    scale_x_continuous(breaks = timegrid) +
    theme(
      axis.text.y = element_text(size = 7),
      axis.text.x = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 1, hjust = 1),
      panel.grid = element_blank()
    )

  # Combine plot and subjects at risk table without row numbers and column values
  patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(3, 1))
}

# Output ----

title_footer <- tab_titles

# Loop through each treatment arm, produce the plot, save it as png, export to rtf.
plots <- list()
index <- 1
for (arm in levels(ana[[trtvar]])) {
  arm_name <- gsub("+", "Plus", x = arm, fixed = TRUE)
  arm_name <- gsub(" ", "", x = arm_name)
  arm_name <- make.names(arm_name)

  ana_arm <- filter(ana, !!rlang::sym(trtvar) == arm)
  png_name <- paste0(tolower(tblid), "_", arm_name, ".png")

  plots[[png_name]] <- plot_survival_curves(ana_arm)
  title <- gsub(
    "[Treatment Group A]",
    arm,
    x = title_footer$title,
    fixed = TRUE
  )

  png_file <- write_path(opath, png_name)
  ggsave(
    png_file,
    plots[[png_name]],
    width = 7,
    height = 5,
    device = png,
    type = "cairo",
    dpi = 300
  )
  print(png_file)

  # Export to rtf.
  rtf_name <- paste0(tblid, letters[index])
  tidytlg::gentlg(
    tlf = "g",
    plotnames = png_file,
    plotheight = 5,
    orientation = "landscape",
    opath = write_path(opath),
    file = rtf_name,
    title = title,
    footers = title_footer$main_footer
  )
  index <- index + 1
}
