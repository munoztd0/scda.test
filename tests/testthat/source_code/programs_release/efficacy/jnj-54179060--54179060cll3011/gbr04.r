###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              gbr04.r
## R Version:                 4.4.2
## junco Version:             0.1.7
## Short Description:         Program to create gbr04: Individual Response Profiles of Subjects
##                            Treated With [Treatment Group A] Over Time (Swimlane Analysis);
##                            [Analysis Set] (Study jjcs - benefit risk)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADAE, ADTTEEFF
## Output:                    GBR04.rtf
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



# Benefit risk ---
source(read_path(cl, "junco_benefit_risk.r"))

# Parameters ----

# Define output ID and file location.
tblid <- "GBR04"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

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

# Optional time origin specification for axis label.
timeorigin <- "since randomization" # Use NULL for omitting origin specification.

# Specify the length of the time intervals (on the time scale, e.g. months).
timeinterval <- 3

# Optional dashed lines time points (in the time unit defined above).
timedash <- c(6, 15) # use c() for omitting the dashed lines

# Define the cutoff date.
cutoff_date <- as.Date("2024-10-31")

# Scores for sorting the patients (higher numbers mean higher up in the plot).
scores <- c(
  "noaepd" = 1,
  "aenopd" = 2,
  "pd" = 3,
  "discfu" = 4,
  "dead" = 5
)
labels <- c(
  "noaepd" = "No AEs/ No PD",
  "aenopd" = "AEs/ No PD",
  "pd" = "Progressive Disease",
  "discfu" = "Discontinued",
  "dead" = "Death"
)
colors <- c(
  "noaepd" = "blue",
  "aenopd" = "orange",
  "pd" = "red",
  "discfu" = "grey",
  "dead" = "black"
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
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
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
  select(USUBJID, ASTDT, AENDT) |>
  group_by(USUBJID) |>
  summarize(AEDAYS = list(all_days(ASTDT, AENDT)))

### Step 4: Classify status per interval ----

ana <- ana_ae |>
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
    ) |>
      mutate(
        SCORE = scores[rev(STATUS[[1]])] |> paste(collapse = "")
      )
  )

# Functions ----

plot_swimlane <- function(data) {
  data <- data |>
    arrange(SCORE) |>
    mutate(YORDER = as.numeric(factor(USUBJID, levels = unique(USUBJID)))) |>
    select(USUBJID, INTERVAL, STATUS, YORDER) |>
    unnest(cols = c(INTERVAL, STATUS)) |>
    mutate(
      STATUS = factor(
        STATUS,
        levels = c("noaepd", "aenopd", "pd", "discfu", "dead")
      )
    )

  # Main swimmer plot.
  p1 <- ggplot(
    data,
    aes(
      xmin = (INTERVAL - 1) * timeinterval,
      xmax = INTERVAL * timeinterval,
      ymin = YORDER - 0.5,
      ymax = YORDER + 0.5,
      fill = STATUS
    )
  ) +
    geom_rect() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(
      breaks = timegrid,
      limits = range(timegrid),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = colors,
      labels = labels,
      drop = FALSE,
      name = "Health Status"
    ) +
    theme_minimal() +
    labs(
      x = paste0(
        "Time",
        ifelse(!is.null(timeorigin), paste0(" ", timeorigin), ""),
        " (",
        timeunit,
        ")"
      ),
      y = "Subjects"
    ) +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom"
    )

  # Optional dashed lines.
  if (length(timedash) > 0) {
    p1 <- p1 +
      geom_vline(xintercept = timedash, linetype = "dashed", color = "white")
  }

  # Table: subject count per category per interval.
  tab_df <- data |>
    mutate(LABEL = labels[STATUS]) |>
    group_by(INTERVAL, LABEL) |>
    summarise(N = n_distinct(USUBJID), .groups = "drop") |>
    tidyr::complete(INTERVAL, LABEL = labels, fill = list(N = 0)) |>
    tidyr::pivot_wider(names_from = INTERVAL, values_from = N) |>
    mutate(LABEL = factor(LABEL, levels = labels)) |>
    arrange(desc(LABEL)) |>
    tidyr::pivot_longer(-LABEL, names_to = "INTERVAL", values_to = "Count") |>
    mutate(INTERVAL = as.numeric(INTERVAL))

  p2 <- ggplot(
    tab_df,
    aes(
      x = (INTERVAL - 0.5) * timeinterval,
      y = LABEL,
      label = Count,
      color = LABEL
    )
  ) +
    geom_text(size = 3.5, na.rm = TRUE) +
    theme_minimal() +
    labs(x = NULL, y = NULL) +
    scale_color_manual(
      values = colors[names(labels)] |> stats::setNames(labels),
      guide = "none"
    ) +
    scale_x_continuous(limits = range(timegrid), expand = c(0, 0)) +
    theme(
      axis.text.x = element_blank(),
      panel.grid = element_blank()
    )

  # Combine plot and table.
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

  plots[[png_name]] <- plot_swimlane(ana_arm)
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
    height = 10,
    width = 15,
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
    plotwidth = 15 / 2,
    plotheight = 10 / 2,
    orientation = "landscape",
    opath = write_path(opath),
    file = rtf_name,
    title = title,
    footers = title_footer$main_footer
  )
  index <- index + 1
}
