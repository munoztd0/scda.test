###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              gefpfs02.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Forest Plot of Subgroup Analysis - HR
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADTTEEF
## Output:                    GEFPFS02.png
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
## Description:
################################################################################

################################################################################
# Prep environment:
################################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)


library(patchwork)

################################################################################
# Define output ID:
################################################################################

tblid <- "GEFPFS02"

################################################################################
# Get titles and footnotes:
################################################################################

# title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

adtte <- haven::read_sas(read_path(a_in, "adtteef.sas7bdat")) |>
  filter(ITTFL == "Y" & PARAMCD == "IRCPFSM") |>
  mutate(
    trt = factor(
      case_when(
        TRT01P == "Dummy A - Tec-Dara" ~ "Dummy A",
        TRT01P == "Dummy B - DPd" ~ "Dummy B",
        TRT01P == "Dummy B - DVd" ~ "Dummy B"
      ),
      levels = c("Dummy B", "Dummy A")
    ),
    is_event = CNSR == 0,
    cnsrc_var = factor(
      case_when(
        CNSR == 0 ~ "Event",
        CNSR == 1 ~ "Censored"
      ),
      levels = c("Event", "Censored")
    )
  ) |>
  mutate(
    AGEGR2 = factor(AGEGR2, levels = c("<65 years", ">=65 years")),
    SEX = factor(SEX, levels = c("M", "F"), labels = c("Male", "Female")),
    RACEGR2 = factor(RACEGR2, levels = c("White", "Asian", "Other")),
    RENGRP1 = factor(RENGRP1, levels = c("<=60 mL/min", ">60 mL/min")),
    ECOGGR1 = factor(ECOGGR1, levels = c("0", ">=1")),
    STRAT01 = factor(STRAT01, levels = c("DPd", "DVd")),
    PRLNGRP = factor(
      PRLNGRP,
      levels = c("1 line", "2 or 3 lines"),
      labels = c("1", "2 or 3")
    ),
    ISSGRP = factor(ISSGRP, levels = c("I", "II", "III")),
    PLASGRP = factor(PLASGRP, levels = c("0", ">=1"), labels = c("No", "Yes")),
    BSLRISK = factor(
      BSLRISK,
      levels = c("High risk", "Standard risk"),
      labels = c("High-risk", "Standard-risk")
    ),
    PLSCGRP = factor(
      PLSCGRP,
      levels = c("<= 30", "> 30 - < 60", ">= 60"),
      labels = c("<= 30", "> 30 to < 60", ">= 60")
    )
  )

# prepare data frame for survival subgroup analysis ---------------

## define subgroup variables
subgrp <- c(
  "AGEGR2",
  "SEX",
  "RACEGR2",
  "RENGRP1",
  "ECOGGR1",
  "STRAT01",
  "PRLNGRP",
  "ISSGRP",
  "PLASGRP",
  "BSLRISK",
  "PLSCGRP"
)

## generate events, median, hazard ration and 95%CI from tern
df_grouped <- extract_survival_subgroups(
  variables = list(
    tte = "AVAL",
    is_event = "is_event",
    arm = "trt",
    subgroups = subgrp
  ),
  data = adtte
)

## format columns in output presentation
surv_df <- df_grouped$survtime
hr_df <- df_grouped$hr

surv_df2 <- surv_df |>
  mutate(
    arm = case_when(
      arm == "Dummy A" ~ "col1",
      arm == "Dummy B" ~ "col2"
    )
  ) |>
  pivot_wider(
    names_from = arm,
    names_glue = "{arm}_{.value}",
    values_from = c(n, n_events, median)
  )

df <- surv_df2 |>
  inner_join(hr_df, by = c("subgroup", "var", "var_label", "row_type")) |>
  select(-arm) |>
  mutate(
    # rounding in SAS rules
    across(
      c(hr, lcl, ucl),
      ~ tidytlg::roundSAS(.x, digits = 2, as_char = TRUE, na_char = "NE")
    ),
    across(
      c(col1_median, col2_median),
      ~ tidytlg::roundSAS(.x, digits = 1, as_char = TRUE, na_char = "NE")
    ),
    # formatting as table presentation
    hr_ci = paste0(hr, " (", lcl, ", ", ucl, ")"),
    col1_evt = paste0(col1_n_events, "/", col1_n),
    col2_evt = paste0(col2_n_events, "/", col2_n)
  ) |>
  select(
    subgroup,
    var,
    var_label,
    col1_evt,
    col1_median,
    col2_evt,
    col2_median,
    hr_ci,
    hr,
    lcl,
    ucl
  )

df$hr <- as.numeric(df$hr)
df$lcl <- as.numeric(df$lcl)
df$ucl <- as.numeric(df$ucl)

## insert rows for each subgroup variable with blank values
subgrp_rows <- data.frame(
  subgroup = "",
  var = subgrp,
  var_label = subgrp,
  col1_evt = NA,
  col1_median = NA,
  col2_evt = NA,
  col2_median = NA,
  hr_ci = NA,
  hr = NA,
  lcl = NA,
  ucl = NA,
  stringsAsFactors = FALSE
)

tbl_df <- rbind(df, subgrp_rows)
tbl_df$ord1 <- match(tbl_df$var, subgrp) # create order for each subgroup variable

## sort table in the defined order
tbl_df <- tbl_df |>
  mutate(
    ord1 = ifelse(is.na(ord1), 0, ord1),
    ord2 = ifelse(subgroup == "", 1, 2)
  ) |>
  arrange(ord1, ord2) |>
  mutate(
    subgroup = ifelse(subgroup == "All Patients", "ALL", subgroup),
    var_label = case_when(
      var == "ALL" ~ "All subjects",
      var == "AGEGR2" ~ "Age",
      var == "SEX" ~ "Sex",
      var == "RACEGR2" ~ "Race",
      var == "RENGRP1" ~ "Baseline renal function",
      var == "ECOGGR1" ~ "ECOG performance score",
      var == "STRAT01" ~ "Investigator's choice of DPd or DVd",
      var == "PRLNGRP" ~ "Number of lines of prior therapy",
      var == "ISSGRP" ~ "Baseline ISS",
      var == "PLASGRP" ~ "Baseline soft-tissue plasmacytomas",
      var == "BSLRISK" ~ "Cytogenetic risk groups",
      var == "PLSCGRP" ~ "Bone marrow % plasma cells"
    )
  ) |>
  # create ord variable for row ordering
  tibble::rownames_to_column(var = "ord") |>
  mutate(ord = as.numeric(ord))

################################################################################
# Generate plot:
################################################################################

# create sub-group label vector for y-axis
ytxt <- tbl_df |>
  mutate(
    row_lbl = ifelse(
      (subgroup == "" | subgroup == "ALL"),
      var_label,
      paste0("  ", subgroup)
    )
  ) |>
  pull(row_lbl)


## generate forest plot ---------------

plot <- ggplot(tbl_df, aes(x = hr, y = ord)) +
  # plot 95%CI
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), linewidth = 0.5, height = 0.5) +
  # plot hazard ratio points
  geom_point(shape = 16, size = 2) +
  # plot the vertical line of hr = 1.0
  geom_vline(xintercept = 1, color = "black") +

  # define title and axis labels
  scale_x_continuous(
    limits = c(0.01, 100),
    breaks = c(0.01, 0.1, 1, 10, 100),
    trans = "log10",
    labels = scales::label_number()
  ) +
  scale_y_continuous(
    breaks = c(1:nrow(tbl_df)),
    labels = ytxt,
    trans = "reverse"
  ) +
  theme_bw() +
  labs(
    title = "Hazard Ratio and 95% CI",
    x = "\u2190 Favor Dummy B             Favor Dummy A \u2192",
    y = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      size = 8,
      margin = margin(b = 0, unit = "cm")
    ),
    axis.title.x = element_text(hjust = 0.3, size = 8),
    axis.text.y = element_text(hjust = 0),
    panel.grid.minor = element_blank()
  )


## generate table plot --------------

tblpt <- ggplot(tbl_df, aes(y = ord)) +
  geom_text(aes(x = 1, label = col1_evt), size = 2.5) +
  geom_text(aes(x = 2, label = col1_median), size = 2.5) +
  geom_text(aes(x = 3, label = col2_evt), size = 2.5) +
  geom_text(aes(x = 4, label = col2_median), size = 2.5) +
  geom_text(aes(x = 5, label = hr_ci), size = 2.5) +
  scale_x_continuous(
    position = "top",
    labels = c(
      "Dummy A \n EVT/N",
      "Dummy A \n Median",
      "Dummy B \n EVT/N",
      "Dummy B \n Median",
      "Hazard Ratio \n (95% CI)"
    )
  ) +
  scale_y_continuous(trans = "reverse") +
  labs(y = "", x = "") +
  theme_minimal() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 7),
    # To expand the right side margin
    plot.margin = unit(c(0, 1, 0, 0), "cm"),
    panel.grid = element_blank()
  ) +
  # To display the complete results
  coord_cartesian(clip = "off")


# compose final object by putting plot, table, legend together --------

final <- (plot | tblpt) +
  plot_layout(widths = c(5, 5))


################################################################################
# Create png and output file:
################################################################################

# create png file and output figure
pname <- paste0(tolower(tblid), ".png")

# png(write_path(opath, pname)
png(
  write_path(opath, pname),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
# print(write_path(opath, pname)) ### print png path and name in log
print(write_path(opath, pname))
print(final)
dev.off()

# if (length( title_footer$main_footer) == 0) {
#   title_footer$main_footer <- NULL
# }

title <- "Forest Plot of Subgroup Analyses on Progression-free Survival based on Independent Review Committee (IRC) Assessment; Intent-to-Treat Analysis Set (Study 64007957MMY3001)"
main_footer <- ""

tidytlg::gentlg(
  tlf = "g",
  plotnames = write_path(opath, pname),
  plotwidth = 8,
  orientation = "landscape",
  opath = write_path(opath),
  file = tblid,
  title = title, # title_footer$title
  footers = main_footer
) # title_footer$main_footer
