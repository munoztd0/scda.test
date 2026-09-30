###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsids01a.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create gsids01a: Kaplan-Meier Plot of Time to
##                            Permanent Discontinuation of Study Treatment
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adttesaf
## Output:                    gsids01a.rtf
## Remarks:
## R-functions:
## R-function Sample Call:
##
## Modification History:
##  Rev #:
##  Modified By:
##  Reporting Effort:
##  Date:
##  Description:
################################################################################

################################################################################
# Prep environment:
################################################################################

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

################################################################################
# Define output ID:
################################################################################

tblid <- "GSIDS01a"

################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

# reading data

adttesaf <- haven::read_sas(envsetup::read_path(a_in, "adttesaf.sas7bdat")) |>
  df_na() |>
  filter(PARAMCD == "TTDCEOTD", SAFFL == "Y") |>
  # add abbreviation to TRT01A
  mutate(
    TRT01A = factor(
      case_when(
        TRT01A == "Xanomeline Low Dose" ~ "Xanomeline (Low)",
        TRT01A == "Xanomeline High Dose" ~ "Xanomeline (High)",
        TRT01A == "Placebo" ~ "Placebo (PBO)"
      ),
      levels = c("Xanomeline (Low)", "Xanomeline (High)", "Placebo (PBO)")
    )
  )

# fit KM model
fit <- survival::survfit(
  survival::Surv(AVAL, 1 - CNSR) ~ TRT01A,
  data = adttesaf
)

## convert model object to tidy data frame for plotting step lines
kmdata <- broom::tidy(fit) |>
  mutate(
    treatment = stringr::str_remove(strata, "TRT01A="),
    treatment = factor(
      treatment,
      levels = c("Xanomeline (Low)", "Xanomeline (High)", "Placebo (PBO)")
    ),
    estimate = estimate * 100
  ) |>
  select(-strata) |>
  arrange(treatment)

## filter censoring times for plotting censored points
cnsr_pts <- kmdata |>
  filter(n.censor > 0) |>
  mutate(cnsr = "censor")

################################################################################
# Generate plot:
################################################################################

# define parameters for plotting:
# assign colorblind friendly palette: black(Xan Low), orange(Xan High), dark blue(PBO)
cbbPalette <- c("#000000", "#E69F00", "#0072B2")

x_breaks <- seq(0, 1200, by = 100)
x_limits <- c(0, 1200)
y_breaks <- seq(0, 100, by = 10)
y_limits <- c(0, 100)
x_label <- "Days From First Dose"
y_label <- "Percentage of Subjects Continuing Treatment"

# present treatment group in y axis table from top to bottom (reverse)
table_text <- c("PBO", "Xan High", "Xan Low")


# survival plot ----------------------------------------------
## create step lines with censor points
plot <- kmdata |>
  ggplot(aes(x = time, y = estimate)) +
  # plot step_lines
  geom_step(aes(linetype = treatment, color = treatment)) +
  # plot censor points
  geom_point(aes(shape = cnsr), data = cnsr_pts) +

  # define censor point shapes
  # Use cross as choices for shape
  scale_shape_manual(values = 3) +

  # assign colorblind friendly palette
  scale_color_manual(values = cbbPalette) +

  # define x-axis breaks and limits
  scale_x_continuous(breaks = x_breaks, limits = x_limits) +
  # define y-axis breaks and limits
  scale_y_continuous(breaks = y_breaks, limits = y_limits) +
  labs(
    x = x_label,
    y = y_label
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 9, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    panel.grid = element_blank()
  )

## get breaks from the plot so they match the table
built_surv_plot <- ggplot_build(plot)
breaks <- built_surv_plot$layout$panel_params[[1]]$x$breaks


# Subjects at risk table -------------------------------------
## create risk table for table plot
risk_table <- summary(fit, times = breaks, extend = TRUE) |>
  with(data.frame(time, strata, n.risk)) |>
  mutate(
    treatment = stringr::str_remove(strata, "TRT01A="),
    treatment = factor(
      treatment,
      levels = c("Xanomeline (Low)", "Xanomeline (High)", "Placebo (PBO)")
    )
  ) |>
  select(-strata) |>
  arrange(treatment)

## create table plot
table_plot <- ggplot(aes(y = treatment), data = risk_table) +
  geom_text(aes(x = time, label = n.risk), size = 3) +

  # reverse y-scale, so the 1st level of trt is on top
  # abbreviate table text label
  scale_y_discrete(limits = rev, labels = table_text) +
  coord_cartesian(xlim = c(0, max(breaks))) +
  labs(title = "Subjects at Risk") +
  theme_minimal() +
  theme(
    title = element_text(size = 8),
    axis.text.y = element_text(size = 8, color = "black", hjust = 0.9),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
  )


# compose final object by putting plot, table together ------------
final <- plot / table_plot + plot_layout(heights = c(7.5, 1.5))


################################################################################
# Create png and output file:
################################################################################

# create png file and output figure
pname <- paste0(tolower(tblid), ".png")

png(
  write_path(opath, pname),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname)) ### print png path and name in log
print(final)
dev.off()

tidytlg::gentlg(
  tlf = "g",
  plotnames = write_path(opath, pname),
  plotwidth = 8,
  orientation = "landscape",
  opath = write_path(opath),
  file = tblid,
  title = title_footer$title,
  footers = title_footer$main_footer
)
