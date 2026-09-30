###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              gefos01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Kaplan-Meier Plot for Overall Survival
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADTTEEF
## Output:                    GEFOS01.png
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

library(survival)

################################################################################
# Define output ID:
################################################################################

tblid <- "GEFOS01"

################################################################################
# Get titles and footnotes:
################################################################################

# title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

# reading data

adtte <- haven::read_sas(read_path(a_in, "adtteef.sas7bdat")) |>
  filter(ITTFL == "Y", PARAMCD == "OSM") |>
  mutate(
    trt = factor(
      case_when(
        TRT01P == "Dummy A - Tec-Dara" ~ "Dummy A",
        TRT01P == "Dummy B - DPd" ~ "Dummy B",
        TRT01P == "Dummy B - DVd" ~ "Dummy B"
      ),
      levels = c("Dummy B", "Dummy A")
    )
  )

# fit KM model
fit <- survfit(Surv(AVAL, 1 - CNSR) ~ trt, data = adtte)

## convert model object to tidy data frame for plotting step lines
kmdata <- broom::tidy(fit) |>
  mutate(
    treatment = stringr::str_remove(strata, "trt="),
    treatment = factor(treatment, levels = c("Dummy A", "Dummy B")),
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
# assign colorblind friendly palette: black(Dummy B), orange(Dummy A)
cbbPalette <- c("#000000", "#E69F00")

x_breaks <- seq(0, 36, by = 3)
x_limits <- c(0, 36)
y_breaks <- seq(0, 100, by = 20)
y_limits <- c(0, 100)
x_label <- "Months From Randomization Date (months)"
y_label <- "% of Subjects Alive"

# present treatment group in y axis table from top to bottom (reverse)
table_text <- c("Dummy B", "Dummy A")


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
    treatment = stringr::str_remove(strata, "trt="),
    treatment = factor(treatment, levels = c("Dummy A", "Dummy B"))
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

title <- "Kaplan-Meier Plot for Overall Survival; Intent-to-Treat Analysis Set (Study 64007957MMY3001)"
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
