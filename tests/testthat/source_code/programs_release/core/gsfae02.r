###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsfae02.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create tsfae20a: Cumulative Incidence Plot of Time to
##                            [Adverse Event of Special Interest]
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adttesaf
## Output:                    gsfae02.rtf
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

tblid <- "GSFAE02"
time_label <- "Days"

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
  filter(PARAMCD == "TTAELPTD", SAFFL == "Y") |>
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


################################################################################
# Generate plot:
################################################################################

# define parameters for plotting:
# assign colorblind friendly palette: black(Xan Low), orange(Xan High), dark blue(PBO)
cbbPalette <- c("#000000", "#E69F00", "#0072B2")

x_label <- paste(time_label, "From First Dose")
y_label <- "Cumulative incidence (%) \n (95% CI)"

# present treatment group in y axis table from top to bottom (reverse)

table_text <- c("Placebo", "Xanomeline (High)", "Xanomeline (Low)")


# survival plot ----------------------------------------------
## convert model object to tidy data frame for plotting step lines
kmdata <- broom::tidy(fit) |>
  mutate(
    treatment = stringr::str_remove(strata, "TRT01A="),
    treatment = factor(
      treatment,
      levels = c("Xanomeline (Low)", "Xanomeline (High)", "Placebo (PBO)")
    )
  ) |>
  select(-strata) |>
  arrange(treatment)

## calculate the percentage of cumulative counts
cum <- kmdata |>
  group_by(treatment) |>
  mutate(cum_count = cumsum(n.event)) |>
  ungroup()

total <- adttesaf |>
  group_by(TRT01A) |>
  summarise(total = n_distinct(USUBJID)) |>
  rename(treatment = TRT01A) |>
  ungroup()

cumdata <- cum |>
  left_join(total, by = "treatment") |>
  mutate(cum = (cum_count / total) * 100)

# x-axis breaks used in the plot
x_breaks <- scales::extended_breaks(n = 10)(cumdata$time)

# keep nearest KM estimate to each x-axis tick per Treatment group
cumdata_err <- bind_rows(
  lapply(x_breaks, function(x) {
    cumdata |>
      group_by(treatment) |>
      slice(which.min(abs(time - x)))
  })
) |>
  distinct(treatment, time, .keep_all = TRUE)

## create step lines plot
plot <- cumdata |>
  ggplot(aes(x = time, y = cum, color = treatment)) +
  # plot step_lines
  geom_step(aes(linetype = treatment)) +

  # CI intervals only at x-axis ticks

  geom_errorbar(
    data = cumdata_err,
    aes(
      ymin = (1 - conf.low) * 100,
      ymax = (1 - conf.high) * 100,
      linetype = treatment
    ),
    width = 0.1
  ) +

  # assign colorblind friendly palette
  scale_color_manual(values = cbbPalette) +

  # Use scales package for automatic breaks calculation
  scale_x_continuous(
    breaks = scales::extended_breaks(n = 10),
    labels = scales::label_number(accuracy = 1, big.mark = "")
  ) +
  # Use scales package for automatic breaks calculation
  scale_y_continuous(
    breaks = scales::extended_breaks(n = 5),
    labels = scales::label_number(accuracy = 1)
  ) +
  labs(
    x = x_label,
    y = y_label
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 9, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    axis.ticks = element_line(),
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
  labs(title = "Number of Subjects at Risk") +
  theme_bw() +
  theme(
    title = element_text(size = 8),
    axis.text.y = element_text(size = 8, color = "black", hjust = 0.9),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
  )

# Cumulative number of Subjects with Event -------------------------------------
## create cumulative number of subjects with event for table plot
cum_table <- summary(fit, times = breaks, extend = TRUE) |>
  with(data.frame(time, strata, n.event)) |>
  mutate(
    treatment = stringr::str_remove(strata, "TRT01A="),
    treatment = factor(
      treatment,
      levels = c("Xanomeline (Low)", "Xanomeline (High)", "Placebo (PBO)")
    )
  ) |>
  select(-strata) |>
  arrange(treatment)

## calculate cumulative counts per each treatment
cum_table <- cum_table |>
  group_by(treatment) |>
  mutate(cum_sum = cumsum(n.event)) |>
  ungroup()

## create table plot
table_plot2 <- ggplot(aes(y = treatment), data = cum_table) +
  geom_text(aes(x = time, label = cum_sum), size = 3) +

  # reverse y-scale, so the 1st level of trt is on top
  # abbreviate table text label
  scale_y_discrete(limits = rev, labels = table_text) +
  coord_cartesian(xlim = c(0, max(breaks))) +
  labs(title = "Cumulative Number of Subjects With Event") +
  theme_bw() +
  theme(
    title = element_text(size = 8),
    axis.text.y = element_text(size = 8, color = "black", hjust = 0.9),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
  )

# compose final object by putting plot, table together ------------
final <- plot /
  table_plot /
  table_plot2 +
  plot_layout(heights = c(7.5, 1.3, 1.3))


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

if (length(title_footer$main_footer) == 0) {
  title_footer$main_footer <- NULL
}

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
