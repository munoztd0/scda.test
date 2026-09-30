###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gefmad06.r
## R Version:                 4.5.2
## Short Description:         Cumulative Response of Percent Change From Baseline to Time Point
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:
## Output:
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


library(ggplot2)

################################################################################
# Define output ID:
################################################################################

tblid <- "GEFXXX06"

################################################################################
# Get titles and footnotes:
################################################################################

# title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

# reading data
admadrs <- haven::read_sas(read_path(a_in, "admadrs.sas7bdat")) |>
  filter(
    FAS1FL == "Y",
    PARAMCD == "MADR0212",
    ANL02FL == "Y",
    AVISIT == "Day 43",
    DTYPE == "",
    APHASE == "Double Blind Phase"
  ) |>
  mutate(
    trt = factor(
      case_when(
        TRT01P == "Placebo" ~ "Placebo",
        TRT01P == "Seltorexant 20 mg" ~ "Selt 20 mg"
      ),
      levels = c("Placebo", "Selt 20 mg")
    ),
    val = -PCHG
  )

# calculate cumulative response of percent change from baseline
admadrs2 <- admadrs |>
  group_by(trt, val) |>
  mutate(count = row_number()) |>
  slice_tail(n = 1) |>
  arrange(trt, desc(val)) |>
  ungroup()

admadrs3 <- admadrs2 |>
  group_by(trt) |>
  mutate(
    cum_freq = cumsum(count),
    cum_pct = (cum_freq / sum(count)) * 100
  ) |>
  ungroup() |>
  select(USUBJID, trt, val, count, cum_freq, cum_pct)

# add (N=) in treatment group label
admadrs_n <- admadrs |>
  group_by(trt) |>
  summarize(total = n()) |>
  ungroup()

admadrs4 <- admadrs3 |>
  left_join(admadrs_n, by = "trt") |>
  mutate(trt = paste0(trt, " (N=", total, ")"))

################################################################################
# Generate plot:
################################################################################

# define parameters for plotting:
# assign colorblind friendly palette: black(Placebo), orange(Selt 20 mg)
cbbPalette <- c("#000000", "#E69F00")

x_breaks <- seq(-40, 100, by = 20)
x_limits <- c(-40, 100)
y_breaks <- seq(0, 100, by = 10)
y_limits <- c(0, 100)
x_label <- "Percent Reduction in MADRS Total Score at Day 43"
y_label <- "Cumulative Percentage of Subjects"


# line plot ----------------------------------------------

pt <- admadrs4 |>
  ggplot(aes(x = val, y = cum_pct, group = trt, color = trt, linetype = trt)) +
  geom_line() +

  # assign colorblind friendly palette
  scale_color_manual(values = cbbPalette) +

  # define x-axis breaks and limits
  scale_x_continuous(breaks = x_breaks, limits = x_limits) +
  # define y-axis breaks and limits
  scale_y_continuous(breaks = y_breaks, limits = y_limits) +
  geom_vline(xintercept = 50, color = "grey", linetype = 2) +
  labs(
    x = x_label,
    y = y_label
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 9, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    axis.ticks = element_blank(),
    legend.position = c(0.95, 0.95), # Position inside the plot
    legend.justification = c("right", "top"), # Align the top right of the legend box
    legend.box.background = element_rect(color = "black", fill = "lightgrey"), # Square box around legend
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    panel.grid = element_blank()
  )


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
print(pt)
dev.off()

title <- "MADRS Total Score: Cumulative Response of Percent Change From Baseline (DB) to Day 43: Observed - Double-blind Phase; FAS1 Analysis Set (Study 42847922MDD3002)"
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
