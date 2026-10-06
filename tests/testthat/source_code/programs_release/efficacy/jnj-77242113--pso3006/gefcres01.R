###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gefcres01.r
## R Version:                 4.5.2
## junco Version:             0.1.6
## Short Description:         Program to create gefcres01:
##                            Subjects in [Response Description] at [Time Point] in
##                            [Treatment Arm] Versus [Control Arm] (Tipping Point
##                            Analysis Based on Multiple Imputation with Bernoulli Draws)
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adpasiri
## Output:                    gefcres01.rtf
## Remarks:
##
## Modification History:
## Rev #:
## Modified By:
## Reporting effort:
## Date:
## Description:
################################################################################

# Environment ----

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)

library(dplyr)
library(rtables)
library(junco)
library(haven)
library(ggplot2)
library(tidytlg)

# Parameters ----

# Test run or real production run? Please set to FALSE for use in the study analysis.
testrun <- TRUE

# Define output ID and file location.
tblid <- "GEFCRES01"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Ustekinumab to JNJ-77242113"

# Define the treatment group label.
trtlab <- "JNJ-77242113"

# Plot axis label for control group.
commonlab <- "Assumed response rate among subjects with missing data\n (after accounting for ICEs) in"
ctrplotlab <- paste(commonlab, "Ustekinumab group")

# Plot axis label for treatment group.
trtplotlab <- paste(commonlab, trtlab, "group")

# Define population flags used.
popfl <- "FASFL"

# Define the stratification variables to use.
stratvar <- c("STRAT01", "STRAT02", "STRAT03")

# Define visit at which to analyze the results.
timepoint <- "Week 28"

# Define response parameter.
# Note that a label is not needed for this plot.
resppar <- "PASI90S"

# Missing data response probability grid.
respprob_grid <- (0:10) / 10

# Number of multiple imputations to use.
n_imputations <- if (testrun) 10 else 200

# Define p-value categories alongside their colors for the tile shading in the
# plot.
pvalcat <- list(
  "<0.001" = c(0, 0.001),
  "0.001 to <0.05" = c(0.001, 0.05),
  "\u22650.05" = c(0.05, 1) # \u2265 is the unicode for >=
)
pvalcols <- c(
  "<0.001" = "white",
  "0.001 to <0.05" = "gray90",
  "\u22650.05" = "gray60" # \u2265 is the unicode for >=
)

# Formats specifications.
formats <- list(
  est_prop = jjcsformat_xx("xx%"), # For the x- and y-axis tick labels.
  est_prop_diff = jjcsformat_xx("xx.x%") # For the tile contents.
)

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(!!rlang::sym(trtvar) %in% c(ctrlab, trtlab)) |>
  mutate(!!trtvar := factor(.data[[trtvar]], levels = c(ctrlab, trtlab))) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(stratvar))

## ADPASIRI ----

adpasiri <- haven::read_sas(read_path(a_in, "adpasiri.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(
    PARAMCD == resppar,
    AVISIT == timepoint
  ) |>
  select(STUDYID, USUBJID, PARAMCD, AVALC) |>
  mutate(
    response = case_when(
      AVALC == "Y" ~ TRUE,
      AVALC == "N" ~ FALSE,
      .default = NA
    )
  )

## Analysis ----

ana <- adpasiri |>
  inner_join(adsl, by = c("STUDYID", "USUBJID"))

# Multiple Imputation ----

# Define 2-dimensional probability grid.
respprob_grid_2d <- expand.grid(respprob_grid, respprob_grid)
colnames(respprob_grid_2d) <- c(ctrlab, trtlab)

# Obtain results.
set.seed(891)
scenario_results <- resp_multiple_imputation(
  dat = ana,
  p_ctrl = respprob_grid_2d[[ctrlab]],
  p_trt = respprob_grid_2d[[trtlab]],
  trtvar = trtvar,
  ctrlab = ctrlab,
  trtlab = trtlab,
  respvar = "response",
  stratvar = stratvar,
  n_imputations = n_imputations,
  pvalcat = pvalcat
)

# Ensure ordered factors for plotting
plot_df <- scenario_results |>
  mutate(
    p_ctrl_f = factor(p_ctrl, levels = respprob_grid),
    p_trt_f = factor(p_trt, levels = respprob_grid),
    p_cat = factor(p_cat, levels = names(pvalcat))
  )

# Plot ----

axis_labels <- respprob_grid |>
  scales::label_percent(accuracy = 1)()

p <- ggplot(plot_df, aes(x = p_trt_f, y = p_ctrl_f, fill = p_cat)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = effect_label), size = 8 / .pt, family = "Arial") +
  scale_fill_manual(values = pvalcols, drop = FALSE, name = "p-value") +
  scale_x_discrete(labels = axis_labels) +
  scale_y_discrete(labels = axis_labels) +
  labs(
    x = trtplotlab,
    y = ctrplotlab
  ) +
  theme_minimal(base_size = 9) +
  theme(
    text = element_text(family = "Arial"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.position = "right",
    legend.key = element_rect(color = "black", linewidth = 0.2)
  )

# Output ----

title_footer <- tab_titles

# Export to png.
png_name <- paste0(tolower(tblid), ".png")
png_file <- write_path(opath, png_name)
ggsave(
  png_file,
  p,
  height = 6,
  width = 8,
  units = "in",
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
  plotwidth = 8,
  plotheight = 6,
  orientation = "landscape",
  opath = write_path(opath),
  file = rtf_name,
  title = title_footer$title,
  footers = title_footer$main_footer
)
