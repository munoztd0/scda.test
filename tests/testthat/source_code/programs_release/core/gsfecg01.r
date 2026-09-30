###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsfecg01.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create gsfecg01: Mean/Mean Change for ECG Parameters Over Time
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adeg
## Output:                    gsfecg01.rtf
## Remarks:                   Include parameters: ECG Mean Heart Rate, PR Interval Aggregate,
##                            RR Interval Aggregate, QRS Duration Aggregate, QT Interval Aggregate,
##                            QTcF Interval Aggregate, QTcB Interval Aggregate
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

tblid <- "GSFECG01"

# Safety population flag (default = SAFFL).
popfl <- "SAFFL"

# Percentage of subjects remain in study as per requirement
proportion <- 0.1

# Actual treatment variable (default = TRT01A).
trtvar <- "TRT01A"

# mean(AVAL)/mean change(CHG)
analvar <- "CHG"

# Actual treatment levels (required to order the columns).
trtlev <- c("Xanomeline Low Dose", "Xanomeline High Dose", "Placebo")

# List of parameter codes to be included.
paramcd_list <- c("EGHRMN", "PRAG", "RRAG", "QRSAG", "QTC", "QTCFAG", "QTCBAG")

# Parameter labels in the same order as paramcd_list.
paramcd_lbls <- c(
  "ECG Mean Heart Rate (beats/min)",
  "PR Interval, Aggregate (ms)",
  "RR Interval, Aggregate (ms)",
  "QRS Duration, Aggregate (ms)",
  "QT Interval, Aggregate (ms)",
  "QTcF Interval, Aggregate (ms)",
  "QTcB Interval, Aggregate (ms)"
)

# Plot titles for each parameter, in the same order as paramcd_list.
# Assign different values here if plot titles should differ from paramcd_lbls.
paramcd_plt_title <- paramcd_lbls

# Split Timepoint into groups if not fit into single page
## User defined timepoint
# timepoint_split <- list(
#   time_point_1 = c("Baseline", "Month 1", "Month 3", "Month 6", "Month 9"),
#   time_point_2 = c("Month 12", "Month 15", "Month 18", "Month 24")
# )

################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################
# Validate paramcd_list and paramcd_lbls are aligned in length and order
if (length(paramcd_list) != length(paramcd_lbls)) {
  stop("paramcd_list and paramcd_lbls must have the same number of elements in the same order.")
}

# Read PARAMCD and PARAMN from input dataset to determine the data-driven parameter order
adeg_meta <- haven::read_sas(envsetup::read_path(a_in, "adeg.sas7bdat")) |>
  df_na() |>
  distinct(PARAMCD, PARAMN) |>
  arrange(PARAMN)

# Preserve user-defined paramcd_list order for reordering labels after PARAMN sort
orig_order <- paramcd_list

# Reorder paramcd_list based on PARAMN from data
paramcd_list <- as.character(adeg_meta$PARAMCD[adeg_meta$PARAMCD %in% paramcd_list])

# Index to reorder paramcd_lbls and paramcd_plt_title to match PARAMN-sorted paramcd_list
idx <- match(paramcd_list, orig_order)
# Reorder parameter labels to match PARAMN-sorted paramcd_list
paramcd_lbls <- paramcd_lbls[idx]
# Reorder plot titles to match PARAMN-sorted paramcd_list
paramcd_plt_title <- paramcd_plt_title[idx]

# Function to read, filter, and summarize ADEG data for a given parameter.
dt <- function(param) {
  adeg <- haven::read_sas(envsetup::read_path(a_in, "adeg.sas7bdat")) |>
    df_na() |>
    filter(PARAMCD == param & toupper(.data[[popfl]]) == "Y") |>
    filter(ABLFL == 'Y' | ANL02FL == 'Y') |>
    mutate(
      !!trtvar := factor(
        .data[[trtvar]],
        levels = trtlev
      ),
      AVISIT = factor(
        .data[["AVISIT"]],
        levels = unique(.data[["AVISIT"]])[order(unique(.data[["AVISITN"]]))]
      )
    ) |>
    mutate(
      !!trtvar := forcats::fct_recode(
        !!sym(trtvar),
        "Xanomeline (Low)" = "Xanomeline Low Dose",
        "Xanomeline (High)" = "Xanomeline High Dose",
        "Placebo (PBO)" = "Placebo"
      )
    ) |>
    group_by(.data[[trtvar]], AVISITN, AVISIT) |>
    summarize(
      n = n(),
      # If analvar is AVAL, compute mean/sd/se from AVAL; if CHG, compute from CHG
      mean = if (analvar == "AVAL") mean(AVAL, na.rm = TRUE) else mean(CHG, na.rm = TRUE),
      sd = if (analvar == "AVAL") sd(AVAL, na.rm = TRUE) else sd(CHG, na.rm = TRUE),
      se = if (analvar == "AVAL") sd(AVAL, na.rm = TRUE) / sqrt(n()) else sd(CHG, na.rm = TRUE) / sqrt(n()),
      # If analvar is CHG, also compute mean/sd/se from AVAL for baseline values
      mean_avl = mean(AVAL, na.rm = TRUE),
      sd_avl = sd(AVAL, na.rm = TRUE),
      se_avl = sd(AVAL, na.rm = TRUE) / sqrt(n())
    ) |>
    ungroup() |>
    # Fill in missing values for visits with no data.
    tidyr::complete(
      .data[[trtvar]],
      tidyr::nesting(AVISIT),
      fill = list(n = 0, mean = 0.0, sd = 0.0, se = 0.0, mean_avl = 0.0, sd_avl = 0.0, se_avl = 0.0)
    ) |>
    mutate(
      meanC = tidytlg::roundSAS(mean, digits = 1, as_char = TRUE),
      meanA = tidytlg::roundSAS(mean_avl, digits = 1, as_char = TRUE),
      # Combined mean change / mean value label for display.
      meanP = paste0(meanC, "/", meanA)
    ) |>
    # Reverse treatment group order for top-to-bottom display in the table.
    mutate(trt_rev = forcats::fct_relevel(!!sym(trtvar), rev(levels(!!sym(trtvar)))))

  # Retain only visits where at least 10% of baseline subjects have data.
  avisit_n <- adeg |>
    group_by(AVISIT) |>
    summarize(n = sum(n)) |>
    ungroup()

  baseline_n <- avisit_n |>
    filter(AVISIT == "Baseline") |>
    select(total = n)

  avisit_n <- cbind(avisit_n, baseline_n) |>
    mutate(prop = n / total) |>
    filter(prop >= proportion) |>
    select(AVISIT)

  time_point <- as.vector(unlist(avisit_n["AVISIT"]))

  adeg |>
    filter(AVISIT %in% time_point)
}

################################################################################
# Generate plot:
################################################################################

# Function to create the line plot with mean value and subject count tables.
g_line <- function(df) {
  # Line plot with 95% CI error bars -------------------------------------------
  plot <- ggplot(
    df,
    aes(
      x = .data$AVISIT,
      y = .data$mean,
      group = .data[[trtvar]],
      color = .data[[trtvar]],
      linetype = .data[[trtvar]],
      shape = .data[[trtvar]]
    )
  ) +
    geom_errorbar(
      aes(
        ymin = .data$mean - 1.96 * .data$se,
        ymax = .data$mean + 1.96 * .data$se
      ),
      width = 0.1,
      position = pd
    ) +
    geom_line(aes(x = AVISIT, y = mean), position = pd) +
    geom_point(position = pd) +
    scale_y_continuous(
      limits = c(df$ymin_l[1], df$ymax_l[1]),
      breaks = scales::extended_breaks(n = 5),
      labels = scales::label_number(accuracy = 1)
    ) +
    # Colorblind-friendly palette.
    scale_color_manual(values = cbbPalette) +
    {
      if (analvar == "CHG") geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) else NULL
    } +
    labs(x = x_label, y = y_label) +
    theme_bw() +
    theme(
      text = element_text(size = 9, color = "black", family = "Arial"),
      axis.text = element_text(size = 9, color = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9, face = "bold")
    )

  # Mean value annotation table ------------------------------------------------
  table_mean <- df |>
    ggplot(aes(x = .data$AVISIT, y = .data$trt_rev, label = if (analvar == "AVAL") .data$meanC else .data$meanP)) +
    geom_text(size = 2.5) +
    scale_y_discrete(labels = table_text) +
    theme_bw() +
    theme(
      title = element_text(size = 8, face = "bold"),
      text = element_text(family = "Arial"),
      axis.text = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    labs(title = if (analvar == "AVAL") "Mean Value" else "Mean Change from Baseline / Mean Value")

  # Subject count annotation table ---------------------------------------------
  table_n <- df |>
    ggplot(aes(x = .data$AVISIT, y = .data$trt_rev, label = .data$n)) +
    geom_text(size = 2.5) +
    scale_y_discrete(labels = table_text) +
    theme_bw() +
    theme(
      title = element_text(size = 8, face = "bold"),
      text = element_text(family = "Arial"),
      axis.text = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    labs(
      title = if (analvar == "AVAL") {
        "Number of Patients with Data"
      } else {
        "Number of Subjects With Both Baseline and Postbaseline Data"
      }
    )

  # Compose final figure: line plot + mean table + subject count table ---------
  final <- plot /
    table_mean /
    table_n +
    plot_layout(heights = c(7.5, 1.3, 1.3)) +
    plot_annotation(
      title = param_title,
      subtitle = ""
    )

  # Apply title theme to suppress patchwork warnings.
  final + theme(plot.title = element_text(size = 10, family = "Arial"))
}

# Initialise list of PNG file names to pass to gentlg().
png_list <- c()

# Loop over each parameter, generate and save a PNG plot.
for (i in seq_along(paramcd_list)) {
  adeg_hr <- dt(paramcd_list[[i]]) |>
    # Compute global y-axis limits across all visits for this parameter.
    mutate(
      ymin_l = floor(min(.data$mean - 1.92 * .data$se, na.rm = TRUE) - 2),
      ymax_l = ceiling(max(.data$mean + 1.92 * .data$se, na.rm = TRUE) + 2)
    )

  # Abbreviated treatment labels for the annotation tables (top to bottom order).
  table_text <- c("PBO", "Xan High", "Xan Low")

  # Colorblind-friendly palette: black (Xan Low), orange (Xan High), dark blue (PBO).
  cbbPalette <- c("#000000", "#E69F00", "#0072B2")

  pd <- position_dodge(0.3)

  x_label <- ""
  y_label <- paste0(
    if (analvar == "AVAL") "Mean Value (95% CI)\n" else "Mean Change From Baseline (95% CI)\n",
    paramcd_lbls[[i]]
  )
  param_title <- paramcd_plt_title[[i]]

  # Determine split indices
  split_active <- exists("timepoint_split") && !is.null(timepoint_split) && length(timepoint_split) > 0
  split_indices <- if (split_active) seq_along(timepoint_split) else 1

  # Create data frames for each split (or single frame if no split)
  if (split_active) {
    for (j in split_indices) {
      assign(
        paste0("adeg_param_", j),
        adeg_hr |> filter(AVISIT %in% timepoint_split[[j]])
      )
    }
  } else {
    adeg_param_1 <- adeg_hr
  }

  # Generate plots for each timepoint split (or single plot if no split)
  for (k in split_indices) {
    assign(
      paste0("pt_param_", k),
      g_line(df = get(paste0("adeg_param_", k)))
    )
  }

  ################################################################################
  # Save PNG and add to output list:
  ################################################################################

  for (pn in split_indices) {
    pname_param <- paste0(tolower(tblid), "_", tolower(paramcd_list[[i]]), "_", pn, ".png")

    png(
      write_path(opath, pname_param),
      width = 22,
      height = 14,
      units = "cm",
      res = 300,
      type = "cairo"
    )
    print(write_path(opath, pname_param)) # Log the PNG file path and name.
    print(get(paste0("pt_param_", pn)))
    dev.off()

    png_list <- c(png_list, pname_param)
  }
}

################################################################################
# Assemble RTF output:
################################################################################

tidytlg::gentlg(
  tlf = "g",
  plotnames = write_path(opath, png_list),
  plotwidth = 8,
  orientation = "landscape",
  opath = write_path(opath),
  file = tblid,
  title = title_footer$title,
  footers = title_footer$main_footer
)
