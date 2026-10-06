###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsfvit02a.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create gsfvit02a: Boxplot of Vital Signs Over Time
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     advs
## Output:                    gsfvit02a.rtf
## Remarks:                   Include parameters: Systolic BP, Diastolic BP
##                            per FDA guidance box plots should not be used if 3 or more groups, so only present 2 treatment groups
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

tblid <- "GSFVIT02a"

# Safety population flag (default=SAFFL).
popfl <- "SAFFL"

# Actual treatment variable (default=TRT01A).
trtvar <- "TRT01A"

# List of Prameters to be included
paramcd_list <- c("SYSBP", "DIABP", "PULSE", "RESP", "TEMP")

# Parameter Labels list as per paramcd_list order
paramcd_lbls <- c(
  "Systolic Blood Pressure (mmHg)",
  "Diastolic Blood Pressure (mmHg)",
  "Pulse Rate (beats/min)",
  "Respiratory Rate (breaths/min)",
  "Temperature (C)"
)

# Plot title list for each parameter as per paramcd_list order
paramcd_plt_title <- paramcd_lbls # User has to Assign Labels if it is different from paramcd_lbls

# Split Timepoint into groups if not fit into single page
## User defined timepoint
timepoint_split <- list(
  time_point_1 = c("Baseline", "Cycle 02", "Cycle 03", "Cycle 04", "Cycle 05", "Cycle 06", "Cycle 07"),
  time_point_2 = c("Cycle 08", "Cycle 09", "Cycle 10", "Cycle 11", "Cycle 12", "Cycle 13", "Cycle 15"),
  time_point_3 = c("Cycle 17", "Cycle 19", "Cycle 21", "Cycle 22", "Cycle 23", "Cycle 25", "Cycle 30", "EOT")
)
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
advs_meta <- haven::read_sas(envsetup::read_path(a_in, "advs.sas7bdat")) |>
  df_na() |>
  distinct(PARAMCD, PARAMN) |>
  arrange(PARAMN)

# Preserve user-defined paramcd_list order for reordering labels after PARAMN sort
orig_order <- paramcd_list

# Reorder paramcd_list based on PARAMN from data
paramcd_list <- as.character(advs_meta$PARAMCD[advs_meta$PARAMCD %in% paramcd_list])

# Index to reorder paramcd_lbls and paramcd_plt_title to match PARAMN-sorted paramcd_list
idx <- match(paramcd_list, orig_order)
# Reorder parameter labels to match PARAMN-sorted paramcd_list
paramcd_lbls <- paramcd_lbls[idx]
# Reorder plot titles to match PARAMN-sorted paramcd_list
paramcd_plt_title <- paramcd_plt_title[idx]

# define function for data processing:
dt <- function(param_cd) {
  df <- haven::read_sas(envsetup::read_path(a_in, "advs.sas7bdat")) |>
    df_na() |>
    filter(PARAMCD == param_cd & toupper(.data[[popfl]]) == "Y" & ANL02FL == "Y") |>
    filter(!grepl("Unscheduled", AVISIT)) |>
    # per FDA guidance box plots should not be used if 3 or more groups
    # so only present 2 treatment groups here
    filter(.data[[trtvar]] %in% c("Xanomeline High Dose", "Placebo")) |>
    mutate(
      !!rlang::sym(trtvar) := factor(
        .data[[trtvar]],
        levels = c("Xanomeline High Dose", "Placebo")
      ),
      AVISIT = factor(
        .data[['AVISIT']],
        levels = unique(.data[['AVISIT']])[order(unique(.data[['AVISITN']]))]
      ),
      # shorten "End of Treatment" as EOT
      AVISIT = forcats::fct_recode(
        AVISIT,
        "EOT" = "End Of Treatment"
      ),
      # add abbreviation to TRT01A
      !!rlang::sym(trtvar) := forcats::fct_recode(
        .data[[trtvar]],
        "Xanomeline (High)" = "Xanomeline High Dose",
        "Placebo (PBO)" = "Placebo"
      ),
      !!rlang::sym(trtvar) := forcats::fct_relevel(
        .data[[trtvar]],
        "Xanomeline (High)",
        "Placebo (PBO)"
      )
    )

  # drop level for trtvar not to be presented in graphic
  df[[trtvar]] <- droplevels(df[[trtvar]])

  df2 <- df |>
    group_by(.data[[trtvar]], AVISITN, AVISIT) |>
    summarize(
      n = n(),
      mean = mean(AVAL, na.rm = TRUE),
      sd = sd(AVAL, na.rm = TRUE),
      se = sd(AVAL, na.rm = TRUE) / sqrt(n())
    ) |>
    ungroup() |>
    # fill in missing values for time points which don't have data
    tidyr::complete(
      !!rlang::sym(trtvar),
      tidyr::nesting(AVISIT),
      fill = list(n = 0, mean = 0.0, sd = 0.0, se = 0.0)
    ) |>
    # filling zero in decimal place for table display
    # mutate(meanC = sprintf('%.1f', mean)) |>
    mutate(meanC = tidytlg::roundSAS(mean, digits = 1, as_char = TRUE)) |>
    # reverse treatment group order: present from top to bottom in table
    mutate(trt_rev = forcats::fct_relevel(.data[[trtvar]], rev(levels(.data[[trtvar]]))))

  # Include time points only for subjects at least 10% left
  avisit_n <- df2 |>
    group_by(AVISIT) |>
    summarize(n = sum(n)) |>
    ungroup()

  base_n <- avisit_n |>
    filter(AVISIT == "Baseline") |>
    select(total = n)

  avisit_n <- cbind(avisit_n, base_n) |>
    mutate(prop = (n / total)) |>
    # remove time points of subjects less than 10%
    filter(prop >= 0.1) |>
    select(AVISIT)

  time_point <- as.vector(unlist(avisit_n["AVISIT"]))

  df <- df |>
    filter(AVISIT %in% time_point)

  df2 <- df2 |>
    filter(AVISIT %in% time_point)

  list(df, df2)
}

################################################################################
# Generate plot:
################################################################################

# define function for creating the body of graphic:

g_box <- function(df, df2) {
  # boxplot -------------------------------------------------------
  plot <- ggplot(
    df,
    aes(x = .data$AVISIT, y = .data$AVAL, fill = .data[[trtvar]])
  ) +
    geom_boxplot(color = "#E69F00") +
    scale_y_continuous(
      limits = c(df$ymin_l[1], df$ymax_l[1]),
      breaks = scales::extended_breaks(n = 5),
      labels = scales::label_number(accuracy = 1)
    ) +

    # assign colorblind friendly palette
    scale_fill_manual(values = cbbPalette) +
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

  # mean value table ----------------------------------------
  table_mean <- df2 |>
    ggplot(aes(x = .data$AVISIT, y = .data$trt_rev, label = .data$meanC)) +
    geom_text(size = 2.5) +

    # abbreviate table text label
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
      legend.position = "none",
    ) +
    labs(title = "Mean Value")

  # number of patients table ----------------------------------------
  table_n <- df2 |>
    ggplot(aes(x = .data$AVISIT, y = .data$trt_rev, label = .data$n)) +
    geom_text(size = 2.5) +

    # abbreviate table text label
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
      legend.position = "none",
    ) +
    labs(title = "Number of Subjects With Data")

  # compose final object by putting plot, table, legend together --------
  final <- plot /
    table_mean /
    table_n +
    plot_layout(heights = c(7.5, 1.3, 1.3)) +
    plot_annotation(
      title = param_title,
      theme = theme(plot.title = element_text(size = 10, family = "Arial")),
      subtitle = ""
    )
}

# Png files lists to add on tidytlg::gentlg()
png_list <- c()

# Creating PNG files as per paramcd_list
for (i in seq_along(paramcd_list)) {
  # generate data set for Parameters
  advs_param <- dt(paramcd_list[[i]])
  advs_param_pt <- advs_param[[1]] |>
    # global limits for Parameters
    mutate(ymin_l = floor(min(.data$AVAL, na.rm = TRUE) - 2), ymax_l = ceiling(max(.data$AVAL, na.rm = TRUE)) + 2)
  advs_param_tb <- advs_param[[2]]

  # define parameters for plotting:
  # present treatment group in y axis table from top to bottom (reverse)
  table_text <- c("PBO", "Xan High")

  # assign colorblind friendly palette to treatment groups: black(Xan Low), dark blue(PBO)
  cbbPalette <- c("#000000", "#0072B2")

  # Define axis labels and title for the plot
  x_label <- ""
  y_label <- paste0(paramcd_lbls[[i]])
  param_title <- paste0(paramcd_plt_title[[i]])

  split_active <- exists("timepoint_split") && !is.null(timepoint_split) && length(timepoint_split) > 0
  split_indices <- if (split_active) seq_along(timepoint_split) else 1

  if (split_active) {
    for (j in split_indices) {
      assign(paste0("advs_param_pt_", j), advs_param_pt |> filter(AVISIT %in% timepoint_split[[j]]))
      assign(paste0("advs_param_tb_", j), advs_param_tb |> filter(AVISIT %in% timepoint_split[[j]]))
    }
  } else {
    advs_param_pt_1 <- advs_param_pt
    advs_param_tb_1 <- advs_param_tb
  }

  for (k in split_indices) {
    assign(
      paste0("pt_param_", k),
      g_box(df = get(paste0("advs_param_pt_", k)), df2 = get(paste0("advs_param_tb_", k)))
    )
  }

  ################################################################################
  # Create png and output file:
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
    print(write_path(opath, pname_param)) ### print png path and name in log
    print(get(paste0("pt_param_", pn)))
    dev.off()

    png_list <- c(png_list, pname_param)
  }
}

tidytlg::gentlg(
  tlf = "g",
  plotnames = write_path(
    opath,
    png_list
  ),
  plotwidth = 8,
  orientation = "landscape",
  opath = write_path(opath),
  file = tblid,
  title = title_footer$title,
  footers = title_footer$main_footer
)
