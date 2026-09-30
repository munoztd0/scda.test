###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsflab01.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create gsflab01: Mean Change From Baseline for
##                            [Laboratory Category] Laboratory Data Over Time
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adlb or adlbc or adlc
## Output:                    gsflab01chm.rtf
##                            gsflab01hem.rtf
## Remarks:                   Include categories: CHEMISTRY, HEMATOLOGY
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

tblid <- "GSFLAB01"

# Safety population flag (default = SAFFL).
popfl <- "SAFFL"

# Percentage of subjects remain in study as per requirement
proportion <- 0.1

# Error bar type for the plot.
# Options: "CI" (95% confidence interval), "SE" (standard error), "SD" (standard deviation).
error_type <- "CI"

# Actual treatment variable (default = TRT01A).
trtvar <- "TRT01A"

# Define parameter categories with their RTF suffix
paramcat <- list(
  chm = "CHEMISTRY",
  hem = "HEMATOLOGY"
)

# List of parameter codes to be included.
# If NULL, all PARAMCD values from adlb data will be used (sorted by PARCAT1, PARCAT3N, PARAM).
paramcd_list <- NULL

# Parameter labels in the same order as paramcd_list.
# If NULL, PARAM values from adlb data will be used.
paramcd_lbls <- NULL

# Plot titles for each parameter, in the same order as paramcd_list.
# If NULL, paramcd_lbls will be used.
paramcd_plt_title <- NULL

# Split Timepoint into groups if not fit into single page
# default = NULL (no split)
timepoint_split <- NULL

## user can comment this part if no split needed
# else User can define timepoints here
# timepoint_split <- list(
#   time_point_1 = c("Baseline", "Cycle 02", "Cycle 03", "Cycle 04", "Cycle 05"),
#   time_point_2 = c("Cycle 07",  "Cycle 09",  "Cycle 11", "Cycle 13", "EOT")
# )

################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

adlb_full <- haven::read_sas(envsetup::read_path(a_in, "adlb.sas7bdat")) |>
  df_na() |>
  filter(PARCAT1 %in% unlist(paramcat), !is.na(PARCAT3), ANL02FL == "Y", toupper(.data[[popfl]]) == "Y")

# define function for data processing:
dt <- function(param_cd) {
  adlb <- adlb_full |>
    filter(PARAMCD == param_cd) |>
    mutate(
      !!trtvar := factor(
        .data[[trtvar]],
        levels = c("Xanomeline Low Dose", "Xanomeline High Dose", "Placebo")
      ),
      AVISIT = factor(
        .data[["AVISIT"]],
        levels = unique(.data[["AVISIT"]])[order(unique(.data[["AVISITN"]]))]
      )
    ) |>
    mutate(
      # shorten "End of Treatment" as EOT
      AVISIT = forcats::fct_recode(AVISIT, "EOT" = "End Of Treatment"),
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
      mean = mean(CHG, na.rm = TRUE),
      sd = sd(CHG, na.rm = TRUE),
      se = sd(CHG, na.rm = TRUE) / sqrt(n()),
      mean_avl = mean(AVAL, na.rm = TRUE),
      sd_avl = sd(AVAL, na.rm = TRUE),
      se_avl = sd(AVAL, na.rm = TRUE) / sqrt(n())
    ) |>
    ungroup() |>
    # fill in missing values for time points which don't have data
    tidyr::complete(
      .data[[trtvar]],
      tidyr::nesting(AVISIT),
      fill = list(n = 0, mean = 0.00, sd = 0.00, se = 0.00, mean_avl = 0.00, sd_avl = 0.00, se_avl = 0.00)
    ) |>
    mutate(
      # User can update decimal places by modifying the digits parameter in roundSAS function
      meanC = tidytlg::roundSAS(mean, digits = 2, as_char = TRUE),
      meanA = tidytlg::roundSAS(mean_avl, digits = 2, as_char = TRUE),
      meanP = paste0(meanC, "/", meanA)
    ) |>
    # reverse treatment group order: present from top to bottom in table
    mutate(trt_rev = forcats::fct_relevel(!!sym(trtvar), rev(levels(!!sym(trtvar))))) |>
    # assign EOT as NA (avoid to connect line to last time point)
    mutate(new_avisit = na_if(AVISIT, "EOT"))

  # Include time points only for subjects at least 10% left
  avisit_n <- adlb |>
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

  adlb |>
    filter(AVISIT %in% time_point)
}

################################################################################
# Generate plot:
################################################################################

# define function for creating the body of graphic:
g_line <- function(df) {
  # line plot -------------------------------------------------------
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
        ymin = .data$mean - err_mult * .data[[err_col]],
        ymax = .data$mean + err_mult * .data[[err_col]]
      ),
      width = 0.1,
      position = pd
    ) +
    geom_line(
      data = df[!is.na(df$new_avisit), ],
      aes(x = new_avisit, y = mean),
      position = pd
    ) +
    geom_point(position = pd) +
    scale_y_continuous(
      limits = if (all(c("ymin_l", "ymax_l") %in% names(df))) c(df$ymin_l[1], df$ymax_l[1]) else NULL,
      breaks = scales::extended_breaks(n = 5),
      # If user needs to display the decimals as per data range on y axis, comment out the labels parameter below
      # labels = scales::label_number(accuracy = 0.1)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
    # assign colorblind friendly palette
    scale_color_manual(values = cbbPalette) +
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
  table_mean <- df |>
    ggplot(aes(x = .data$AVISIT, y = .data$trt_rev, label = .data$meanP)) +
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
    labs(title = "Mean Change From Baseline / Mean Value")

  # number of patients table ----------------------------------------
  table_n <- df |>
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
    labs(title = "Number of Subjects With Both Baseline and Postbaseline Data")

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

################################################################################
# Generate PNGs and RTF per category:
################################################################################

for (cat_suffix in names(paramcat)) {
  cat_value <- paramcat[[cat_suffix]]

  # Derive paramcd_list/lbls/title scoped to this category
  adlb_meta <- adlb_full |>
    filter(PARCAT1 == cat_value) |>
    distinct(PARAMCD, PARAM, PARCAT1, PARCAT3N) |>
    arrange(PARCAT1, PARCAT3N, PARAM)

  cat_paramcd_list <- if (!is.null(paramcd_list)) {
    paramcd_list[paramcd_list %in% adlb_meta$PARAMCD]
  } else {
    adlb_meta$PARAMCD
  }
  cat_paramcd_lbls <- if (!is.null(paramcd_lbls)) {
    paramcd_lbls[match(cat_paramcd_list, paramcd_list)]
  } else {
    adlb_meta$PARAM[match(cat_paramcd_list, adlb_meta$PARAMCD)]
  }
  cat_paramcd_title <- if (!is.null(paramcd_plt_title)) {
    paramcd_plt_title[match(cat_paramcd_list, paramcd_list)]
  } else {
    cat_paramcd_lbls
  }

  png_list <- c()

  table_text <- c("PBO", "Xan High", "Xan Low")

  cbbPalette <- c("#000000", "#E69F00", "#0072B2")

  pd <- position_dodge(0.3)

  err_col <- switch(error_type, CI = "se", SE = "se", SD = "sd")
  err_mult <- switch(error_type, CI = 1.96, SE = 1, SD = 1)
  err_label <- switch(error_type, CI = "95% CI", SE = "SE", SD = "SD")

  for (i in seq_along(cat_paramcd_list)) {
    adlb_param <- dt(cat_paramcd_list[[i]])

    x_label <- ""
    y_label <- paste0(
      "Mean Change From Baseline (",
      err_label,
      ")\n",
      stringr::str_wrap(cat_paramcd_lbls[[i]], width = 30)
    )
    param_title <- paste0("Laboratory test: ", cat_paramcd_title[[i]])

    split_active <- exists("timepoint_split") && !is.null(timepoint_split) && length(timepoint_split) > 0
    split_indices <- if (split_active) seq_along(timepoint_split) else 1

    if (split_active) {
      ymin_l <- floor(min(adlb_param$mean - err_mult * adlb_param[[err_col]], na.rm = TRUE) - 0.5)
      ymax_l <- ceiling(max(adlb_param$mean + err_mult * adlb_param[[err_col]], na.rm = TRUE)) + 0.5
      for (j in split_indices) {
        assign(
          paste0("adlb_param_", j),
          adlb_param |>
            filter(AVISIT %in% timepoint_split[[j]]) |>
            mutate(ymin_l = ymin_l, ymax_l = ymax_l)
        )
      }
    } else {
      adlb_param_1 <- adlb_param
    }

    for (k in split_indices) {
      assign(paste0("pt_param_", k), g_line(df = get(paste0("adlb_param_", k))))
    }

    for (pn in split_indices) {
      pname_param <- paste0(tolower(tblid), "_", cat_suffix, "_", tolower(cat_paramcd_list[[i]]), "_", pn, ".png")

      png(
        write_path(opath, pname_param),
        width = 22,
        height = 14,
        units = "cm",
        res = 300,
        type = "cairo"
      )
      print(write_path(opath, pname_param))
      print(get(paste0("pt_param_", pn)))
      dev.off()

      png_list <- c(png_list, pname_param)
    }
  }

  tidytlg::gentlg(
    tlf = "g",
    plotnames = write_path(opath, png_list),
    plotwidth = 8,
    orientation = "landscape",
    opath = write_path(opath),
    file = paste0(tblid, cat_suffix),
    title = title_footer$title,
    footers = title_footer$main_footer
  )
}
