###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsflab05.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create gsflab05: Subjects Remaining in Study
##                            at Each Visit by Availability of [Liver Function Laboratory Test] Results
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adlb, adsl
## Output:                    gsflab05.rtf
## Remarks:                   Include parameters: Liver Function laboratory tests-ALT,AST,ALP,BILI,INR,CREAT
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

tblid <- "GSFLAB05"

# Safety population flag (default=SAFFL).
popfl <- "SAFFL"

# Actual treatment variable (default=TRT01A).
trtvar <- "TRT01A"

# List of Parameters to be included
paramcd_list <- c("ALT", "AST", "ALP", "BILI", "GGT", "INR", "CREAT", "GFRCRT")

# Plot title list for each parameter as per paramcd_list order
paramcd_plt_title <- c(
  "Alanine Aminotransferase (U/L)",
  "Aspartate Aminotransferase (U/L)",
  "Alkaline Phosphatase (U/L)",
  "Bilirubin (umol/L)",
  "Gamma Glutamyl Transferase",
  "Prothrombin Intl. Normalized Ratio (RATIO)",
  "Creatinine (umol/L)",
  "GFR from Creatinine"
)

# Split Timepoint into groups if not fit into single page
# default = NULL (no split)
# timepoint_split <- NULL

## user can comment this part if no split needed
# else User can define timepoints here
timepoint_split <- list(
  time_point_1 = c(
    "Baseline",
    "Cycle 02",
    "Cycle 03",
    "Cycle 04",
    "Cycle 05"
  ),
  time_point_2 = c(
    "Cycle 12",
    "Cycle 16",
    "Cycle 20",
    "Cycle 24",
    "EOT"
  )
)

################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

adlb_full <- haven::read_sas(read_path(a_in, "adlb.sas7bdat")) |>
  df_na() |>
  filter(toupper(.data[[popfl]]) == "Y" & ANL02FL == "Y")

# Filter paramcd_list to only parameters existing in data
existing_paramcd <- paramcd_list %in% unique(adlb_full$PARAMCD)
paramcd_list <- paramcd_list[existing_paramcd]
paramcd_plt_title <- paramcd_plt_title[existing_paramcd]

# define function for data processing:
dt <- function(param_cd) {
  df <- adlb_full |>
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
      # add abbreviation to TRT01A
      !!trtvar := forcats::fct_recode(
        !!sym(trtvar),
        "Xanomeline (Low)" = "Xanomeline Low Dose",
        "Xanomeline (High)" = "Xanomeline High Dose",
        "Placebo (PBO)" = "Placebo"
      )
    )

  df_1 <- df |>
    filter(PARAMCD == param_cd) |>
    group_by(.data[[trtvar]], AVISITN, AVISIT) |>
    # subjects with existing data for parameter: LBSTRESC^=.
    # assume each subject collecting only one record in each time point
    summarize(n = length(!is.na(LBSTRESC))) |>
    mutate(type = "with_data") |>
    ungroup()

  df_2 <- df |>
    group_by(.data[[trtvar]], AVISITN, AVISIT) |>
    # subjects remaining with missing data for parameter
    # assume subjects collecting any lab tests in each time points remaining in study
    summarize(n = n_distinct(USUBJID)) |>
    mutate(type = "in_trial") |>
    ungroup()

  df_3 <- df_1 |>
    select(-type) |>
    rename(e = n) |>
    right_join(df_2) |>
    # calculate subjects in trial but no data to be presented as white in bar
    mutate(n = n - e) |>
    select(-e)

  df <- bind_rows(df_1, df_3) |>
    mutate(visit = as.numeric(factor(AVISIT))) |>
    select(all_of(trtvar), AVISIT, visit, type, n) |>
    # fill in missing values for time points which don't have data
    tidyr::complete(
      .data[[trtvar]],
      tidyr::nesting(AVISIT, visit, type),
      fill = list(n = 0)
    )

  # get total subjects per each treatment group
  adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
    df_na() |>
    filter(toupper(.data[[popfl]]) == "Y") |>
    mutate(
      !!trtvar := forcats::fct_recode(
        !!sym(trtvar),
        "Xanomeline (Low)" = "Xanomeline Low Dose",
        "Xanomeline (High)" = "Xanomeline High Dose",
        "Placebo (PBO)" = "Placebo"
      )
    ) |>
    group_by(.data[[trtvar]]) |>
    summarize(total = n()) |>
    ungroup()

  df <- df |>
    left_join(adsl, by = trtvar) |>
    # merge in total subjects to calculate percent of subjects
    mutate(n = tidytlg::roundSAS(n / total * 100), digits = 0) |>
    select(-total, -digits)

  df2 <- bind_rows(df_1, df_2) |>
    mutate(visit = as.numeric(factor(AVISIT))) |>
    tidyr::pivot_wider(
      names_from = type,
      values_from = n
    ) |>
    # fill in missing values for time points which don't have data
    tidyr::complete(
      .data[[trtvar]],
      tidyr::nesting(AVISIT, visit),
      fill = list(with_data = 0, in_trial = 0)
    ) |>
    mutate(prop = paste0(with_data, "/", in_trial)) |>
    select(all_of(trtvar), AVISIT, visit, prop) |>
    # reverse treatment group order: present from top to bottom in table
    mutate(trt_rev = forcats::fct_relevel(!!sym(trtvar), rev(levels(!!sym(trtvar)))))

  list(df, df2)
}

################################################################################
# Generate plot:
################################################################################

# define function for creating the body of graphic:

g_bar <- function(df, df2) {
  plot <- ggplot() +
    geom_bar(
      data = df$data[[1]],
      aes(
        x = .data$visit,
        y = .data$n,
        fill = factor(
          ifelse(.data$type == "with_data", "Xanomeline (High)", "InTrial"),
          levels = c("InTrial", "Xanomeline (High)")
        )
      ),
      position = "stack",
      stat = "identity",
      width = barwidth,
      color = "black"
    ) +
    geom_bar(
      data = df$data[[2]],
      aes(
        x = .data$visit + barwidth,
        y = .data$n,
        fill = factor(
          ifelse(.data$type == "with_data", "Xanomeline (Low)", "InTrial"),
          levels = c("InTrial", "Xanomeline (Low)")
        )
      ),
      position = "stack",
      stat = "identity",
      width = barwidth,
      color = "black"
    ) +
    geom_bar(
      data = df$data[[3]],
      aes(
        x = .data$visit + barwidth + barwidth,
        y = .data$n,
        fill = factor(
          ifelse(.data$type == "with_data", "Placebo (PBO)", "InTrial"),
          levels = c("InTrial", "Placebo (PBO)")
        )
      ),
      position = "stack",
      stat = "identity",
      width = barwidth,
      color = "black"
    ) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = scales::extended_breaks(n = 5),
      labels = scales::label_number(accuracy = 1)
    ) +

    # assign colorblind friendly palette to treatment groups: black(Xan Low), orange(Xan High), dark blue(PBO)
    scale_fill_manual(
      values = c(
        "Xanomeline (High)" = "#000000",
        "Xanomeline (Low)" = "#E69F00",
        "Placebo (PBO)" = "#0072B2",
        "InTrial" = "#FFFFFF"
      ),
      breaks = c("Xanomeline (High)", "Xanomeline (Low)", "Placebo (PBO)")
    ) +
    labs(x = x_label, y = y_label) +
    theme_bw() +
    theme(
      text = element_text(size = 9, color = "black", family = "Arial"),
      axis.text = element_text(size = 9, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9, face = "bold")
    )

  ## get x limit from the plot so they match the table
  built_plot <- ggplot_build(plot)
  plot_limit <- built_plot$layout$panel_params[[1]]$x$limits

  df2 <- df2 |>
    mutate(visit = visit + barwidth)

  table_prop <- df2 |>
    ggplot(aes(x = .data$visit, y = .data$trt_rev, label = .data$prop)) +
    geom_text(size = 2.5) +

    # abbreviate table text label
    scale_y_discrete(labels = table_text) +
    coord_cartesian(xlim = plot_limit) +
    theme_bw() +
    theme(
      text = element_text(family = "Arial"),
      title = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 9),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
    ) +
    labs(
      title = stringr::str_wrap(
        "Number of Subjects With Data (Solid Bar)/Number of Subjects Remaining in Study (Bar Height)",
        width = 70
      )
    )

  final <- plot /
    table_prop +
    plot_layout(heights = c(7.5, 1.5)) +
    plot_annotation(
      title = param_title,
      theme = theme(plot.title = element_text(size = 10, family = "Arial"))
    )
}

# Png files lists to add on tidytlg::gentlg()
png_list <- c()

# Creating PNG files as per paramcd_list
for (i in seq_along(paramcd_list)) {
  # generate data set for Parameters
  adlb_param <- dt(paramcd_list[[i]])
  adlb_param_pt <- adlb_param[[1]]
  adlb_param_tb <- adlb_param[[2]]

  # define parameters for plotting:
  # present treatment group in y axis table from top to bottom (reverse)
  table_text <- c("PBO", "Xan High", "Xan Low")

  x_label <- ""
  y_label <- "Percent of Subjects"

  barwidth <- 0.3

  # Define title for the plot
  param_title <- paste0("Laboratory Test: ", paramcd_plt_title[[i]])

  split_active <- exists("timepoint_split") && !is.null(timepoint_split) && length(timepoint_split) > 0
  split_indices <- if (split_active) seq_along(timepoint_split) else 1

  visit_offset <- 0
  if (split_active) {
    for (j in split_indices) {
      split_visits <- timepoint_split[[j]][timepoint_split[[j]] %in% as.character(adlb_param_tb$AVISIT)]
      assign(
        paste0("adlb_param_pt_", j),
        adlb_param_pt |>
          filter(AVISIT %in% split_visits) |>
          mutate(visit = match(as.character(AVISIT), split_visits) + visit_offset) |>
          select(-AVISIT) |>
          tidyr::nest(data = c(n, type, visit))
      )
      assign(
        paste0("adlb_param_tb_", j),
        adlb_param_tb |>
          filter(AVISIT %in% split_visits) |>
          mutate(visit = match(as.character(AVISIT), split_visits) + visit_offset)
      )
      visit_offset <- visit_offset + length(split_visits)
    }
  } else {
    visits_all <- as.character(unique(adlb_param_tb$AVISIT))
    adlb_param_pt_1 <- adlb_param_pt |>
      mutate(visit = match(as.character(AVISIT), visits_all)) |>
      select(-AVISIT) |>
      tidyr::nest(data = c(n, type, visit))
    adlb_param_tb_1 <- adlb_param_tb
  }

  for (k in split_indices) {
    pt_df <- get(paste0("adlb_param_pt_", k))
    tb_df <- get(paste0("adlb_param_tb_", k))
    x_labels <- if (split_active) {
      timepoint_split[[k]][timepoint_split[[k]] %in% as.character(tb_df$AVISIT)]
    } else {
      visits_all
    }
    x_breaks <- (min(tb_df$visit) + barwidth):(max(tb_df$visit) + barwidth)
    assign(paste0("pt_param_", k), g_bar(df = pt_df, df2 = tb_df))
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
