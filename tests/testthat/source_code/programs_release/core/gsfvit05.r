###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsfvit05.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create gsfvit05: Subjects Remaining in Study at Each Visit
##                            by Availability of Vital Sign Results
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     advs
## Output:                    gsfvit05.rtf
## Remarks:                   Include parameters: Systolic BP, Diastolic BP
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

tblid <- "GSFVIT05"

# Safety population flag (default=SAFFL).
popfl <- "SAFFL"

# Actual treatment variable (default=TRT01A).
trtvar <- "TRT01A"

# List of Prameters to be included
paramcd_list <- c("SYSBP", "DIABP")

# Parameter Labels list as per paramcd_list order
paramcd_lbls <- c(
  "Systolic Blood Pressure (mmHg)",
  "Diastolic Blood Pressure (mmHg)"
)

# Plot title list for each parameter as per paramcd_list order
paramcd_plt_title <- paramcd_lbls # User has to Assign Labels if it is different from paramcd_lbls

# Split Timepoint into groups if not fit into single page
## User defined timepoint
timepoint_split <- list(
  time_point_1 = c("Baseline", "Cycle 02", "Cycle 03", "Cycle 04", "Cycle 05", "Cycle 06", "Cycle 07", "Cycle 08"),
  time_point_2 = c("Cycle 09", "Cycle 10", "Cycle 11", "Cycle 12", "Cycle 13", "Cycle 15", "Cycle 17"),
  time_point_3 = c("Cycle 19", "Cycle 21", "Cycle 23", "Cycle 25", "Cycle 29", "EOT")
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
    filter(toupper(.data[[popfl]]) == "Y") |>
    filter(!grepl("Unscheduled", AVISIT)) |>
    # to add levels
    mutate(
      !!rlang::sym(trtvar) := factor(
        .data[[trtvar]],
        levels = c("Xanomeline Low Dose", "Xanomeline High Dose", "Placebo")
      ),
      AVISIT = factor(
        .data[['AVISIT']],
        levels = unique(.data[['AVISIT']])[order(unique(.data[['AVISITN']]))]
      )
    ) |>
    mutate(
      # shorten "End of Treatment" as EOT
      AVISIT = forcats::fct_recode(
        AVISIT,
        "EOT" = "End Of Treatment"
      ),

      # add abbreviation to trtvar
      !!rlang::sym(trtvar) := forcats::fct_recode(
        .data[[trtvar]],
        "Xanomeline (Low)" = "Xanomeline Low Dose",
        "Xanomeline (High)" = "Xanomeline High Dose",
        "Placebo (PBO)" = "Placebo"
      )
    )

  #### perform check on unique record per subject/param/timepoint
  check_unique <- df |>
    group_by(USUBJID, PARAMCD, .data[[trtvar]], AVISITN, AVISIT) |>
    mutate(n_recsub = n()) |>
    filter(n_recsub > 1)

  if (nrow(check_unique) > 0) {
    stop(
      "Your input dataset needs extra attention, as some subjects have more than one record per parameter/visit"
    )
  }

  df_1 <- df |>
    filter(PARAMCD == param_cd & ANL02FL == "Y") |>
    group_by(.data[[trtvar]], AVISITN, AVISIT) |>
    # subjects with existing data for parameter: AVAL^=.
    # assume each subject collecting only one record in each time point
    summarize(n = length(!is.na(AVAL))) |>
    mutate(type = "with_data") |>
    ungroup()

  df_2 <- df |>
    group_by(.data[[trtvar]], AVISITN, AVISIT) |>
    # subjects remaining with missing data for parameter
    # assume subjects collecting any vital sign parameters in each time points remaining in study
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
  adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) |>
    df_na() |>
    filter(toupper(.data[[popfl]]) == "Y") |>
    mutate(
      !!rlang::sym(trtvar) := forcats::fct_recode(
        .data[[trtvar]],
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
    mutate(trt_rev = forcats::fct_relevel(.data[[trtvar]], rev(levels(.data[[trtvar]]))))

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
      data = df$data[[2]],
      aes(
        x = .data$visit + barwidth,
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
      breaks = c("Xanomeline (Low)", "Xanomeline (High)", "Placebo (PBO)")
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
      title = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
    ) +
    labs(
      title = "Number of Subjects With Data (Solid Bar)/Number of Subjects Remaining in Study (Bar Height)"
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
  advs_param <- dt(paramcd_list[[i]])
  advs_param_pt <- advs_param[[1]]
  advs_param_tb <- advs_param[[2]]

  # define parameters for plotting:
  # present treatment group in y axis table from top to bottom (reverse)
  table_text <- c("PBO", "Xan High", "Xan Low")

  x_label <- ""
  y_label <- "Percent of Subjects"

  barwidth <- 0.3

  # Define title for the plot
  param_title <- paste0("Vital Sign: ", paramcd_plt_title[[i]])

  split_active <- exists("timepoint_split") && !is.null(timepoint_split) && length(timepoint_split) > 0
  split_indices <- if (split_active) seq_along(timepoint_split) else 1

  if (split_active) {
    for (j in split_indices) {
      tp <- timepoint_split[[j]]
      assign(
        paste0("advs_param_pt_", j),
        advs_param_pt |>
          filter(AVISIT %in% tp) |>
          select(-AVISIT) |>
          tidyr::nest(data = c(n, type, visit))
      )
      assign(paste0("advs_param_tb_", j), advs_param_tb |> filter(AVISIT %in% tp))
    }
  } else {
    advs_param_pt_1 <- advs_param_pt |>
      select(-AVISIT) |>
      tidyr::nest(data = c(n, type, visit))
    advs_param_tb_1 <- advs_param_tb
  }

  for (k in split_indices) {
    tp_k <- if (split_active) timepoint_split[[k]] else unique(advs_param_tb$AVISIT)
    offset <- if (split_active) sum(lengths(timepoint_split[seq_len(k - 1)])) else 0
    n_tp <- length(tp_k)
    x_breaks <- seq(offset + 1, offset + n_tp) + barwidth
    x_labels <- tp_k
    assign(
      paste0("pt_param_", k),
      g_bar(df = get(paste0("advs_param_pt_", k)), df2 = get(paste0("advs_param_tb_", k)))
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
