###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsfvit04.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create gsfvit04: Baseline vs. Minimum On-treatment [Vital Sign Parameter]
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     advs
## Output:                    gsfvit04.rtf
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

tblid <- "GSFVIT04"

# Safety population flag (default=SAFFL).
popfl <- "SAFFL"

# Actual treatment variable (default=TRT01A).
trtvar <- "TRT01A"

# List of Prameters to be included
paramcd_list <- c("SYSBP", "DIABP")

# Parameter Labels list as per paramcd_list order
paramcd_lbls <- c("Systolic Blood Pressure (mmHg)", "Diastolic Blood Pressure (mmHg)")

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

# define function for data processing
dt <- function(param_cd) {
  advs <- haven::read_sas(envsetup::read_path(a_in, "advs.sas7bdat")) |>
    df_na() |>
    filter(PARAMCD == param_cd & toupper(.data[[popfl]]) == "Y") |>
    # Add treatment levels in the specified order for consistent factor ordering
    mutate(
      !!rlang::sym(trtvar) := factor(
        .data[[trtvar]],
        levels = c("Xanomeline Low Dose", "Xanomeline High Dose", "Placebo")
      )
    ) |>
    # add abbreviation to TRT01A
    mutate(
      !!rlang::sym(trtvar) := forcats::fct_recode(
        .data[[trtvar]],
        "Xanomeline (Low)" = "Xanomeline Low Dose",
        "Xanomeline (High)" = "Xanomeline High Dose",
        "Placebo (PBO)" = "Placebo"
      )
    )

  baseline <- advs |>
    filter(ABLFL == "Y") |>
    select(USUBJID, all_of(trtvar), PARAMCD, AVAL) |>
    rename(baseline = AVAL)

  post_baseline <- advs |>
    filter(ANL06FL == "Y" & APOBLFL == "Y") |>
    select(USUBJID, all_of(trtvar), PARAMCD, AVAL) |>
    rename(postvalue = AVAL)

  vs <- full_join(baseline, post_baseline) |>
    filter(!is.na(baseline) & !is.na(postvalue))
}

################################################################################
# Generate plot:
################################################################################

# define function for multiple calls:
g_scatter <- function(df, x_label, y_label) {
  # scatter plot ---------------------------------------------------
  plot <- ggplot(
    df,
    aes(
      x = .data$baseline,
      y = .data$postvalue,
      shape = .data[[trtvar]],
      color = .data[[trtvar]],
      linetype = .data[[trtvar]]
    )
  ) +

    # Make each dot partially transparent, with 0.5 opacity on alpha
    geom_point(alpha = 0.5, size = 1.2) +
    # alternative - apply jitter or dodge point with random noise as follows
    # geom_point(alpha = 0.75,
    #          position = position_dodge(width = 2)) +

    # Use a hollow circle, triangle, and cross as choices for shape
    scale_shape_manual(values = c(1, 2, 3)) +

    # Add linear regression line, w/o confidence region
    geom_smooth(method = lm, se = FALSE) +

    # Add dotted reference line for no increase
    geom_abline(intercept = 0, slope = 1, linetype = 3) +
    scale_x_continuous(
      limits = c(df$axis_min[1], df$axis_max[1]),
      breaks = scales::extended_breaks(n = 5),
      labels = scales::label_number(accuracy = 1)
    ) +
    scale_y_continuous(
      limits = c(df$axis_min[1], df$axis_max[1]),
      breaks = scales::extended_breaks(n = 5),
      labels = scales::label_number(accuracy = 1)
    ) +

    # assign colorblind friendly palette
    scale_color_manual(values = cbbPalette) +
    labs(x = x_label, y = y_label) +
    theme_bw() +
    theme(
      text = element_text(size = 9, color = "black", family = "Arial"),
      axis.text = element_text(size = 9, color = "black"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9, face = "bold")
    )
}

# Png files lists to add on tidytlg::gentlg()
png_list <- c()

# Creating PNG files as per paramcd_list
for (i in seq_along(paramcd_list)) {
  # generate data set for Parameters
  advs_param <- dt(paramcd_list[[i]]) |>
    mutate(
      axis_min = floor(min(.data$baseline, .data$postvalue, na.rm = TRUE) - 5),
      axis_max = ceiling(max(.data$baseline, .data$postvalue, na.rm = TRUE) + 5)
    )

  # define parameters for plotting:
  # assign colorblind friendly palette: black(Xan Low), orange(Xan High), dark blue(PBO)
  cbbPalette <- c("#000000", "#E69F00", "#0072B2")

  # Generate scatter plot
  pt_param <- g_scatter(
    df = advs_param,
    x_label = paste0("Baseline ", paramcd_lbls[[i]]),
    y_label = paste0("Minimum On-treatment ", paramcd_lbls[[i]])
  )

  ################################################################################
  # Create png and output file:
  ################################################################################

  # create png file and output figure
  pname_param_min <- paste0(tolower(tblid), "_", tolower(paramcd_list[[i]]), "_min.png")

  png(
    write_path(opath, pname_param_min),
    width = 16,
    height = 18,
    units = "cm",
    res = 300,
    type = "cairo"
  )
  print(write_path(opath, pname_param_min)) ### print png path and name in log
  print(pt_param)
  dev.off()

  png_list <- c(png_list, pname_param_min)
}

tidytlg::gentlg(
  tlf = "g",
  plotnames = write_path(opath, png_list),
  # plotwidth   = 10,
  orientation = "portrait",
  opath = write_path(opath),
  file = tblid,
  title = title_footer$title,
  footers = title_footer$main_footer
)
