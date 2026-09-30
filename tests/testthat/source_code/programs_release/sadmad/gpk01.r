###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gpk01.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Mean and [95% Confidence Interval/Standard Deviation] of [Matrix] [Active Study
##                            Agent/Analyte] Concentrations ([units]) Over Time – [SAD/MAD][Part 1]
## Disclaimer:                This script is a direct copy of the corresponding Core Standard output identifier. For
##                            SAD/MAD specific changes, refer to tsfvit02b.r, lsidm05.r, and gsfvit02.r for examples of
##                            STUDYPRT filtering, COHORT handling, treatment column structure modifications, pooled
##                            placebo derivations, combined treatment columns, dose-level updates, and other
##                            output-specific structural differences as applicable.
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adpc
## Output:                    gpk01.rtf
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

tblid <- "gpk01"

# PK population flag
popfl <- "PKFL"

# Actual treatment variable
trtvar <- "TRT01A"

# PARAMCD to filter adpc
paramcd <- "XAN"

# Note 2: x-axis label — choose one: "Planned Study Day" / "Planned Study Week" / "Planned Study Visit"
x_axis_label <- "Days"

# Note 3: y-axis scale range for semi-log — specify [value to value]
# set to NULL to use automatic scale
y_semilog_limits <- NULL

# Note 4: LLOQ [x MRD] imputation value — specify per study
lloq_mrd <- 0.01

# flag: if TRUE, time point is concatenation of AVISIT and ATPT; if FALSE, AVISIT only
use_atpt <- TRUE

################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

# define function for data processing
dt <- function() {
  adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) |>
    df_na() |>
    filter(toupper(.data[[popfl]]) == "Y") |>
    filter(.data[[trtvar]] != "Placebo") |>
    select(STUDYID, USUBJID, all_of(trtvar)) |>
    mutate(
      !!rlang::sym(trtvar) := factor(
        .data[[trtvar]],
        levels = c(
          "Xanomeline Low Dose",
          "Xanomeline High Dose"
        )
      )
    )

  adpc <- haven::read_sas(envsetup::read_path(a_in, "adpc.sas7bdat")) |>
    df_na() |>
    filter(PARAMCD == paramcd) |>
    select(STUDYID, USUBJID, AVISIT, AVISITN, any_of(c("ATPT", "ATPTN")), AVAL) |>
    inner_join(adsl, by = c("STUDYID", "USUBJID")) |>
    mutate(
      # Note 4: values below LLOQ x MRD are treated as 0
      AVAL = ifelse(AVAL < lloq_mrd, 0, AVAL),
      AVISIT_ATPT = factor(
        if (use_atpt) paste(AVISIT, ATPT, sep = ", ") else as.character(AVISIT),
        levels = {
          if (use_atpt) {
            arrange(distinct(adpc, AVISIT, ATPT, AVISITN, ATPTN), AVISITN, ATPTN) |>
              mutate(lv = paste(AVISIT, ATPT, sep = ", ")) |>
              pull(lv)
          } else {
            arrange(distinct(adpc, AVISIT, AVISITN), AVISITN) |> pull(AVISIT)
          }
        }
      ),
      TRT01A = forcats::fct_recode(
        .data[[trtvar]],
        "Xanomeline (Low)" = "Xanomeline Low Dose",
        "Xanomeline (High)" = "Xanomeline High Dose"
      )
    ) |>
    group_by(TRT01A, AVISIT_ATPT) |>
    summarise(
      mean = mean(AVAL, na.rm = TRUE),
      sd = sd(AVAL, na.rm = TRUE),
      se = sd(AVAL, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) |>
    tidyr::complete(
      TRT01A,
      tidyr::nesting(AVISIT_ATPT),
      fill = list(mean = 0.0, sd = 0.0, se = 0.0)
    ) |>
    mutate(
      # Note 6: any summary statistic below LLOQ x MRD is replaced with LLOQ x MRD
      raw_mean = mean,
      mean = pmax(mean, lloq_mrd),
      ci_lo = pmax(mean - 1.96 * se, lloq_mrd),
      ci_hi = pmax(mean + 1.96 * se, lloq_mrd)
    ) |>
    #Converting this as per PK format
    mutate(
      mean = sapply(mean, format_sigfig_j(3, format = "xx")),
      sd = sapply(sd, format_sigfig_j(3, format = "xx"))
    ) |>
    mutate(
      mean = as.numeric(mean),
      sd = as.numeric(sd)
    )
}

adpc_sum <- dt()


# Note 7: for any visit/timepoint where mean is zero, mark its x-axis label with '<0.xxx'
n_dec <- nchar(sub(".*\\.", "", as.character(lloq_mrd)))
lloq_label <- paste0("<", formatC(lloq_mrd, format = "f", digits = n_dec))

zero_visits <- adpc_sum |>
  group_by(AVISIT_ATPT) |>
  summarise(any_zero = any(raw_mean == 0), .groups = "drop") |>
  filter(any_zero) |>
  pull(AVISIT_ATPT) |>
  as.character()

x_axis_labels <- levels(adpc_sum$AVISIT_ATPT)

################################################################################
# Generate plot:
################################################################################

# define function for creating the body of graphic
g_line <- function(df, scale_type = c("linear", "semilog")) {
  scale_type <- match.arg(scale_type)

  # line plot -------------------------------------------------------
  plot <- ggplot(
    df,
    aes(
      x = .data$AVISIT_ATPT,
      y = .data$mean,
      group = .data$TRT01A,
      color = .data$TRT01A,
      linetype = .data$TRT01A,
      shape = .data$TRT01A
    )
  ) +
    geom_errorbar(
      aes(ymin = .data$ci_lo, ymax = .data$ci_hi),
      width = 0.5,
      position = pd
    ) +
    geom_line(position = pd) +
    geom_point(position = pd) +
    # Note 5: horizontal dotted line at LLOQ x MRD, value noted on top left
    geom_hline(yintercept = lloq_mrd, linetype = "dotted", color = "black") +
    annotate(
      "text",
      x = 0.5,
      y = lloq_mrd,
      label = paste0("LLOQxMRD = ", lloq_mrd),
      hjust = 0,
      vjust = -1.2,
      size = 8 / .pt,
      family = "Arial"
    ) +
    scale_color_manual(values = cbbPalette) +
    scale_x_discrete(labels = x_axis_labels) +
    labs(x = paste("Planned Study", x_axis_label), y = y_label) +
    theme_bw() +
    theme(
      text = element_text(size = 9, color = "black", family = "Arial"),
      axis.text = element_text(size = 9, color = "black", family = "Arial"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title.x = element_text(face = "bold", size = 9, family = "Arial"),
      axis.title.y = element_text(face = "bold", size = 9, family = "Arial"),
      plot.title = element_text(size = 10, family = "Arial"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9, face = "bold")
    )

  # Note 3: apply linear or semi-log y-axis scale
  if (scale_type == "linear") {
    plot <- plot +
      scale_y_continuous(
        limits = y_semilog_limits,
        breaks = scales::extended_breaks(n = 10),
        labels = scales::label_number(drop0trailing = TRUE)
      ) +
      ggtitle("Linear Scale")
  } else {
    plot <- plot +
      scale_y_log10(
        limits = y_semilog_limits,
        breaks = scales::log_breaks(),
        labels = scales::label_number(drop0trailing = TRUE)
      ) +
      ggtitle("Semi-Log Scale")
  }

  plot
}


# assign colorblind friendly palette: black (Xan Low), orange (Xan High)
cbbPalette <- c("#000000", "#E69F00")

pd <- position_dodge(0.3)

y_label <- paste("Active Study Agent/\nAnalyte Concentrations (\u00b5g/mL)")

################################################################################
# Create png and output file:
################################################################################

png_list <- c()

# linear scale
pname_lin <- paste0(tolower(tblid), "_linear.png")

png(
  write_path(opath, pname_lin),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_lin))
print(g_line(adpc_sum, scale_type = "linear"))
dev.off()

png_list <- c(png_list, pname_lin)

# semi-log scale
pname_log <- paste0(tolower(tblid), "_semilog.png")

png(
  write_path(opath, pname_log),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_log))
print(g_line(adpc_sum, scale_type = "semilog"))
dev.off()

png_list <- c(png_list, pname_log)

if (length(title_footer$main_footer) == 0) {
  title_footer$main_footer <- NULL
}

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
