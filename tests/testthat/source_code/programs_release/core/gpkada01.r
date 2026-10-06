###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gpkada01.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create gpkada01: [Mean/Median] [Matrix]
##                            [Active Study Agent] Concentrations ([units]) by
##                            Treatment-emergent Antibody Status
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adpc, adishum
## Output:                    gpkada01.rtf
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

tblid <- "GPKADA01"

# PK population flag
popfl <- "PKFL"

# Actual treatment variable
trtvar <- "TRT01A"

# PARAMCD to filter adpc
paramcd <- "XAN"

# Treatment label used in ADA status descriptions
treatment_label <- "Antibodies to Active Study Agent"

# Note 1: x-axis label — choose one: "Planned Study Day" / "Planned Study Week" / "Planned Study Visit"
x_axis_label <- "Planned Study Days"

# Note 2: LLOQ [x MRD] imputation value — specify per study
lloq_mrd <- 0.01

# flag: if TRUE, time point is concatenation of AVISIT and ATPT; if FALSE, AVISIT only
use_atpt <- TRUE

# Summary statistic for the line plot — choose one: "median" or "mean"
stat <- "mean"

################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

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
    mutate(
      # Note 2: values below LLOQ x MRD are treated as 0
      AVAL = ifelse(AVAL < lloq_mrd, 0, AVAL)
    )

  # ADA status: Negative = ADANTRE + AVALC=Y, Positive = ADATRE + AVALC=Y
  ada_status_levels <- c(
    paste("Positive for Treatment-emergent", treatment_label),
    paste("Negative for Treatment-emergent", treatment_label)
  )

  adishum_ada <- haven::read_sas(envsetup::read_path(a_in, "adishum.sas7bdat")) |>
    df_na() |>
    filter(PARAMCD %in% c("ADANTRE", "ADATRE") & IMEVFL == "Y") |>
    select(STUDYID, USUBJID, IMEVFL, PARAMCD, AVALC) |>
    pivot_wider(
      id_cols = c(STUDYID, USUBJID, IMEVFL),
      names_from = PARAMCD,
      values_from = AVALC,
      values_fn = dplyr::first
    ) |>
    mutate(
      ADA_STATUS = factor(
        case_when(
          ADATRE == "Y" ~ paste("Positive for Treatment-emergent", treatment_label),
          ADANTRE == "Y" ~ paste("Negative for Treatment-emergent", treatment_label),
          TRUE ~ NA_character_
        ),
        levels = ada_status_levels
      )
    )

  # join adpc with ADA status and adsl
  adpc_ada <- adpc |>
    inner_join(adishum_ada, by = c("STUDYID", "USUBJID")) |>
    inner_join(adsl, by = c("STUDYID", "USUBJID")) |>
    mutate(
      AVISIT_ATPT = factor(
        if (use_atpt) paste(AVISIT, ATPT, sep = ", ") else as.character(AVISIT),
        levels = {
          if (use_atpt) {
            arrange(distinct(adpc_ada, AVISIT, ATPT, AVISITN, ATPTN), AVISITN, ATPTN) |>
              mutate(lv = paste(AVISIT, ATPT, sep = ", ")) |>
              pull(lv)
          } else {
            arrange(distinct(adpc_ada, AVISIT, AVISITN), AVISITN) |> pull(AVISIT)
          }
        }
      ),
      ADA_STATUS = factor(as.character(ADA_STATUS), levels = ada_status_levels)
    ) |>
    select(-trtvar)

  adpc_ada
}

adpc_ada <- dt()

ada_status_levels <- levels(adpc_ada$ADA_STATUS)

# summarise median, mean, SD, n per ADA status x visit
# n = unique subjects to match titer table subject counts
adpc_sum <- adpc_ada |>
  group_by(ADA_STATUS, AVISIT_ATPT) |>
  summarise(
    n = n_distinct(USUBJID),
    mean = mean(AVAL, na.rm = TRUE),
    sd = sd(AVAL, na.rm = TRUE),
    median = median(AVAL, na.rm = TRUE),
    .groups = "drop"
  ) |>
  tidyr::complete(
    ADA_STATUS,
    tidyr::nesting(AVISIT_ATPT),
    fill = list(n = 0, mean = 0.0, sd = 0.0, median = 0.0)
  ) |>
  mutate(
    # clamp summary stats below LLOQ x MRD
    median = pmax(median, lloq_mrd),
    mean = pmax(mean, lloq_mrd)
  ) |>
  #Converting this as per PK format
  mutate(
    mean_f = sapply(mean, format_sigfig_j(3, format = "xx")),
    sd_f = sapply(sd, format_sigfig_j(3, format = "xx")),
    median_f = sapply(median, format_sigfig_j(3, format = "xx"))
  ) |>
  mutate(
    mean = as.numeric(mean_f),
    sd = as.numeric(sd_f),
    median = as.numeric(median_f)
  ) |>
  # formatted strings for below-plot tables
  mutate(
    mean_sd_fmt = paste0(mean_f, " (", sd_f, ")"),
    median_fmt = median_f,
    n_fmt = as.character(n),
    #reverse ADA_STATUS for y-axis ordering in tables (top = first level)
    ada_rev = forcats::fct_rev(ADA_STATUS)
  )

# abbreviated labels for table y-axis — symbols match graph shapes (○ = Positive, □ = Negative)
# ada_rev = fct_rev(Positive, Negative) = (Negative, Positive) — so map accordingly
table_text <- setNames(
  c("□ N -", "○ N +"),
  rev(ada_status_levels)
)

################################################################################
# Generate plot:
################################################################################

g_ada <- function(df_sum) {
  # line plot
  plot <- ggplot(
    df_sum,
    aes(
      x = .data$AVISIT_ATPT,
      y = .data[[stat]],
      group = .data$ADA_STATUS,
      color = .data$ADA_STATUS,
      linetype = .data$ADA_STATUS,
      shape = .data$ADA_STATUS
    )
  ) +
    geom_line(position = pd) +
    geom_point(position = pd) +
    geom_hline(yintercept = lloq_mrd, linetype = "dotted", color = "black") +
    annotate(
      "text",
      x = 0.5,
      y = lloq_mrd,
      label = paste0("LLOQxMRD = ", lloq_mrd),
      hjust = 0,
      vjust = -1.5,
      size = 8 / .pt,
      family = "Arial"
    ) +
    scale_color_manual(
      values = c("black", "black"),
      labels = c("ADA +", "ADA -")
    ) +
    scale_linetype_manual(
      values = c("solid", "solid"),
      labels = c("ADA +", "ADA -")
    ) +
    scale_shape_manual(
      values = c(1, 0), # open circle = Positive Ab, open square = Negative Ab
      labels = c("ADA +", "ADA -")
    ) +
    scale_y_continuous(
      breaks = scales::extended_breaks(n = 10),
      labels = scales::label_number(drop0trailing = TRUE)
    ) +
    scale_x_discrete() +
    labs(
      x = x_axis_label,
      y = y_label
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 9, color = "black", family = "Arial"),
      axis.text = element_text(size = 9, color = "black", family = "Arial"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title.x = element_text(face = "bold", size = 9, family = "Arial"),
      axis.title.y = element_text(face = "bold", size = 9, family = "Arial"),
      plot.title = element_text(size = 10, family = "Arial"),
      legend.position = c(0.80, 0.70),
      legend.title = element_blank(),
      legend.text = element_text(size = 9, face = "bold"),
      #plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )

  # below-plot table: N
  tbl_n <- df_sum |>
    ggplot(aes(x = .data$AVISIT_ATPT, y = .data$ada_rev, label = .data$n_fmt)) +
    geom_text(size = 9 / .pt) +
    scale_y_discrete(labels = table_text) +
    theme_bw() +
    theme(
      title = element_text(size = 9, face = "bold", family = "Arial"),
      text = element_text(family = "Arial"),
      axis.text = element_text(size = 9, face = "bold"),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    labs(title = "Number of Subjects")

  # compose: line plot + three tables stacked
  list(
    plot = plot,
    tables = tbl_n + plot_layout(heights = c(1.5))
  )
}

pd <- position_dodge(0.3)
y_label <- paste0(stringr::str_to_title(stat), " Matrix Active Study Agent\nConcentration (\u00b5g/mL)")

################################################################################
# Create png and output file:
################################################################################

png_list <- c()

out <- g_ada(adpc_sum)

pname_plot <- paste0(tolower(tblid), "_plot.png")
png(
  write_path(opath, pname_plot),
  width = 22,
  height = 14,
  units = "cm",
  res = 600,
  type = "cairo"
)
print(write_path(opath, pname_plot))
print(out$plot)
dev.off()
png_list <- c(png_list, pname_plot)

pname_tbl <- paste0(tolower(tblid), "_tables.png")
png(
  write_path(opath, pname_tbl),
  width = 40,
  height = 3,
  units = "cm",
  res = 600,
  type = "cairo"
)

print(write_path(opath, pname_tbl))
print(out$tables)
dev.off()
png_list <- c(png_list, pname_tbl)

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
