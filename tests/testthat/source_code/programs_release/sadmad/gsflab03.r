###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsflab03.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Cholestatic Drug-induced Liver Injury Screening Plot – [SAD/MAD][Part 1]
## Disclaimer:                This script is a direct copy of the corresponding Core Standard output identifier. For
##                            SAD/MAD specific changes, refer to tsfvit02b.r, lsidm05.r, and gsfvit02.r for examples of
##                            STUDYPRT filtering, COHORT handling, treatment column structure modifications, pooled
##                            placebo derivations, combined treatment columns, dose-level updates, and other
##                            output-specific structural differences as applicable.
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     addili
## Output:                    gsflab03.rtf
## Remarks:                   Include parameters: ALP
##
## Modification History:
## Rev #:
## Modified By:
## Reporting effort:
## Date:
## Description:
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

tblid <- "gsflab03"
use_log_scale <- FALSE

# Safety population flag (default = SAFFL).
popfl <- "SAFFL"

# Actual treatment variable (default = TRT01A).
trtvar <- "TRT01A"

################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

sas_vs_rds_check("addili", a_in)
addili <- readRDS(read_path(a_in, "addili.rds")) |>
  filter(toupper(.data[[popfl]]) == "Y", !is.na(.data[[trtvar]])) |>
  mutate(
    !!trtvar := factor(.data[[trtvar]], levels = c("Xanomeline Low Dose", "Xanomeline High Dose", "Placebo")),
    !!trtvar := forcats::fct_recode(
      .data[[trtvar]],
      "Xanomeline (High)" = "Xanomeline High Dose",
      "Xanomeline (Low)" = "Xanomeline Low Dose",
      "Placebo (PBO)" = "Placebo"
    )
  )

# Extract max
addili_max <- addili |>
  filter(ANL03FL == "Y", PARAMCD %in% c("ALP", "BILI")) |>
  select(USUBJID, all_of(trtvar), PARAMCD, R2ANRHI) |>
  tidyr::pivot_wider(names_from = PARAMCD, values_from = R2ANRHI) |>
  filter(!is.na(ALP) & !is.na(BILI))

# Data points to have red circle: PARAMCD='CDILI' and CRIT1FL='Y'.
hylaw <- addili |>
  filter(PARAMCD == "CDILI", CRIT1FL == "Y") |>
  distinct(USUBJID, .data[[trtvar]]) |>
  mutate(flag = "Y")

addili_max <- addili_max |>
  left_join(hylaw, by = c("USUBJID", trtvar))


################################################################################
# Generate plot:
################################################################################

# define parameters for plotting:
# assign colorblind friendly palette: black(Xan Low), orange(Xan High), dark blue(PBO)
cbbPalette <- c("#000000", "#E69F00", "#0072B2")


# call for Hepatocellular Drug-induced Liver Injury Screening Plot
# ALP vs BILI
x_label <- "Maximum On-treatment Alkaline Phosphatase (x ULN)"
y_label <- "Maximum On-treatment Total Bilirubin (x ULN)"


# eDISH plot ---------------------------------------------------
plot <- ggplot(
  addili_max,
  aes(x = ALP, y = BILI, color = TRT01A, shape = TRT01A)
) +

  # Make each dot partially transparent, with 0.5 opacity on alpha
  geom_point(
    alpha = 0.5,
    size = 1.2,
    # or apply dodge point with random noise to avoid overlapping
    # position = position_dodge(width = 0.05),
    na.rm = T
  ) +

  # encircling points for Hy's Law cases (Cholestatic)
  geom_point(
    data = addili_max |> filter(flag == "Y"),
    # position = position_dodge(width = 0.05),
    pch = 21,
    size = 4,
    colour = "red"
  ) +

  # Use a hollow circle, triangle, and cross as choices for shape
  scale_shape_manual(values = c(1, 2, 3)) +
  scale_x_continuous(
    trans = ifelse(use_log_scale, "log10", "identity"),
    labels = scales::label_number(accuracy = 1)
  ) +
  scale_y_continuous(
    trans = ifelse(use_log_scale, "log10", "identity"),
    labels = scales::label_number(accuracy = 1)
  ) +
  geom_hline(yintercept = 2, color = "grey", linetype = 2) +
  geom_vline(xintercept = 2, color = "grey", linetype = 2) +

  # assign colorblind friendly palette
  scale_color_manual(values = cbbPalette) +
  labs(x = x_label, y = y_label) +
  theme_bw() +
  theme(
    text = element_text(size = 9, color = "black", family = 'arial'),
    axis.text = element_text(size = 9, color = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold")
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
print(plot)
dev.off()

tidytlg::gentlg(
  tlf = "g",
  plotnames = write_path(opath, pname),
  # plotwidth   = 10,
  # plotheight = 4,
  orientation = "landscape",
  opath = write_path(opath),
  file = tblid,
  title = title_footer$title,
  footers = title_footer$main_footer
)
