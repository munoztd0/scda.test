###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsflab02.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Hepatocellular Drug-induced Liver Injury Screening Plot – [SAD/MAD][Part 1]
## Disclaimer:                This script is a direct copy of the corresponding Core Standard output identifier. For
##                            SAD/MAD specific changes, refer to tsfvit02b.r, lsidm05.r, and gsfvit02.r for examples of
##                            STUDYPRT filtering, COHORT handling, treatment column structure modifications, pooled
##                            placebo derivations, combined treatment columns, dose-level updates, and other
##                            output-specific structural differences as applicable.
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     addili
## Output:                    gsflab02.rtf
## Remarks:                   Include parameters: ALT or AST
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

###############################################################################
# Script level parameters
###############################################################################

tblid <- "gsflab02"

# Safety population flag (default = SAFFL).
popfl <- "SAFFL"

# Actual treatment variable (default = TRT01A).
trtvar <- "TRT01A"

alt_ast <- "both" # Enzyme selection option: "both", "alt", or "ast"
use_log_scale <- FALSE
labels <- c("Cholestasis", "Potential Hy's Law", "Temple's Corollary")

################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

enzyme_label <- switch(alt_ast, alt = "ALT", ast = "AST", "ALT or AST")

addili <- haven::read_sas(read_path(a_in, "addili.sas7bdat")) |>
  df_na() |>
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

#Extract max
addili_max <- addili |>
  filter(ANL02FL == "Y") |>
  filter(
    (alt_ast == "both" & PARCAT1 == "ALT or AST") |
      (alt_ast == "alt" & PARAMCD == "ALT") |
      (alt_ast == "ast" & PARAMCD == "AST") |
      PARAMCD == "BILI"
  ) |>
  mutate(
    PLOT_PARAM = if_else(PARAMCD == "BILI", "BILI", "ALTAST")
  ) |>
  select(USUBJID, all_of(trtvar), PLOT_PARAM, R2ANRHI) |>
  tidyr::pivot_wider(names_from = PLOT_PARAM, values_from = R2ANRHI) |>
  filter(!is.na(ALTAST) & !is.na(BILI))

# Data points to have red circle: PARAMCD='HDILI' and CRIT1FL=Y.
hylaw <- addili |>
  filter(PARAMCD == "HDILI", CRIT1FL == "Y") |>
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
# Selected enzyme vs BILI
x_label <- paste0("Maximum On-treatment ", enzyme_label, " (x ULN)")
y_label <- "Maximum On-treatment Total Bilirubin (x ULN)"
xintercept <- 3
yintercept <- 2


# find limits
min_x <- min(addili_max$ALTAST, na.rm = TRUE)
min_y <- min(addili_max$BILI, na.rm = TRUE)
max_x <- max(addili_max$ALTAST, na.rm = TRUE)
max_y <- max(addili_max$BILI, na.rm = TRUE)

# plot limits beyond the thresholds (3 and 2)
plot_min_x <- if (use_log_scale) max(min_x, 0.1) else min(min_x, 0)
plot_min_y <- if (use_log_scale) max(min_y, 0.1) else min(min_y, 0)
plot_max_x <- max(max_x, xintercept + 2)
plot_max_y <- max(max_y, yintercept + 2)

if (use_log_scale) {
  x_left <- 10^((log10(plot_min_x) + log10(3)) / 2)
  x_right <- 10^((log10(3) + log10(plot_max_x)) / 2)
  y_bot <- 10^((log10(plot_min_y) + log10(2)) / 2)
  y_top <- 10^((log10(2) + log10(plot_max_y)) / 2)
} else {
  x_left <- (plot_min_x + 3) / 2
  x_right <- (3 + plot_max_x) / 2
  y_bot <- (plot_min_y + 2) / 2
  y_top <- (2 + plot_max_y) / 2
}

# build labels
anno_df <- data.frame(
  x_pos = c(x_left, x_right, x_right),
  y_pos = c(y_top, y_top, y_bot),
  label = labels
)

# eDISH plot ---------------------------------------------------
plot <- ggplot(
  addili_max,
  aes(x = ALTAST, y = BILI, color = .data[[trtvar]], shape = .data[[trtvar]])
) +

  # Make each dot partially transparent, with 0.5 opacity on alpha
  geom_point(
    alpha = 0.5,
    size = 1.2,
    na.rm = T
  ) +

  # encircling points for Hy's Law cases
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
  geom_hline(yintercept = yintercept, color = "grey", linetype = 2) +
  geom_vline(xintercept = xintercept, color = "grey", linetype = 2) +

  # assign colorblind friendly palette
  scale_color_manual(values = cbbPalette) +

  # annotate text in quadrants
  geom_text(
    data = anno_df,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 8 / .pt,
    family = "arial",
    color = "black"
  ) +
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
  plotheight = 4,
  orientation = "landscape",
  opath = write_path(opath),
  file = tblid,
  title = title_footer$title,
  footers = title_footer$main_footer
)
