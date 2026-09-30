###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsiex01.r
## R version:                 4.5.2
## junco version:             0.1.3
## Short Description:         Program to create gsiex01: Duration of Treatment
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adexsum
## Output:                    gsiex01.rtf
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

tblid <- "GSIEX01"
trtvar <- "TRT01A"
popfl <- "SAFFL"

################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

# reading data

adexsum <- haven::read_sas(envsetup::read_path(a_in, "adexsum.sas7bdat")) |>
  df_na() |>
  filter(PARAMCD == "TRTDURM", !!rlang::sym(popfl) == "Y") |>
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c(
        "Xanomeline Low Dose",
        "Xanomeline High Dose",
        "Placebo"
      )
    )
  ) |>
  mutate(
    cat = forcats::fct_reorder(AVALCAT1, AVALCA1N),
    # abbreviate treatment group
    trt_abb = factor(case_when(
      .data[[trtvar]] == "Xanomeline Low Dose" ~ "Xan Low",
      .data[[trtvar]] == "Xanomeline High Dose" ~ "Xan High",
      .data[[trtvar]] == "Placebo" ~ "PBO"
    ))
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

# drop levels not associated with TRTDURM
adexsum$cat <- droplevels(adexsum$cat)

# calculate percentage of count per each treatment group
df1 <- adexsum |>
  group_by(cat, !!rlang::sym(trtvar), trt_abb) |>
  summarise(n = n()) |>
  ungroup() |>
  # fill in missing values for time points which don't have data
  tidyr::complete(cat, tidyr::nesting(!!rlang::sym(trtvar), trt_abb), fill = list(n = 0))

df2 <- df1 |>
  group_by(!!rlang::sym(trtvar)) |>
  summarise(total = sum(n)) |>
  ungroup()

df3 <- df1 |>
  left_join(df2, by = trtvar) |>
  mutate(
    pern = tidytlg::roundSAS((n / total) * 100, digits = 1),
    perc = tidytlg::roundSAS((n / total) * 100, digits = 1, as_char = TRUE)
  )

# check month group having no subjects in any treatment group

ck <- df3 |>
  group_by(cat) |>
  summarise(total = sum(pern)) |>
  filter(total == 0)

df <- anti_join(df3, ck, by = "cat")

# split data in two graphics
time_point_1 <- c(
  "0 to <3 months",
  "3 to <6 months",
  "6 to <9 months",
  "9 to <12 months",
  "12 to <15 months",
  "15 to <18 months",
  "18 to <21 months"
)

time_point_2 <- c(
  "21 to <24 months",
  "24 to <27 months",
  "27 to <30 months",
  "30 to <33 months",
  "33 to <36 months",
  "36 to <39 months"
)

adex_1 <- df |>
  filter(cat %in% time_point_1)

adex_2 <- df |>
  filter(cat %in% time_point_2)


################################################################################
# Generate plot:
################################################################################

# define parameters for plotting:
# assign colorblind friendly palette: black(Xan Low), orange(Xan High), dark blue(PBO)
cbbPalette <- c("#000000", "#E69F00", "#0072B2")

x_label <- " "
y_label <- "Percentage of Subjects"
param_title <- "Safety Analysis Set"


# Bar plot ---------------------------------------------------

g_facet <- function(df) {
  # bar plot in facet -------------------------------------
  plot1 <- df |>
    ggplot(aes(x = .data$trt_abb, y = .data$pern, fill = .data$TRT01A)) +
    geom_col(position = position_dodge(0.5)) +
    geom_text(aes(label = perc), position = position_dodge(0.5), vjust = -0.5, family = "arial", size = 8 / .pt) +
    facet_wrap(~ .data$cat, nrow = 1, strip.position = "bottom") +

    # assign colorblind friendly palette
    scale_fill_manual(values = cbbPalette) +
    labs(x = x_label, y = y_label) +
    theme_bw() +
    theme(
      text = element_text(size = 9, color = "black", family = "arial"),
      strip.background = element_rect(fill = NA, color = "black"),
      strip.placement = "outside",
      strip.text = element_text(size = 8),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(5, "mm"),
      legend.box.background = element_rect(colour = "black", linewidth = 0.5),
      legend.text = element_text(size = 8),
      panel.spacing.x = unit(0, "line"),
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      axis.text.x = element_text(angle = 320, vjust = -1),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 9)
    ) +
    plot_annotation(
      title = param_title,
      theme = theme(plot.title = element_text(size = 10, family = "Arial", hjust = 0.5)),
      subtitle = ""
    )
}

plot1 <- g_facet(adex_1)
plot2 <- g_facet(adex_2)


################################################################################
# Create png and output file:
################################################################################

# create png file and output figure
pname_1 <- paste0(tolower(tblid), "_1", ".png")

png(
  write_path(opath, pname_1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_1)) ### print png path and name in log
print(plot1)
dev.off()

pname_2 <- paste0(tolower(tblid), "_2", ".png")

png(
  write_path(opath, pname_2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_2)) ### print png path and name in log
print(plot2)
dev.off()

if (length(title_footer$main_footer) == 0) {
  title_footer$main_footer <- NULL
}

tidytlg::gentlg(
  tlf = "g",
  plotnames = write_path(opath, c(pname_1, pname_2)),
  plotwidth = 8,
  orientation = "landscape",
  opath = write_path(opath),
  file = tblid,
  title = title_footer$title,
  footers = title_footer$main_footer
)
