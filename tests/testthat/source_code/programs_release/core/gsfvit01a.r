###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsfvit01a.r
## R version:                 4.2.1
## Short Description:         Mean and 95% Confidence Interval of Blood Pressure Over Time by Subgroup
## Author:                    Technology Solutions
## Date:                      Feb 9, 2024
## Input:
## Output:
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

tblid <- "GSFVIT01a"

################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

# define function for data processing:
dt <- function(subgrp, param_cd) {
  sas_vs_rds_check("adsl", a_in)
adsl <- readRDS(read_path(a_in, "adsl.rds")) %>%
    select(USUBJID, {{ subgrp }})

  sas_vs_rds_check("advs", a_in)
advs <- readRDS(read_path(a_in, "advs.rds")) %>%
    left_join(adsl, by = "USUBJID") %>%
    filter(PARAMCD == param_cd & SAFFL == "Y" & ANL01FL == "Y") %>%
    filter(!grepl("Unscheduled", AVISIT)) %>%
    mutate(
      # shorten "End of Treatment" as EOT
      AVISIT = forcats::fct_recode(
        AVISIT,
        "EOT" = "End Of Treatment",
        "Baseline" = "Screening"
      ),

      # add abbreviation to TRT01A
      TRT01A = forcats::fct_recode(
        TRT01A,
        "Xanomeline (High)" = "Xanomeline High Dose",
        "Xanomeline (Low)" = "Xanomeline Low Dose",
        "Placebo (PBO)" = "Placebo"
      )
    ) %>%
    group_by({{ subgrp }}, TRT01A, AVISITN, AVISIT) %>%
    summarize(
      n = n(),
      mean = mean(AVAL, na.rm = TRUE),
      sd = sd(AVAL, na.rm = TRUE),
      se = sd(AVAL, na.rm = TRUE) / sqrt(n())
    ) %>%
    ungroup() %>%
    # fill in missing values for time points which don't have data
    tidyr::complete(
      TRT01A,
      tidyr::nesting(AVISIT),
      fill = list(n = 0, mean = 0.0, sd = 0.0, se = 0.0)
    ) %>%
    # filling zero in decimal place for table display
    # mutate(meanC = sprintf('%.1f', mean)) %>%
    mutate(meanC = tidytlg::roundSAS(mean, digits = 1, as_char = TRUE)) %>%
    # reverse treatment group order: present from top to bottom in table
    mutate(trt_rev = forcats::fct_relevel(TRT01A, rev(levels(TRT01A)))) %>%
    # assing EOT as NA (avoid to connect line to last time point)
    mutate(new_avisit = na_if(AVISIT, "EOT"))

  # Include time points only for subjects at least 10% left
  avisit_n <- advs %>%
    group_by(AVISIT) %>%
    summarize(n = sum(n)) %>%
    ungroup()

  base_n <- avisit_n %>%
    filter(AVISIT == "Baseline") %>%
    select(total = n)

  avisit_n <- cbind(avisit_n, base_n) %>%
    mutate(prop = (n / total)) %>%
    # remove time points of subjects less than 10%
    filter(prop >= 0.1) %>%
    select(AVISIT)

  time_point <- as.vector(unlist(avisit_n["AVISIT"]))

  advs <- advs %>%
    filter(AVISIT %in% time_point)
}


# generate data set for Systolic blood pressure
advs_sys <- dt(AGEGR1, "SYSBP")

# split data in two graphics
time_point_1 <- c(
  "Baseline",
  "Cycle 02",
  "Cycle 03",
  "Cycle 04",
  "Cycle 05",
  "Cycle 06",
  "Cycle 07",
  "Cycle 08",
  "Cycle 09",
  "Cycle 10",
  "Cycle 11",
  "Cycle 12",
  "Cycle 13"
)

time_point_2 <- c(
  "Cycle 15",
  "Cycle 17",
  "Cycle 19",
  "Cycle 21",
  "Cycle 23",
  "Cycle 25",
  "Cycle 29",
  "EOT"
)


advs_sys_1 <- advs_sys %>%
  filter(AVISIT %in% time_point_1) %>%
  group_by(AGEGR1) %>%
  tidyr::nest()

advs_sys_2 <- advs_sys %>%
  filter(AVISIT %in% time_point_2) %>%
  group_by(AGEGR1) %>%
  tidyr::nest()


# generate data set for Diastolic blood pressure
advs_dia <- dt(AGEGR1, "DIABP")

advs_dia_1 <- advs_dia %>%
  filter(AVISIT %in% time_point_1) %>%
  group_by(AGEGR1) %>%
  tidyr::nest()

advs_dia_2 <- advs_dia %>%
  filter(AVISIT %in% time_point_2) %>%
  group_by(AGEGR1) %>%
  tidyr::nest()


################################################################################
# Generate plot:
################################################################################

# define function for creating the body of graphic:
g_line <- function(df, subgrp_title) {
  # line plot -------------------------------------------------------
  plot <- ggplot(
    df,
    aes(
      x = .data$AVISIT,
      y = .data$mean,
      group = .data$TRT01A,
      color = .data$TRT01A,
      linetype = .data$TRT01A,
      shape = .data$TRT01A
    )
  ) +
    geom_errorbar(
      aes(
        ymin = .data$mean - 1.96 * .data$se,
        ymax = .data$mean + 1.96 * .data$se
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
      breaks = scales::extended_breaks(n = 5),
      labels = scales::label_number(accuracy = 1)
    ) +

    # assign colorblind friendly palette
    scale_color_manual(values = cbbPalette) +
    labs(x = x_label, y = y_label) +
    theme_bw() +
    theme(
      text = element_text(size = 9, color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9, face = "bold")
    )

  # mean value table ----------------------------------------
  table_mean <- df %>%
    ggplot(aes(x = .data$AVISIT, y = .data$trt_rev, label = .data$meanC)) +
    geom_text(size = 2.5) +

    # abbreviate table text label
    scale_y_discrete(labels = table_text) +
    theme_bw() +
    theme(
      title = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
    ) +
    labs(title = "Mean Value")

  # number of patients table ----------------------------------------
  table_n <- df %>%
    ggplot(aes(x = .data$AVISIT, y = .data$trt_rev, label = .data$n)) +
    geom_text(size = 2.5) +

    # abbreviate table text label
    scale_y_discrete(labels = table_text) +
    theme_bw() +
    theme(
      title = element_text(size = 8, face = "bold"),
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
      title = subgrp_title,
      subtitle = param_title,
      theme = theme(
        plot.title = element_text(size = 9),
        plot.subtitle = element_text(size = 9)
      )
    )
}


# define parameters for plotting:
# present treatment group in y axis table from top to bottom (reverse)
table_text <- c("PBO", "Xan High", "Xan Low")

# assign colorblind friendly palette to treatment groups: black(Xan Low), orange(Xan High), dark blue(PBO)
cbbPalette <- c("#000000", "#E69F00", "#0072B2")

pd <- position_dodge(0.3)


# Systolic blood pressure
x_label <- ""
y_label <- "Mean Value (95% CI) \n Systolic Blood Pressure (mmHg)"
param_title <- "Blood pressure: Systolic (mmHg)"

# call for Systolic blood pressure part 1
pt_sys_1_grp_1 <- g_line(
  df = advs_sys_1$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt_sys_1_grp_2 <- g_line(
  df = advs_sys_1$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt_sys_1_grp_3 <- g_line(
  df = advs_sys_1$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# call for Systolic blood pressure part 2
pt_sys_2_grp_1 <- g_line(
  df = advs_sys_2$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt_sys_2_grp_2 <- g_line(
  df = advs_sys_2$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt_sys_2_grp_3 <- g_line(
  df = advs_sys_2$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# Diastolic blood pressure
x_label <- ""
y_label <- "Mean Value (95% CI) \n Diastolic Blood Pressure (mmHg)"
param_title <- "Blood pressure: Diastolic (mmHg)"

# call for Diastolic blood pressure part 1
pt_dia_1_grp_1 <- g_line(
  df = advs_dia_1$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt_dia_1_grp_2 <- g_line(
  df = advs_dia_1$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt_dia_1_grp_3 <- g_line(
  df = advs_dia_1$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# call for Diastolic blood pressure part 2
pt_dia_2_grp_1 <- g_line(
  df = advs_dia_2$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt_dia_2_grp_2 <- g_line(
  df = advs_dia_2$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt_dia_2_grp_3 <- g_line(
  df = advs_dia_2$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


################################################################################
# Create png and output file:
################################################################################

# create png file and output figure
# systolic BP
# subgroup 1
pname_sys_1_grp_1 <- paste0(tolower(tblid), "_sys_1_grp_1", ".png")

png(
  write_path(opath, pname_sys_1_grp_1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sys_1_grp_1)) ### print png path and name in log
print(pt_sys_1_grp_1)
dev.off()

pname_sys_2_grp_1 <- paste0(tolower(tblid), "_sys_2_grp_1", ".png")

png(
  write_path(opath, pname_sys_2_grp_1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sys_2_grp_1)) ### print png path and name in log
print(pt_sys_2_grp_1)
dev.off()

# subgroup 2
pname_sys_1_grp_2 <- paste0(tolower(tblid), "_sys_1_grp_2", ".png")

png(
  write_path(opath, pname_sys_1_grp_2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sys_1_grp_2)) ### print png path and name in log
print(pt_sys_1_grp_2)
dev.off()

pname_sys_2_grp_2 <- paste0(tolower(tblid), "_sys_2_grp_2", ".png")

png(
  write_path(opath, pname_sys_2_grp_2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sys_2_grp_2)) ### print png path and name in log
print(pt_sys_2_grp_2)
dev.off()

# subgroup 3
pname_sys_1_grp_3 <- paste0(tolower(tblid), "_sys_1_grp_3", ".png")

png(
  write_path(opath, pname_sys_1_grp_3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sys_1_grp_3)) ### print png path and name in log
print(pt_sys_1_grp_3)
dev.off()

pname_sys_2_grp_3 <- paste0(tolower(tblid), "_sys_2_grp_3", ".png")

png(
  write_path(opath, pname_sys_2_grp_3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sys_2_grp_3)) ### print png path and name in log
print(pt_sys_2_grp_3)
dev.off()


# Diastolic BP
# subgroup 1
pname_dia_1_grp_1 <- paste0(tolower(tblid), "_dia_1_grp_1", ".png")

png(
  write_path(opath, pname_dia_1_grp_1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_dia_1_grp_1)) ### print png path and name in log
print(pt_dia_1_grp_1)
dev.off()

pname_dia_2_grp_1 <- paste0(tolower(tblid), "_dia_2_grp_1", ".png")

png(
  write_path(opath, pname_dia_2_grp_1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_dia_2_grp_1)) ### print png path and name in log
print(pt_dia_2_grp_1)
dev.off()

# subgroup 2
pname_dia_1_grp_2 <- paste0(tolower(tblid), "_dia_1_grp_2", ".png")

png(
  write_path(opath, pname_dia_1_grp_2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_dia_1_grp_2)) ### print png path and name in log
print(pt_dia_1_grp_2)
dev.off()

pname_dia_2_grp_2 <- paste0(tolower(tblid), "_dia_2_grp_2", ".png")

png(
  write_path(opath, pname_dia_2_grp_2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_dia_2_grp_2)) ### print png path and name in log
print(pt_dia_2_grp_2)
dev.off()

# subgroup 3
pname_dia_1_grp_3 <- paste0(tolower(tblid), "_dia_1_grp_3", ".png")

png(
  write_path(opath, pname_dia_1_grp_3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_dia_1_grp_3)) ### print png path and name in log
print(pt_dia_1_grp_3)
dev.off()

pname_dia_2_grp_3 <- paste0(tolower(tblid), "_dia_2_grp_3", ".png")

png(
  write_path(opath, pname_dia_2_grp_3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_dia_2_grp_3)) ### print png path and name in log
print(pt_dia_2_grp_3)
dev.off()


tidytlg::gentlg(
  tlf = "g",
  plotnames = write_path(
    opath,
    c(
      pname_sys_1_grp_1,
      pname_sys_2_grp_1,
      pname_sys_1_grp_2,
      pname_sys_2_grp_2,
      pname_sys_1_grp_3,
      pname_sys_2_grp_3,
      pname_dia_1_grp_1,
      pname_dia_2_grp_1,
      pname_dia_1_grp_2,
      pname_dia_2_grp_2,
      pname_dia_1_grp_3,
      pname_dia_2_grp_3
    )
  ),
  plotwidth = 8,
  orientation = "landscape",
  opath = write_path(opath),
  file = tblid,
  title = title_footer$title,
  footers = title_footer$main_footer
)
