###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsflab01a.r
## R version:                 4.2.1
## Short Description:         Mean Change From Baseline for [Laboratory Category]
##                            Laboratory Data Over Time by Subgroup
## Author:                    Technology Solutions
## Date:                      April 12, 2024
## Input:
## Output:
## Remarks:                   Include categories: GC, KF, LV, LP, HM
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

tblid <- "GSFLAB01a"

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

  sas_vs_rds_check("adlb", a_in)
adlb <- readRDS(read_path(a_in, "adlb.rds")) %>%
    left_join(adsl, by = "USUBJID") %>%
    filter(
      PARAMCD == param_cd,
      !is.na(PARCAT3),
      ANL02FL == "Y",
      SAFFL == "Y"
    ) %>%
    filter(!grepl("Unscheduled", AVISIT)) %>%
    filter(!(AVISIT %in% c("Endpoint", "Screening"))) %>%
    mutate(
      # shorten "End of Treatment" as EOT
      AVISIT = forcats::fct_recode(AVISIT, "EOT" = "End Of Treatment"),

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
      m = mean(AVAL, na.rm = TRUE),
      mean = mean(CHG, na.rm = TRUE),
      sd = sd(CHG, na.rm = TRUE),
      se = sd(CHG, na.rm = TRUE) / sqrt(n())
    ) %>%
    ungroup() %>%
    # fill in missing values for time points which don't have data
    tidyr::complete(
      TRT01A,
      tidyr::nesting(AVISIT),
      fill = list(n = 0, m = 0.0, mean = 0.0, sd = 0.0, se = 0.0)
    ) %>%
    # filling zero in decimal place for table display
    # mutate(meanC = sprintf('%.1f', mean)) %>%
    mutate(
      meanC = paste0(
        tidytlg::roundSAS(mean, digits = 1, as_char = TRUE),
        "/",
        tidytlg::roundSAS(m, digits = 1, as_char = TRUE)
      )
    ) %>%
    # reverse treatment group order: present from top to bottom in table
    mutate(trt_rev = forcats::fct_relevel(TRT01A, rev(levels(TRT01A)))) %>%
    # assing EOT as NA (avoid to connect line to last time point)
    mutate(new_avisit = na_if(AVISIT, "EOT"))

  # Include time points only for subjects at least 10% left
  avisit_n <- adlb %>%
    group_by(AVISIT) %>%
    summarize(n = sum(n)) %>%
    ungroup()

  baseline_n <- avisit_n %>%
    filter(AVISIT == "Baseline") %>%
    select(total = n)

  avisit_n <- cbind(avisit_n, baseline_n) %>%
    mutate(prop = (n / total)) %>%
    # remove time points of subjects less than 10%
    filter(prop >= 0.1) %>%
    select(AVISIT)

  time_point <- as.vector(unlist(avisit_n["AVISIT"]))

  adlb <- adlb %>%
    filter(AVISIT %in% time_point)
}


# split data in two graphics
time_point_1 <- c(
  "Baseline",
  "Cycle 02",
  "Cycle 03",
  "Cycle 04",
  "Cycle 05",
  "Cycle 06"
)

time_point_2 <- c(
  "Cycle 07",
  "Cycle 08",
  "Cycle 09",
  "Cycle 10",
  "Cycle 11",
  "Cycle 12",
  "Cycle 13"
)

time_point_3 <- c(
  "Cycle 15",
  "Cycle 17",
  "Cycle 19",
  "Cycle 21",
  "Cycle 23",
  "Cycle 25",
  "Cycle 29",
  "EOT"
)


# list of laboratory tests
lb_param <- c(
  "SODIUM",
  "K",
  "GLUC",
  "CA",
  "PROT",
  "ALB",
  "CREAT",
  "ALP",
  "ALT",
  "AST",
  "BILI",
  "LDL",
  "WBC",
  "HGB",
  "PLAT",
  "NEUT"
)


# generate split data sets

adlb1 <- list()
adlb2 <- list()
adlb3 <- list()

for (param in lb_param) {
  df <- dt(param_cd = param, subgrp = AGEGR1)

  adlb1[[param]] <- df %>%
    filter(AVISIT %in% time_point_1) %>%
    group_by(AGEGR1) %>%
    tidyr::nest()

  adlb2[[param]] <- df %>%
    filter(AVISIT %in% time_point_2) %>%
    group_by(AGEGR1) %>%
    tidyr::nest()

  adlb3[[param]] <- df %>%
    filter(AVISIT %in% time_point_3) %>%
    group_by(AGEGR1) %>%
    tidyr::nest()
}

# generate data sets for below of parameters
# only 4 time points - baseline, cycle 12, cycle 25, cycle 29
adlb_chol <- dt(param_cd = "CHOL", subgrp = AGEGR1) %>%
  group_by(AGEGR1) %>%
  tidyr::nest()

adlb_hdl <- dt(param_cd = "HDL", subgrp = AGEGR1) %>%
  group_by(AGEGR1) %>%
  tidyr::nest()

adlb_trig <- dt(param_cd = "TRIG", subgrp = AGEGR1) %>%
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
      labels = scales::label_number(accuracy = 0.1)
    ) +
    geom_hline(yintercept = 0, color = "grey", linetype = 2) +

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
    labs(title = "Change From Baseline/Mean Value")

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

x_label <- ""

# General chemistry *************************
# Sodium (mmol/L) - SODIUM ---------
y_label <- "Mean Change From Baseline (95% CI) \n Sodium (mmol/L)"
param_title <- "Laboratory test: Sodium (mmol/L) \n  \n"

pt1_sodium_g1 <- g_line(
  df = adlb1$SODIUM$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt1_sodium_g2 <- g_line(
  df = adlb1$SODIUM$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt1_sodium_g3 <- g_line(
  df = adlb1$SODIUM$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt2_sodium_g1 <- g_line(
  df = adlb2$SODIUM$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt2_sodium_g2 <- g_line(
  df = adlb2$SODIUM$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt2_sodium_g3 <- g_line(
  df = adlb2$SODIUM$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt3_sodium_g1 <- g_line(
  df = adlb3$SODIUM$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt3_sodium_g2 <- g_line(
  df = adlb3$SODIUM$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt3_sodium_g3 <- g_line(
  df = adlb3$SODIUM$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# Potassium (mmol/L) - K ------------
y_label <- "Mean Change From Baseline (95% CI) \n Potassium (mmol/L)"
param_title <- "Laboratory test: Potassium (mmol/L) \n \n "

pt1_k_g1 <- g_line(
  df = adlb1$K$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt1_k_g2 <- g_line(
  df = adlb1$K$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt1_k_g3 <- g_line(
  df = adlb1$K$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt2_k_g1 <- g_line(
  df = adlb2$K$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt2_k_g2 <- g_line(
  df = adlb2$K$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt2_k_g3 <- g_line(
  df = adlb2$K$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt3_k_g1 <- g_line(
  df = adlb3$K$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt3_k_g2 <- g_line(
  df = adlb3$K$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt3_k_g3 <- g_line(
  df = adlb3$K$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# Glucose (mmol/L) - GLUC ------------
y_label <- "Mean Change From Baseline (95% CI) \n Glucose (mmol/L)"
param_title <- "Laboratory test: Glucose (mmol/L) \n \n "

# too few subjects
# no call setup

# Calcium (mmol/L) - CA ------------
y_label <- "Mean Change From Baseline (95% CI) \n Calcium (mmol/L)"
param_title <- "Laboratory test: Calcium (mmol/L) \n \n "

pt1_ca_g1 <- g_line(
  df = adlb1$CA$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt1_ca_g2 <- g_line(
  df = adlb1$CA$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt1_ca_g3 <- g_line(
  df = adlb1$CA$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt2_ca_g1 <- g_line(
  df = adlb2$CA$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt2_ca_g2 <- g_line(
  df = adlb2$CA$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt2_ca_g3 <- g_line(
  df = adlb2$CA$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt3_ca_g1 <- g_line(
  df = adlb3$CA$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt3_ca_g2 <- g_line(
  df = adlb3$CA$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt3_ca_g3 <- g_line(
  df = adlb3$CA$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# Protein (g/L) - PROT ------------
y_label <- "Mean Change From Baseline (95% CI) \n Protein (g/L)"
param_title <- "Laboratory test: Protein (g/L) \n \n "

pt1_prot_g1 <- g_line(
  df = adlb1$PROT$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt1_prot_g2 <- g_line(
  df = adlb1$PROT$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt1_prot_g3 <- g_line(
  df = adlb1$PROT$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt2_prot_g1 <- g_line(
  df = adlb2$PROT$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt2_prot_g2 <- g_line(
  df = adlb2$PROT$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt2_prot_g3 <- g_line(
  df = adlb2$PROT$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt3_prot_g1 <- g_line(
  df = adlb3$PROT$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt3_prot_g2 <- g_line(
  df = adlb3$PROT$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt3_prot_g3 <- g_line(
  df = adlb3$PROT$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# Albumin (g/L) - ALB -------------
y_label <- "Mean Change From Baseline (95% CI) \n Albumin (g/L)"
param_title <- "Laboratory test: Albumin (g/L) \n \n "

# too few subjects
# no call setup

# Kidney function ****************************
# Creatinine (umol/L) - CREAT -------------
y_label <- "Mean Change From Baseline (95% CI) \n Creatinine (umol/L)"
param_title <- "Laboratory test: Creatinine (umol/L) \n \n "

pt1_creat_g1 <- g_line(
  df = adlb1$CREAT$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt1_creat_g2 <- g_line(
  df = adlb1$CREAT$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt1_creat_g3 <- g_line(
  df = adlb1$CREAT$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt2_creat_g1 <- g_line(
  df = adlb2$CREAT$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt2_creat_g2 <- g_line(
  df = adlb2$CREAT$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt2_creat_g3 <- g_line(
  df = adlb2$CREAT$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt3_creat_g1 <- g_line(
  df = adlb3$CREAT$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt3_creat_g2 <- g_line(
  df = adlb3$CREAT$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt3_creat_g3 <- g_line(
  df = adlb3$CREAT$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# Liver biochemistry *************************
# Alkaline Phosphatase (U/L) - ALP ---------
y_label <- "Mean Change From Baseline (95% CI) \n Alkaline Phosphatase (U/L)"
param_title <- "Laboratory test: Alkaline Phosphatase (U/L) \n \n "

pt1_alp_g1 <- g_line(
  df = adlb1$ALP$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt1_alp_g2 <- g_line(
  df = adlb1$ALP$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt1_alp_g3 <- g_line(
  df = adlb1$ALP$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt2_alp_g1 <- g_line(
  df = adlb2$ALP$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt2_alp_g2 <- g_line(
  df = adlb2$ALP$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt2_alp_g3 <- g_line(
  df = adlb2$ALP$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt3_alp_g1 <- g_line(
  df = adlb3$ALP$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt3_alp_g2 <- g_line(
  df = adlb3$ALP$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt3_alp_g3 <- g_line(
  df = adlb3$ALP$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# Alanine Aminotransferase (U/L) - ALT -----
y_label <- "Mean Change From Baseline (95% CI) \n Alanine Aminotransferase (U/L)"
param_title <- "Laboratory test: Alanine Aminotransferase (U/L) \n \n "

pt1_alt_g1 <- g_line(
  df = adlb1$ALT$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt1_alt_g2 <- g_line(
  df = adlb1$ALT$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt1_alt_g3 <- g_line(
  df = adlb1$ALT$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt2_alt_g1 <- g_line(
  df = adlb2$ALT$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt2_alt_g2 <- g_line(
  df = adlb2$ALT$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt2_alt_g3 <- g_line(
  df = adlb2$ALT$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt3_alt_g1 <- g_line(
  df = adlb3$ALT$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt3_alt_g2 <- g_line(
  df = adlb3$ALT$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt3_alt_g3 <- g_line(
  df = adlb3$ALT$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# Aspartate Aminotransferase (U/L) - AST ---
y_label <- "Mean Change From Baseline (95% CI) \n Aspartate Aminotransferase (U/L)"
param_title <- "Laboratory test: Aspartate Aminotransferase (U/L) \n \n  "

pt1_ast_g1 <- g_line(
  df = adlb1$AST$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt1_ast_g2 <- g_line(
  df = adlb1$AST$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt1_ast_g3 <- g_line(
  df = adlb1$AST$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt2_ast_g1 <- g_line(
  df = adlb2$AST$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt2_ast_g2 <- g_line(
  df = adlb2$AST$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt2_ast_g3 <- g_line(
  df = adlb2$AST$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt3_ast_g1 <- g_line(
  df = adlb3$AST$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt3_ast_g2 <- g_line(
  df = adlb3$AST$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt3_ast_g3 <- g_line(
  df = adlb3$AST$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# Bilirubin (umol/L) - BILI ----------------
y_label <- "Mean Change From Baseline (95% CI) \n Bilirubin (umol/L)"
param_title <- "Laboratory test: Bilirubin (umol/L) \n \n "

pt1_bili_g1 <- g_line(
  df = adlb1$BILI$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt1_bili_g2 <- g_line(
  df = adlb1$BILI$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt1_bili_g3 <- g_line(
  df = adlb1$BILI$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt2_bili_g1 <- g_line(
  df = adlb2$BILI$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt2_bili_g2 <- g_line(
  df = adlb2$BILI$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt2_bili_g3 <- g_line(
  df = adlb2$BILI$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt3_bili_g1 <- g_line(
  df = adlb3$BILI$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt3_bili_g2 <- g_line(
  df = adlb3$BILI$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt3_bili_g3 <- g_line(
  df = adlb3$BILI$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# Lipids *************************************
# Cholesterol (mmol/L) - CHOL -------------
y_label <- "Mean Change From Baseline (95% CI) \n Cholesterol (mmol/L)"
param_title <- "Laboratory test: Cholesterol (mmol/L) \n \n "

pt_chol_g1 <- g_line(
  df = adlb_chol$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt_chol_g2 <- g_line(
  df = adlb_chol$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt_chol_g3 <- g_line(
  df = adlb_chol$data[[3]],
  subgrp_title = "Subgroup: Age >=75  years"
)

# HDL Cholesterol (mmol/L) - HDL ----------
y_label <- "Mean Change From Baseline (95% CI) \n HDL Cholesterol (mmol/L)"
param_title <- "Laboratory test: HDL Cholesterol (mmol/L) \n \n "

pt_hdl_g1 <- g_line(
  df = adlb_hdl$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt_hdl_g2 <- g_line(
  df = adlb_hdl$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt_hdl_g3 <- g_line(
  df = adlb_hdl$data[[3]],
  subgrp_title = "Subgroup: Age >=75  years"
)

# LDL Cholesterol (mmol/L) - LDL ----------
y_label <- "Mean Change From Baseline (95% CI) \n LDL Cholesterol (mmol/L)"
param_title <- "Laboratory test: LDL Cholesterol (mmol/L) \n \n"

# too few subjects
# no call setup

# Triglycerides (mmol/L) - TRIG -----------
y_label <- "Mean Change From Baseline (95% CI) \n Triglycerides (mmol/L)"
param_title <- "Laboratory test: Triglycerides (mmol/L) \n \n "

pt_trig_g1 <- g_line(
  df = adlb_trig$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt_trig_g2 <- g_line(
  df = adlb_trig$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt_trig_g3 <- g_line(
  df = adlb_trig$data[[3]],
  subgrp_title = "Subgroup: Age >=75  years"
)

# Complete blood count ************************
# Leukocytes (x10E9/L) - WBC -------------
y_label <- "Mean Change From Baseline (95% CI) \n Leukocytes (x10E9/L)"
param_title <- "Laboratory test: Leukocytes (x10E9/L) \n \n "

pt1_wbc_g1 <- g_line(
  df = adlb1$WBC$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt1_wbc_g2 <- g_line(
  df = adlb1$WBC$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt1_wbc_g3 <- g_line(
  df = adlb1$WBC$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt2_wbc_g1 <- g_line(
  df = adlb2$WBC$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt2_wbc_g2 <- g_line(
  df = adlb2$WBC$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt2_wbc_g3 <- g_line(
  df = adlb2$WBC$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt3_wbc_g1 <- g_line(
  df = adlb3$WBC$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt3_wbc_g2 <- g_line(
  df = adlb3$WBC$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt3_wbc_g3 <- g_line(
  df = adlb3$WBC$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# Hemoglobin (g/L) - HGB -----------------
y_label <- "Mean Change From Baseline (95% CI) \n Hemoglobin (g/L)"
param_title <- "Laboratory test: Hemoglobin (g/L) \n \n  "

pt1_hgb_g1 <- g_line(
  df = adlb1$HGB$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt1_hgb_g2 <- g_line(
  df = adlb1$HGB$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt1_hgb_g3 <- g_line(
  df = adlb1$HGB$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt2_hgb_g1 <- g_line(
  df = adlb2$HGB$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt2_hgb_g2 <- g_line(
  df = adlb2$HGB$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt2_hgb_g3 <- g_line(
  df = adlb2$HGB$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt3_hgb_g1 <- g_line(
  df = adlb3$HGB$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt3_hgb_g2 <- g_line(
  df = adlb3$HGB$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt3_hgb_g3 <- g_line(
  df = adlb3$HGB$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# Platelets (x10E9/L) - PLAT -------------
y_label <- "Mean Change From Baseline (95% CI) \n Platelets (x10E9/L)"
param_title <- "Laboratory test: Platelets (x10E9/L) \n \n "

pt1_plat_g1 <- g_line(
  df = adlb1$PLAT$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt1_plat_g2 <- g_line(
  df = adlb1$PLAT$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt1_plat_g3 <- g_line(
  df = adlb1$PLAT$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt2_plat_g1 <- g_line(
  df = adlb2$PLAT$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt2_plat_g2 <- g_line(
  df = adlb2$PLAT$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt2_plat_g3 <- g_line(
  df = adlb2$PLAT$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)

pt3_plat_g1 <- g_line(
  df = adlb3$PLAT$data[[1]],
  subgrp_title = "Subgroup: Age >=18 to <65 years"
)

pt3_plat_g2 <- g_line(
  df = adlb3$PLAT$data[[2]],
  subgrp_title = "Subgroup: Age >=65 to <75 years"
)

pt3_plat_g3 <- g_line(
  df = adlb3$PLAT$data[[3]],
  subgrp_title = "Subgroup: Age >=75 years"
)


# WBC differential ***************************
# Neutrophils (x10E9/L) - NEUT -----------
y_label <- "Mean Change From Baseline (95% CI) \n Neutrophils (x10E9/L)"
param_title <- "Laboratory test: Neutrophils (x10E9/L) \n \n "

# too few subjects
# no call setup

################################################################################
# Create png and output file:
################################################################################

# create png file and output figure
pname_sodium_1_g1 <- paste0(tolower(tblid), "_sodium_1_g1", ".png")

png(
  write_path(opath, pname_sodium_1_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sodium_1_g1)) ### print png path and name in log
print(pt1_sodium_g1)
dev.off()

pname_sodium_2_g1 <- paste0(tolower(tblid), "_sodium_2_g1", ".png")

png(
  write_path(opath, pname_sodium_2_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sodium_2_g1)) ### print png path and name in log
print(pt2_sodium_g1)
dev.off()

pname_sodium_3_g1 <- paste0(tolower(tblid), "_sodium_3_g1", ".png")

png(
  write_path(opath, pname_sodium_3_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sodium_3_g1)) ### print png path and name in log
print(pt3_sodium_g1)
dev.off()

pname_sodium_1_g2 <- paste0(tolower(tblid), "_sodium_1_g2", ".png")

png(
  write_path(opath, pname_sodium_1_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sodium_1_g2)) ### print png path and name in log
print(pt1_sodium_g2)
dev.off()

pname_sodium_2_g2 <- paste0(tolower(tblid), "_sodium_2_g2", ".png")

png(
  write_path(opath, pname_sodium_2_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sodium_2_g2)) ### print png path and name in log
print(pt2_sodium_g2)
dev.off()

pname_sodium_3_g2 <- paste0(tolower(tblid), "_sodium_3_g2", ".png")

png(
  write_path(opath, pname_sodium_3_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sodium_3_g2)) ### print png path and name in log
print(pt3_sodium_g2)
dev.off()

pname_sodium_1_g3 <- paste0(tolower(tblid), "_sodium_1_g3", ".png")

png(
  write_path(opath, pname_sodium_1_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sodium_1_g3)) ### print png path and name in log
print(pt1_sodium_g3)
dev.off()

pname_sodium_2_g3 <- paste0(tolower(tblid), "_sodium_2_g3", ".png")

png(
  write_path(opath, pname_sodium_2_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sodium_2_g3)) ### print png path and name in log
print(pt2_sodium_g3)
dev.off()

pname_sodium_3_g3 <- paste0(tolower(tblid), "_sodium_3_g3", ".png")

png(
  write_path(opath, pname_sodium_3_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_sodium_3_g3)) ### print png path and name in log
print(pt3_sodium_g3)
dev.off()

pname_k_1_g1 <- paste0(tolower(tblid), "_k_1_g1", ".png")

png(
  write_path(opath, pname_k_1_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_k_1_g1)) ### print png path and name in log
print(pt1_k_g1)
dev.off()

pname_k_2_g1 <- paste0(tolower(tblid), "_k_2_g1", ".png")

png(
  write_path(opath, pname_k_2_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_k_2_g1)) ### print png path and name in log
print(pt2_k_g1)
dev.off()

pname_k_3_g1 <- paste0(tolower(tblid), "_k_31_g1", ".png")

png(
  write_path(opath, pname_k_3_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_k_3_g1)) ### print png path and name in log
print(pt3_k_g1)
dev.off()

pname_k_1_g2 <- paste0(tolower(tblid), "_k_1_g2", ".png")

png(
  write_path(opath, pname_k_1_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_k_1_g2)) ### print png path and name in log
print(pt1_k_g2)
dev.off()

pname_k_2_g2 <- paste0(tolower(tblid), "_k_2_g2", ".png")

png(
  write_path(opath, pname_k_2_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_k_2_g2)) ### print png path and name in log
print(pt2_k_g2)
dev.off()

pname_k_3_g2 <- paste0(tolower(tblid), "_k_3_g2", ".png")

png(
  write_path(opath, pname_k_3_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_k_3_g2)) ### print png path and name in log
print(pt3_k_g2)
dev.off()

pname_k_1_g3 <- paste0(tolower(tblid), "_k_1_g3", ".png")

png(
  write_path(opath, pname_k_1_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_k_1_g3)) ### print png path and name in log
print(pt1_k_g3)
dev.off()

pname_k_2_g3 <- paste0(tolower(tblid), "_k_2_g3", ".png")

png(
  write_path(opath, pname_k_2_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_k_2_g3)) ### print png path and name in log
print(pt2_k_g3)
dev.off()

pname_k_3_g3 <- paste0(tolower(tblid), "_k_3_g3", ".png")

png(
  write_path(opath, pname_k_3_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_k_3_g3)) ### print png path and name in log
print(pt3_k_g3)
dev.off()

pname_ca_1_g1 <- paste0(tolower(tblid), "_ca_1_g1", ".png")

png(
  write_path(opath, pname_ca_1_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ca_1_g1)) ### print png path and name in log
print(pt1_ca_g1)
dev.off()

pname_ca_2_g1 <- paste0(tolower(tblid), "_ca_2_g1", ".png")

png(
  write_path(opath, pname_ca_2_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ca_2_g1)) ### print png path and name in log
print(pt2_ca_g1)
dev.off()

pname_ca_3_g1 <- paste0(tolower(tblid), "_ca_3_g1", ".png")

png(
  write_path(opath, pname_ca_3_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ca_3_g1)) ### print png path and name in log
print(pt3_ca_g1)
dev.off()

pname_ca_1_g2 <- paste0(tolower(tblid), "_ca_1_g2", ".png")

png(
  write_path(opath, pname_ca_1_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ca_1_g2)) ### print png path and name in log
print(pt1_ca_g2)
dev.off()

pname_ca_2_g2 <- paste0(tolower(tblid), "_ca_2_g2", ".png")

png(
  write_path(opath, pname_ca_2_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ca_2_g2)) ### print png path and name in log
print(pt2_ca_g2)
dev.off()

pname_ca_3_g2 <- paste0(tolower(tblid), "_ca_3_g2", ".png")

png(
  write_path(opath, pname_ca_3_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ca_3_g2)) ### print png path and name in log
print(pt3_ca_g2)
dev.off()

pname_ca_1_g3 <- paste0(tolower(tblid), "_ca_1_g3", ".png")

png(
  write_path(opath, pname_ca_1_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ca_1_g3)) ### print png path and name in log
print(pt1_ca_g3)
dev.off()

pname_ca_2_g3 <- paste0(tolower(tblid), "_ca_2_g3", ".png")

png(
  write_path(opath, pname_ca_2_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ca_2_g3)) ### print png path and name in log
print(pt2_ca_g3)
dev.off()

pname_ca_3_g3 <- paste0(tolower(tblid), "_ca_3_g3", ".png")

png(
  write_path(opath, pname_ca_3_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ca_3_g3)) ### print png path and name in log
print(pt3_ca_g3)
dev.off()

pname_prot_1_g1 <- paste0(tolower(tblid), "_prot_1_g1", ".png")

png(
  write_path(opath, pname_prot_1_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_prot_1_g1)) ### print png path and name in log
print(pt1_prot_g1)
dev.off()

pname_prot_2_g1 <- paste0(tolower(tblid), "_prot_2_g1", ".png")

png(
  write_path(opath, pname_prot_2_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_prot_2_g1)) ### print png path and name in log
print(pt2_prot_g1)
dev.off()

pname_prot_3_g1 <- paste0(tolower(tblid), "_prot_3_g1", ".png")

png(
  write_path(opath, pname_prot_3_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_prot_3_g1)) ### print png path and name in log
print(pt3_prot_g1)
dev.off()

pname_prot_1_g2 <- paste0(tolower(tblid), "_prot_1_g2", ".png")

png(
  write_path(opath, pname_prot_1_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_prot_1_g2)) ### print png path and name in log
print(pt1_prot_g2)
dev.off()

pname_prot_2_g2 <- paste0(tolower(tblid), "_prot_2_g2", ".png")

png(
  write_path(opath, pname_prot_2_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_prot_2_g2)) ### print png path and name in log
print(pt2_prot_g2)
dev.off()

pname_prot_3_g2 <- paste0(tolower(tblid), "_prot_3_g2", ".png")

png(
  write_path(opath, pname_prot_3_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_prot_3_g2)) ### print png path and name in log
print(pt3_prot_g2)
dev.off()

pname_prot_1_g3 <- paste0(tolower(tblid), "_prot_1_g3", ".png")

png(
  write_path(opath, pname_prot_1_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_prot_1_g3)) ### print png path and name in log
print(pt1_prot_g3)
dev.off()

pname_prot_2_g3 <- paste0(tolower(tblid), "_prot_2_g3", ".png")

png(
  write_path(opath, pname_prot_2_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_prot_2_g3)) ### print png path and name in log
print(pt2_prot_g3)
dev.off()

pname_prot_3_g3 <- paste0(tolower(tblid), "_prot_3_g3", ".png")

png(
  write_path(opath, pname_prot_3_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_prot_3_g3)) ### print png path and name in log
print(pt3_prot_g3)
dev.off()

pname_creat_1_g1 <- paste0(tolower(tblid), "_creat_1_g1", ".png")

png(
  write_path(opath, pname_creat_1_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_creat_1_g1)) ### print png path and name in log
print(pt1_creat_g1)
dev.off()

pname_creat_2_g1 <- paste0(tolower(tblid), "_creat_2_g1", ".png")

png(
  write_path(opath, pname_creat_2_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_creat_2_g1)) ### print png path and name in log
print(pt2_creat_g1)
dev.off()

pname_creat_3_g1 <- paste0(tolower(tblid), "_creat_3_g1", ".png")

png(
  write_path(opath, pname_creat_3_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_creat_3_g1)) ### print png path and name in log
print(pt3_creat_g1)
dev.off()

pname_creat_1_g2 <- paste0(tolower(tblid), "_creat_1_g2", ".png")

png(
  write_path(opath, pname_creat_1_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_creat_1_g2)) ### print png path and name in log
print(pt1_creat_g2)
dev.off()

pname_creat_2_g2 <- paste0(tolower(tblid), "_creat_2_g2", ".png")

png(
  write_path(opath, pname_creat_2_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_creat_2_g2)) ### print png path and name in log
print(pt2_creat_g2)
dev.off()

pname_creat_3_g2 <- paste0(tolower(tblid), "_creat_3_g2", ".png")

png(
  write_path(opath, pname_creat_3_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_creat_3_g2)) ### print png path and name in log
print(pt3_creat_g2)
dev.off()

pname_creat_1_g3 <- paste0(tolower(tblid), "_creat_1_g3", ".png")

png(
  write_path(opath, pname_creat_1_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_creat_1_g3)) ### print png path and name in log
print(pt1_creat_g3)
dev.off()

pname_creat_2_g3 <- paste0(tolower(tblid), "_creat_2_g3", ".png")

png(
  write_path(opath, pname_creat_2_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_creat_2_g3)) ### print png path and name in log
print(pt2_creat_g3)
dev.off()

pname_creat_3_g3 <- paste0(tolower(tblid), "_creat_3_g3", ".png")

png(
  write_path(opath, pname_creat_3_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_creat_3_g3)) ### print png path and name in log
print(pt3_creat_g3)
dev.off()

pname_alp_1_g1 <- paste0(tolower(tblid), "_alp_1_g1", ".png")

png(
  write_path(opath, pname_alp_1_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alp_1_g1)) ### print png path and name in log
print(pt1_alp_g1)
dev.off()

pname_alp_2_g1 <- paste0(tolower(tblid), "_alp_2_g1", ".png")

png(
  write_path(opath, pname_alp_2_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alp_2_g1)) ### print png path and name in log
print(pt2_alp_g1)
dev.off()

pname_alp_3_g1 <- paste0(tolower(tblid), "_alp_3_g1", ".png")

png(
  write_path(opath, pname_alp_3_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alp_3_g1)) ### print png path and name in log
print(pt3_alp_g1)
dev.off()

pname_alp_1_g2 <- paste0(tolower(tblid), "_alp_1_g2", ".png")

png(
  write_path(opath, pname_alp_1_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alp_1_g2)) ### print png path and name in log
print(pt1_alp_g2)
dev.off()

pname_alp_2_g2 <- paste0(tolower(tblid), "_alp_2_g2", ".png")

png(
  write_path(opath, pname_alp_2_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alp_2_g2)) ### print png path and name in log
print(pt2_alp_g2)
dev.off()

pname_alp_3_g2 <- paste0(tolower(tblid), "_alp_3_g2", ".png")

png(
  write_path(opath, pname_alp_3_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alp_3_g2)) ### print png path and name in log
print(pt3_alp_g2)
dev.off()

pname_alp_1_g3 <- paste0(tolower(tblid), "_alp_1_g3", ".png")

png(
  write_path(opath, pname_alp_1_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alp_1_g3)) ### print png path and name in log
print(pt1_alp_g3)
dev.off()

pname_alp_2_g3 <- paste0(tolower(tblid), "_alp_2_g3", ".png")

png(
  write_path(opath, pname_alp_2_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alp_2_g3)) ### print png path and name in log
print(pt2_alp_g3)
dev.off()

pname_alp_3_g3 <- paste0(tolower(tblid), "_alp_3_g3", ".png")

png(
  write_path(opath, pname_alp_3_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alp_3_g3)) ### print png path and name in log
print(pt3_alp_g3)
dev.off()

pname_alt_1_g1 <- paste0(tolower(tblid), "_alt_1_g1", ".png")

png(
  write_path(opath, pname_alt_1_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alt_1_g1)) ### print png path and name in log
print(pt1_alt_g1)
dev.off()

pname_alt_2_g1 <- paste0(tolower(tblid), "_alt_2_g1", ".png")

png(
  write_path(opath, pname_alt_2_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alt_2_g1)) ### print png path and name in log
print(pt2_alt_g1)
dev.off()

pname_alt_3_g1 <- paste0(tolower(tblid), "_alt_3_g1", ".png")

png(
  write_path(opath, pname_alt_3_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alt_3_g1)) ### print png path and name in log
print(pt3_alt_g1)
dev.off()

pname_alt_1_g2 <- paste0(tolower(tblid), "_alt_1_g2", ".png")

png(
  write_path(opath, pname_alt_1_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alt_1_g2)) ### print png path and name in log
print(pt1_alt_g2)
dev.off()

pname_alt_2_g2 <- paste0(tolower(tblid), "_alt_2_g2", ".png")

png(
  write_path(opath, pname_alt_2_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alt_2_g2)) ### print png path and name in log
print(pt2_alt_g2)
dev.off()

pname_alt_3_g2 <- paste0(tolower(tblid), "_alt_3_g2", ".png")

png(
  write_path(opath, pname_alt_3_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alt_3_g2)) ### print png path and name in log
print(pt3_alt_g2)
dev.off()

pname_alt_1_g3 <- paste0(tolower(tblid), "_alt_1_g3", ".png")

png(
  write_path(opath, pname_alt_1_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alt_1_g3)) ### print png path and name in log
print(pt1_alt_g3)
dev.off()

pname_alt_2_g3 <- paste0(tolower(tblid), "_alt_2_g3", ".png")

png(
  write_path(opath, pname_alt_2_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alt_2_g3)) ### print png path and name in log
print(pt2_alt_g3)
dev.off()

pname_alt_3_g3 <- paste0(tolower(tblid), "_alt_3_g3", ".png")

png(
  write_path(opath, pname_alt_3_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_alt_3_g3)) ### print png path and name in log
print(pt3_alt_g3)
dev.off()

pname_ast_1_g1 <- paste0(tolower(tblid), "_ast_1_g1", ".png")

png(
  write_path(opath, pname_ast_1_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ast_1_g1)) ### print png path and name in log
print(pt1_ast_g1)
dev.off()

pname_ast_2_g1 <- paste0(tolower(tblid), "_ast_2_g1", ".png")

png(
  write_path(opath, pname_ast_2_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ast_2_g1)) ### print png path and name in log
print(pt2_ast_g1)
dev.off()

pname_ast_3_g1 <- paste0(tolower(tblid), "_ast_3_g1", ".png")

png(
  write_path(opath, pname_ast_3_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ast_3_g1)) ### print png path and name in log
print(pt3_ast_g1)
dev.off()

pname_ast_1_g2 <- paste0(tolower(tblid), "_ast_1_g2", ".png")

png(
  write_path(opath, pname_ast_1_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ast_1_g2)) ### print png path and name in log
print(pt1_ast_g2)
dev.off()

pname_ast_2_g2 <- paste0(tolower(tblid), "_ast_2_g2", ".png")

png(
  write_path(opath, pname_ast_2_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ast_2_g2)) ### print png path and name in log
print(pt2_ast_g2)
dev.off()

pname_ast_3_g2 <- paste0(tolower(tblid), "_ast_3_g2", ".png")

png(
  write_path(opath, pname_ast_3_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ast_3_g2)) ### print png path and name in log
print(pt3_ast_g2)
dev.off()

pname_ast_1_g3 <- paste0(tolower(tblid), "_ast_1_g3", ".png")

png(
  write_path(opath, pname_ast_1_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ast_1_g3)) ### print png path and name in log
print(pt1_ast_g3)
dev.off()

pname_ast_2_g3 <- paste0(tolower(tblid), "_ast_2_g3", ".png")

png(
  write_path(opath, pname_ast_2_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ast_2_g3)) ### print png path and name in log
print(pt2_ast_g3)
dev.off()

pname_ast_3_g3 <- paste0(tolower(tblid), "_ast_3_g3", ".png")

png(
  write_path(opath, pname_ast_3_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_ast_3_g3)) ### print png path and name in log
print(pt3_ast_g3)
dev.off()

pname_bili_1_g1 <- paste0(tolower(tblid), "_bili_1_g1", ".png")

png(
  write_path(opath, pname_bili_1_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_bili_1_g1)) ### print png path and name in log
print(pt1_bili_g1)
dev.off()

pname_bili_2_g1 <- paste0(tolower(tblid), "_bili_2_g1", ".png")

png(
  write_path(opath, pname_bili_2_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_bili_2_g1)) ### print png path and name in log
print(pt2_bili_g1)
dev.off()

pname_bili_3_g1 <- paste0(tolower(tblid), "_bili_3_g1", ".png")

png(
  write_path(opath, pname_bili_3_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_bili_3_g1)) ### print png path and name in log
print(pt3_bili_g1)
dev.off()

pname_bili_1_g2 <- paste0(tolower(tblid), "_bili_1_g2", ".png")

png(
  write_path(opath, pname_bili_1_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_bili_1_g2)) ### print png path and name in log
print(pt1_bili_g2)
dev.off()

pname_bili_2_g2 <- paste0(tolower(tblid), "_bili_2_g2", ".png")

png(
  write_path(opath, pname_bili_2_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_bili_2_g2)) ### print png path and name in log
print(pt2_bili_g2)
dev.off()

pname_bili_3_g2 <- paste0(tolower(tblid), "_bili_3_g2", ".png")

png(
  write_path(opath, pname_bili_3_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_bili_3_g2)) ### print png path and name in log
print(pt3_bili_g2)
dev.off()

pname_bili_1_g3 <- paste0(tolower(tblid), "_bili_1_g3", ".png")

png(
  write_path(opath, pname_bili_1_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_bili_1_g3)) ### print png path and name in log
print(pt1_bili_g3)
dev.off()

pname_bili_2_g3 <- paste0(tolower(tblid), "_bili_2_g3", ".png")

png(
  write_path(opath, pname_bili_2_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_bili_2_g3)) ### print png path and name in log
print(pt2_bili_g3)
dev.off()

pname_bili_3_g3 <- paste0(tolower(tblid), "_bili_3_g3", ".png")

png(
  write_path(opath, pname_bili_3_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_bili_3_g3)) ### print png path and name in log
print(pt3_bili_g3)
dev.off()

pname_chol_g1 <- paste0(tolower(tblid), "_chol_g1", ".png")

png(
  write_path(opath, pname_chol_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_chol_g1)) ### print png path and name in log
print(pt_chol_g1)
dev.off()

pname_chol_g2 <- paste0(tolower(tblid), "_chol_g2", ".png")

png(
  write_path(opath, pname_chol_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_chol_g2)) ### print png path and name in log
print(pt_chol_g2)
dev.off()

pname_chol_g3 <- paste0(tolower(tblid), "_chol_g3", ".png")

png(
  write_path(opath, pname_chol_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_chol_g3)) ### print png path and name in log
print(pt_chol_g3)
dev.off()

pname_hdl_g1 <- paste0(tolower(tblid), "_hdl_g1", ".png")

png(
  write_path(opath, pname_hdl_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_hdl_g1)) ### print png path and name in log
print(pt_hdl_g1)
dev.off()

pname_hdl_g2 <- paste0(tolower(tblid), "_hdl_g2", ".png")

png(
  write_path(opath, pname_hdl_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_hdl_g2)) ### print png path and name in log
print(pt_hdl_g2)
dev.off()

pname_hdl_g3 <- paste0(tolower(tblid), "_hdl_g3", ".png")

png(
  write_path(opath, pname_hdl_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_hdl_g3)) ### print png path and name in log
print(pt_hdl_g3)
dev.off()

pname_trig_g1 <- paste0(tolower(tblid), "_trig_g1", ".png")

png(
  write_path(opath, pname_trig_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_trig_g1)) ### print png path and name in log
print(pt_trig_g1)
dev.off()

pname_trig_g2 <- paste0(tolower(tblid), "_trig_g2", ".png")

png(
  write_path(opath, pname_trig_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_trig_g2)) ### print png path and name in log
print(pt_trig_g2)
dev.off()

pname_trig_g3 <- paste0(tolower(tblid), "_trig_g3", ".png")

png(
  write_path(opath, pname_trig_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_trig_g3)) ### print png path and name in log
print(pt_trig_g3)
dev.off()

pname_wbc_1_g1 <- paste0(tolower(tblid), "_wbc_1_g1", ".png")

png(
  write_path(opath, pname_wbc_1_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_wbc_1_g1)) ### print png path and name in log
print(pt1_wbc_g1)
dev.off()

pname_wbc_2_g1 <- paste0(tolower(tblid), "_wbc_2_g1", ".png")

png(
  write_path(opath, pname_wbc_2_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_wbc_2_g1)) ### print png path and name in log
print(pt2_wbc_g1)
dev.off()

pname_wbc_3_g1 <- paste0(tolower(tblid), "_wbc_3_g1", ".png")

png(
  write_path(opath, pname_wbc_3_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_wbc_3_g1)) ### print png path and name in log
print(pt3_wbc_g1)
dev.off()

pname_wbc_1_g2 <- paste0(tolower(tblid), "_wbc_1_g2", ".png")

png(
  write_path(opath, pname_wbc_1_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_wbc_1_g2)) ### print png path and name in log
print(pt1_wbc_g2)
dev.off()

pname_wbc_2_g2 <- paste0(tolower(tblid), "_wbc_2_g2", ".png")

png(
  write_path(opath, pname_wbc_2_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_wbc_2_g2)) ### print png path and name in log
print(pt2_wbc_g2)
dev.off()

pname_wbc_3_g2 <- paste0(tolower(tblid), "_wbc_3_g2", ".png")

png(
  write_path(opath, pname_wbc_3_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_wbc_3_g2)) ### print png path and name in log
print(pt3_wbc_g2)
dev.off()

pname_wbc_1_g3 <- paste0(tolower(tblid), "_wbc_1_g3", ".png")

png(
  write_path(opath, pname_wbc_1_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_wbc_1_g3)) ### print png path and name in log
print(pt1_wbc_g3)
dev.off()

pname_wbc_2_g3 <- paste0(tolower(tblid), "_wbc_2_g3", ".png")

png(
  write_path(opath, pname_wbc_2_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_wbc_2_g3)) ### print png path and name in log
print(pt2_wbc_g3)
dev.off()

pname_wbc_3_g3 <- paste0(tolower(tblid), "_wbc_3_g3", ".png")

png(
  write_path(opath, pname_wbc_3_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_wbc_3_g3)) ### print png path and name in log
print(pt3_wbc_g3)
dev.off()

pname_hgb_1_g1 <- paste0(tolower(tblid), "_hgb_1_g1", ".png")

png(
  write_path(opath, pname_hgb_1_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_hgb_1_g1)) ### print png path and name in log
print(pt1_hgb_g1)
dev.off()

pname_hgb_2_g1 <- paste0(tolower(tblid), "_hgb_2_g1", ".png")

png(
  write_path(opath, pname_hgb_2_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_hgb_2_g1)) ### print png path and name in log
print(pt2_hgb_g1)
dev.off()

pname_hgb_3_g1 <- paste0(tolower(tblid), "_hgb_3_g1", ".png")

png(
  write_path(opath, pname_hgb_3_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_hgb_3_g1)) ### print png path and name in log
print(pt3_hgb_g1)
dev.off()

pname_hgb_1_g2 <- paste0(tolower(tblid), "_hgb_1_g2", ".png")

png(
  write_path(opath, pname_hgb_1_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_hgb_1_g2)) ### print png path and name in log
print(pt1_hgb_g2)
dev.off()

pname_hgb_2_g2 <- paste0(tolower(tblid), "_hgb_2_g2", ".png")

png(
  write_path(opath, pname_hgb_2_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_hgb_2_g2)) ### print png path and name in log
print(pt2_hgb_g2)
dev.off()

pname_hgb_3_g2 <- paste0(tolower(tblid), "_hgb_3_g2", ".png")

png(
  write_path(opath, pname_hgb_3_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_hgb_3_g2)) ### print png path and name in log
print(pt3_hgb_g2)
dev.off()

pname_hgb_1_g3 <- paste0(tolower(tblid), "_hgb_1_g3", ".png")

png(
  write_path(opath, pname_hgb_1_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_hgb_1_g3)) ### print png path and name in log
print(pt1_hgb_g3)
dev.off()

pname_hgb_2_g3 <- paste0(tolower(tblid), "_hgb_2_g3", ".png")

png(
  write_path(opath, pname_hgb_2_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_hgb_2_g3)) ### print png path and name in log
print(pt2_hgb_g3)
dev.off()

pname_hgb_3_g3 <- paste0(tolower(tblid), "_hgb_3_g3", ".png")

png(
  write_path(opath, pname_hgb_3_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_hgb_3_g3)) ### print png path and name in log
print(pt3_hgb_g3)
dev.off()

pname_plat_1_g1 <- paste0(tolower(tblid), "_plat_1_g1", ".png")

png(
  write_path(opath, pname_plat_1_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_plat_1_g1)) ### print png path and name in log
print(pt1_plat_g1)
dev.off()

pname_plat_2_g1 <- paste0(tolower(tblid), "_plat_2_g1", ".png")

png(
  write_path(opath, pname_plat_2_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_plat_2_g1)) ### print png path and name in log
print(pt2_plat_g1)
dev.off()

pname_plat_3_g1 <- paste0(tolower(tblid), "_plat_3_g1", ".png")

png(
  write_path(opath, pname_plat_3_g1),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_plat_3_g1)) ### print png path and name in log
print(pt3_plat_g1)
dev.off()

pname_plat_1_g2 <- paste0(tolower(tblid), "_plat_1_g2", ".png")

png(
  write_path(opath, pname_plat_1_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_plat_1_g2)) ### print png path and name in log
print(pt1_plat_g2)
dev.off()

pname_plat_2_g2 <- paste0(tolower(tblid), "_plat_2_g2", ".png")

png(
  write_path(opath, pname_plat_2_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_plat_2_g2)) ### print png path and name in log
print(pt2_plat_g2)
dev.off()

pname_plat_3_g2 <- paste0(tolower(tblid), "_plat_3_g2", ".png")

png(
  write_path(opath, pname_plat_3_g2),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_plat_3_g2)) ### print png path and name in log
print(pt3_plat_g2)
dev.off()

pname_plat_1_g3 <- paste0(tolower(tblid), "_plat_1_g3", ".png")

png(
  write_path(opath, pname_plat_1_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_plat_1_g3)) ### print png path and name in log
print(pt1_plat_g3)
dev.off()

pname_plat_2_g3 <- paste0(tolower(tblid), "_plat_2_g3", ".png")

png(
  write_path(opath, pname_plat_2_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_plat_2_g3)) ### print png path and name in log
print(pt2_plat_g3)
dev.off()

pname_plat_3_g3 <- paste0(tolower(tblid), "_plat_3_g3", ".png")

png(
  write_path(opath, pname_plat_3_g3),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname_plat_3_g3)) ### print png path and name in log
print(pt3_plat_g3)
dev.off()

tidytlg::gentlg(
  tlf = "g",
  plotnames = write_path(
    opath,
    c(
      pname_sodium_1_g1,
      pname_sodium_2_g1,
      pname_sodium_3_g1,
      pname_sodium_1_g2,
      pname_sodium_2_g2,
      pname_sodium_3_g2,
      pname_sodium_1_g3,
      pname_sodium_3_g3,
      pname_k_1_g1,
      pname_k_2_g1,
      pname_k_3_g1,
      pname_k_1_g2,
      pname_k_2_g2,
      pname_k_3_g2,
      pname_k_1_g3,
      pname_k_2_g3,
      pname_k_3_g3,
      pname_ca_1_g1,
      pname_ca_2_g1,
      pname_ca_3_g1,
      pname_ca_1_g2,
      pname_ca_2_g2,
      pname_ca_3_g2,
      pname_ca_1_g3,
      pname_ca_2_g3,
      pname_ca_3_g3,
      pname_prot_1_g1,
      pname_prot_2_g1,
      pname_prot_3_g1,
      pname_prot_1_g2,
      pname_prot_2_g2,
      pname_prot_3_g2,
      pname_prot_1_g3,
      pname_prot_2_g3,
      pname_prot_3_g3,
      pname_creat_1_g1,
      pname_creat_2_g1,
      pname_creat_3_g1,
      pname_creat_1_g2,
      pname_creat_2_g2,
      pname_creat_3_g2,
      pname_creat_1_g3,
      pname_creat_2_g3,
      pname_creat_3_g3,
      pname_alp_1_g1,
      pname_alp_2_g1,
      pname_alp_3_g1,
      pname_alp_1_g2,
      pname_alp_2_g2,
      pname_alp_3_g2,
      pname_alp_1_g3,
      pname_alp_2_g3,
      pname_alp_3_g3,
      pname_alt_1_g1,
      pname_alt_2_g1,
      pname_alt_3_g1,
      pname_alt_1_g2,
      pname_alt_2_g2,
      pname_alt_3_g2,
      pname_alt_1_g3,
      pname_alt_2_g3,
      pname_alt_3_g3,
      pname_ast_1_g1,
      pname_ast_2_g1,
      pname_ast_3_g1,
      pname_ast_1_g2,
      pname_ast_2_g2,
      pname_ast_3_g2,
      pname_ast_1_g3,
      pname_ast_2_g3,
      pname_ast_3_g3,
      pname_bili_1_g1,
      pname_bili_2_g1,
      pname_bili_3_g1,
      pname_bili_1_g2,
      pname_bili_2_g2,
      pname_bili_3_g2,
      pname_bili_1_g3,
      pname_bili_2_g3,
      pname_bili_3_g3,
      pname_chol_g1,
      pname_chol_g2,
      pname_chol_g3,
      pname_hdl_g1,
      pname_hdl_g2,
      pname_hdl_g3,
      pname_trig_g1,
      pname_trig_g2,
      pname_trig_g3,
      pname_wbc_1_g1,
      pname_wbc_2_g1,
      pname_wbc_3_g1,
      pname_wbc_1_g2,
      pname_wbc_2_g2,
      pname_wbc_3_g2,
      pname_wbc_1_g3,
      pname_wbc_2_g3,
      pname_wbc_3_g3,
      pname_hgb_1_g1,
      pname_hgb_2_g1,
      pname_hgb_3_g1,
      pname_hgb_1_g2,
      pname_hgb_2_g2,
      pname_hgb_3_g2,
      pname_hgb_1_g3,
      pname_hgb_2_g3,
      pname_hgb_3_g3,
      pname_plat_1_g1,
      pname_plat_2_g1,
      pname_plat_3_g1,
      pname_plat_1_g2,
      pname_plat_2_g2,
      pname_plat_3_g2,
      pname_plat_1_g3,
      pname_plat_2_g3,
      pname_plat_3_g3
    )
  ),
  plotwidth = 8,
  orientation = "landscape",
  opath = write_path(opath),
  file = tblid,
  title = title_footer$title,
  footers = title_footer$main_footer
)
