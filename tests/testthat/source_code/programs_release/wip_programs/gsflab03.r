###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsflab03.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Hepatocellular Drug-induced Liver Injury Screening Plot
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:
## Output:
## Remarks:                   Include parameters: ALT or AST vs TBL
##
## Modification History:
## Rev #:
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
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

tblid <- "GSFLAB03"

################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

# reading data

sas_vs_rds_check("adlb", a_in)
adlb <- readRDS(read_path(a_in, "adlb.rds")) %>%
  filter(PARAMCD %in% c("ALT", "AST", "BILI"), APOBLFL == "Y", SAFFL == "Y") %>%
  filter(!grepl("Unscheduled", AVISIT)) %>%
  filter(!(AVISIT %in% c("Endpoint"))) %>%
  # add abbreviation to TRT01A
  mutate(
    TRT01A = forcats::fct_recode(
      TRT01A,
      "Xanomeline (High)" = "Xanomeline High Dose",
      "Xanomeline (Low)" = "Xanomeline Low Dose",
      "Placebo (PBO)" = "Placebo"
    )
  ) %>%
  filter(!is.na(ANRHI)) %>% # data to be cleaned
  mutate(
    ratio = AVAL / ANRHI,
    PARAMCD = case_when(
      PARAMCD %in% c("ALT", "AST") ~ "ALTAST",
      TRUE ~ as.character(PARAMCD)
    )
  )


# create maximum post-baseline data points

adlb_max <- adlb %>%
  group_by(USUBJID, TRT01A, PARAMCD) %>%
  summarise(ratio = max(ratio)) %>%
  tidyr::pivot_wider(
    values_from = ratio,
    names_from = PARAMCD
  )

# define Hy's Law cases

bili <- adlb %>%
  filter(ratio >= 2, PARAMCD == "BILI") %>%
  select(USUBJID, TRT01A, ADY, BILI = ratio)

altast <- adlb %>%
  filter(ratio >= 3, PARAMCD == "ALTAST") %>%
  select(USUBJID, TRT01A, startdy = ADY)

check <- inner_join(bili, altast) %>%
  filter(ADY - startdy <= 30)

hylaw <- bili %>%
  semi_join(check) %>%
  select(-ADY) %>%
  mutate(flag = "Y")

adlb_max <- adlb_max %>%
  left_join(hylaw, by = c("USUBJID", "TRT01A", "BILI"))


################################################################################
# Generate plot:
################################################################################

# define parameters for plotting:
# assign colorblind friendly palette: black(Xan Low), orange(Xan High), dark blue(PBO)
cbbPalette <- c("#000000", "#E69F00", "#0072B2")


# call for Hepatocellular Drug-induced Liver Injury Screening Plot
# ALT vs BILI
x_label <- "Maximum post-baseline ALT or AST (value/ULN)"
y_label <- "Maximum post-baseline billrubin (value/ULN)"


# eDISH plot ---------------------------------------------------
plot <- ggplot(
  adlb_max,
  aes(x = ALTAST, y = BILI, color = TRT01A, shape = TRT01A)
) +

  # Make each dot partially transparent, with 0.5 opacity on alpha
  geom_point(
    alpha = 0.5,
    size = 1.2,
    # or apply dodge point with random noise to avoid overlapping
    # position = position_dodge(width = 0.05),
    na.rm = T
  ) +

  # encircling points for Hy's Law cases
  geom_point(
    data = adlb_max %>% filter(flag == "Y"),
    # position = position_dodge(width = 0.05),
    pch = 21,
    size = 4,
    colour = "red"
  ) +

  # Use a hollow circle, triangle, and cross as choices for shape
  scale_shape_manual(values = c(1, 2, 3)) +
  scale_x_continuous(
    breaks = scales::extended_breaks(n = 10),
    labels = scales::label_number(accuracy = 1)
  ) +
  scale_y_continuous(
    breaks = scales::extended_breaks(n = 5),
    labels = scales::label_number(accuracy = 1)
  ) +
  geom_hline(yintercept = 2, color = "grey", linetype = 2) +
  geom_vline(xintercept = 3, color = "grey", linetype = 2) +

  # assign colorblind friendly palette
  scale_color_manual(values = cbbPalette) +

  # annotate text in quadrants
  annotate("text", x = 1, y = 12, size = 4, label = "Cholestasis") +
  annotate("text", x = 12, y = 12, size = 4, label = "Potential Hy's Law") +
  annotate("text", x = 12, y = 0, size = 4, label = "Temple's corollary") +
  labs(x = x_label, y = y_label) +
  theme_bw() +
  theme(
    text = element_text(size = 9, color = "black"),
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
