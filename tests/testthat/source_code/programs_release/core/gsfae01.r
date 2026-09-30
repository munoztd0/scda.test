###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsfae01.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create gsfae01: Subjects With Treatment
##                            -emergent Adverse Events by System Organ Class
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adae
## Output:                    gsfae01.rtf
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

tblid <- "GSFAE01"
popfl <- "SAFFL"

################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

# prepare data frame for frequency count -----------------------

# get subjects count of each treatment group
adsl_n <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) |>
  df_na() |>
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(
    TRT01A = factor(
      TRT01A,
      levels = c(
        "Xanomeline Low Dose",
        "Xanomeline High Dose",
        "Placebo"
      )
    )
  ) |>
  group_by(TRT01A) |>
  summarize(total = n()) |>
  ungroup()

# get subjects count of each AEBODSYS
adae_n <- haven::read_sas(envsetup::read_path(a_in, "adae.sas7bdat")) |>
  mutate(
    AEBODSYS = case_when(
      AEBODSYS == "" ~ "Uncoded",
      .default = AEBODSYS
    )
  ) |>
  df_na() |>
  filter(TRTEMFL == "Y", !!rlang::sym(popfl) == "Y") |>
  mutate(
    TRT01A = factor(
      TRT01A,
      levels = c(
        "Xanomeline Low Dose",
        "Xanomeline High Dose",
        "Placebo"
      )
    )
  ) |>
  mutate(
    AEBODSYS = forcats::fct_recode(
      AEBODSYS,
      "Neoplasms benign, malignant and unspecified" = "Neoplasms benign, malignant and unspecified (incl cysts and polyps)",
      "General disorders, administration site conditions" = "General disorders and administration site conditions"
    )
  ) |>
  group_by(TRT01A, AEBODSYS) |>
  summarise(n = n_distinct(USUBJID)) |>
  ungroup() |>
  # fill in missing values for AEBODSYS which don't have incidence subjects
  tidyr::complete(TRT01A, tidyr::nesting(AEBODSYS), fill = list(n = 0))

# merge data to derive frequency count
adae_freq <- adae_n |>
  left_join(adsl_n, by = "TRT01A") |>
  mutate(
    prop = tidytlg::roundSAS((n / total) * 100, digits = 1),
    propC = tidytlg::roundSAS((n / total) * 100, digits = 1, as_char = TRUE)
  )
# propC = paste0(tidytlg::roundSAS((n/total)*100, digits = 1, as_char = TRUE), "%"))

dt_freq <- function(trt) {
  df <- adae_freq |>
    # select compared active and placebo groups
    filter(TRT01A %in% c(trt, "Placebo")) |>
    # flag higher/lower numbers to present on right/left side each
    group_by(AEBODSYS) |>
    mutate(flag = ifelse(prop == max(prop), 1, 0)) |>
    # add total number to present in legend label
    rename(TRT01A_ = TRT01A) |>
    mutate(TRT01A = paste0(TRT01A_, " (N=", total, ")")) |>
    ungroup()

  # Restore levels
  df$TRT01A <- factor(
    df$TRT01A,
    levels = c(
      unique(df$TRT01A[df$TRT01A_ == "Xanomeline Low Dose"]),
      unique(df$TRT01A[df$TRT01A_ == "Xanomeline High Dose"]),
      unique(df$TRT01A[df$TRT01A_ == "Placebo"])
    )
  )

  # create sorting order by decreasing incidence of System Organ Class
  ae_order <- adae_freq |>
    filter(TRT01A == trt) |>
    arrange(desc(prop)) |>
    mutate(srt = row_number()) |>
    select(AEBODSYS, srt)

  df <- df |>
    left_join(ae_order, by = "AEBODSYS") |>
    select(-TRT01A_)
}

df_freq_1 <- dt_freq("Xanomeline High Dose")
df_freq_2 <- dt_freq("Xanomeline Low Dose")


# prepare data frame for risk difference and 95% CI -------------

adsl_pt <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) |>
  df_na() |>
  filter(SAFFL == "Y") |>
  select(USUBJID, TRT01A)

adae_pt <- haven::read_sas(envsetup::read_path(a_in, "adae.sas7bdat")) |>
  mutate(
    AEBODSYS = case_when(
      AEBODSYS == "" ~ "Uncoded",
      .default = AEBODSYS
    )
  ) |>
  df_na() |>
  filter(TRTEMFL == "Y", SAFFL == "Y") |>
  mutate(
    AEBODSYS = forcats::fct_recode(
      AEBODSYS,
      "Neoplasms benign, malignant and unspecified" = "Neoplasms benign, malignant and unspecified (incl cysts and polyps)",
      "General disorders, administration site conditions" = "General disorders and administration site conditions"
    )
  ) |>
  distinct(USUBJID, TRT01A, AEBODSYS) |>
  mutate(rsp = TRUE)

# create each row response for subjects and AEBODSYS
adae_rsp <- tidyr::expand_grid(
  distinct(adsl_pt, USUBJID, TRT01A),
  distinct(adae_pt, AEBODSYS)
) |>
  left_join(adae_pt, by = c("USUBJID", "TRT01A", "AEBODSYS")) |>
  mutate(rsp = ifelse(is.na(rsp), FALSE, rsp))


# get AEBODSYS unique terms (used below in for loop)
ae_term <- as.vector(unlist(distinct(adae_rsp, AEBODSYS)))


dt_diff <- function(trt) {
  # select compared active and placebo groups
  df <- adae_rsp |>
    filter(TRT01A %in% c(trt, "Placebo"))

  # keep only two levels for modeling and adjust level order in appropriate
  df$TRT01A <- droplevels(df$TRT01A)
  df$TRT01A <- factor(df$TRT01A, levels = c("Placebo", trt))

  # apply proportion difference method to
  # calculate risk difference and 95% CI of each AEBODSYS

  datalist <- list()

  for (i in ae_term) {
    ae <- subset(df, AEBODSYS == i)
    # here call down to prop_diff_wald, check in ?prop_diff_wald to select desired method
    res <- prop_diff_wald(rsp = ae$rsp, grp = ae$TRT01A)
    datalist[[i]] <- as.data.frame(res) |>
      mutate(
        ci = row_number(),
        ci = case_when(
          ci == 1 ~ "lower",
          ci == 2 ~ "upper"
        )
      ) |>
      tidyr::pivot_wider(
        id_cols = diff,
        names_from = ci,
        values_from = diff_ci
      ) |>
      mutate(
        diff = tidytlg::roundSAS(diff * 100, digits = 1),
        upper = tidytlg::roundSAS(upper * 100, digits = 1),
        lower = tidytlg::roundSAS(lower * 100, digits = 1)
      )
  }

  dt <- do.call(rbind, datalist)
  dt$AEBODSYS <- row.names(dt)

  # create sorting order by decreasing incidence of System Organ Class
  ae_order <- adae_freq |>
    filter(TRT01A == trt) |>
    arrange(desc(prop)) |>
    mutate(srt = row_number()) |>
    select(AEBODSYS, srt)

  dt <- dt |>
    left_join(ae_order, by = "AEBODSYS")

  return(dt)
}

df_diff_1 <- dt_diff("Xanomeline High Dose")
df_diff_2 <- dt_diff("Xanomeline Low Dose")


################################################################################
# Generate plot:
################################################################################

# define function for creating the body of graphic:
g_diff <- function(df_freq, df_diff) {
  # frequency plot ----------------------------------
  # Calculate padding for lower limit to allow space for labels
  lower_padding <- -4

  plot_freq <- ggplot(
    df_freq,
    aes(x = reorder(.data$AEBODSYS, -.data$srt), y = .data$prop)
  ) +
    geom_point(aes(shape = .data$TRT01A, color = .data$TRT01A)) +
    geom_text(
      data = df_freq |> filter(.data$flag == 0),
      aes(label = .data$propC),
      hjust = 1.3,
      size = 3
    ) +
    geom_text(
      data = df_freq |> filter(.data$flag == 1),
      aes(label = .data$propC),
      hjust = -0.3,
      size = 3
    ) +
    coord_flip() +
    # Use scales package for automatic breaks calculation with a lower limit padding
    scale_y_continuous(
      limits = c(lower_padding, NA), # Only set lower limit, let upper be automatic
      breaks = scales::extended_breaks(n = 5), # Automatically calculate appropriate breaks
      labels = scales::label_number(accuracy = 1) # Format labels
    ) +
    # assign colorblind friendly palette
    scale_color_manual(values = cbbPalette) +
    labs(x = "", y = "Frequency (%)") +
    theme_bw() +
    theme(
      # text = element_text(size = 9, color = "black"),
      axis.text.y = element_text(size = 7, color = "black"),
      axis.ticks = element_blank(),
      axis.title.x = element_text(face = "bold", size = 7, vjust = -2),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 7, face = "bold")
    )

  # risk difference plot ----------------------------
  plot_diff <- ggplot(
    df_diff,
    aes(x = reorder(.data$AEBODSYS, -.data$srt), y = .data$diff)
  ) +
    geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), width = .4) +
    geom_line(aes()) +
    geom_point(aes(), size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip() +
    # Use scales package for automatic breaks calculation
    scale_y_continuous(
      breaks = scales::extended_breaks(n = 5), # Automatically calculate appropriate breaks
      labels = scales::label_number(accuracy = 1), # Format labels
      expand = expansion(mult = c(0.1, 0.1)) # Add 10% padding on both sides
    ) +
    labs(x = " ", y = "Risk Difference With 95% CI") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_text(face = "bold", size = 7, vjust = -2),
      legend.position = "none",
    )

  # compose final object by putting plot, table, legend together --------
  final <- (plot_freq | plot_diff) +
    plot_layout(widths = c(6.5, 3.5))
}


# define parameters for plotting:
# assign colorblind friendly palette: black(Xan Low), orange(Xan High), dark blue(PBO)
cbbPalette <- c("#000000", "#E69F00", "#0072B2")


# call for Xan High vs PBO
pt_xanhigh <- g_diff(
  df_freq = df_freq_1,
  df_diff = df_diff_1
)

# call for Xan Low vs PBO
pt_xan <- g_diff(
  df_freq = df_freq_2,
  df_diff = df_diff_2
)


################################################################################
# Create png and output file:
################################################################################

# create png file and output figure
pname_1 <- paste0(tolower(tblid), "_1", ".png")

png(
  write_path(opath, pname_1),
  width = 40,
  height = 14,
  units = "cm",
  res = 600,
  type = "cairo"
)
print(write_path(opath, pname_1)) ### print png path and name in log
print(pt_xan)
dev.off()

pname_2 <- paste0(tolower(tblid), "_2", ".png")

png(
  write_path(opath, pname_2),
  width = 40,
  height = 14,
  units = "cm",
  res = 600,
  type = "cairo"
)
print(write_path(opath, pname_2)) ### print png path and name in log
print(pt_xanhigh)
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
