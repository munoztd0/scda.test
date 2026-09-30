###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsfae04.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create gsfae04: Subjects With Treatment-
##                            emergent Adverse Events With Frequency ≥[xx]% in
##                            [Any Treatment Group] by OCMQ (Narrow)
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:
## Output:
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

tblid <- "GSFAE04"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
popfl <- "SAFFL"
trtvar <- "TRT01A"

################################################################################
# Process data:
################################################################################

# prepare data frame for frequency count -----------------------

adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) |>
  df_na() |>
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c(
        "Xanomeline Low Dose",
        "Xanomeline High Dose",
        "Placebo"
      )
    )
  )

adaeocmq <- haven::read_sas(envsetup::read_path(a_in, "adaeocmq.sas7bdat")) |>
  df_na() |>
  mutate(
    TRT01A = factor(
      TRT01A,
      levels = c(
        "Xanomeline Low Dose",
        "Xanomeline High Dose",
        "Placebo"
      )
    )
  )

# get subjects count of each treatment group
adsl_n <- adsl |>
  filter(!!sym(popfl) == "Y") |>
  group_by(!!sym(trtvar)) |>
  summarize(total = n()) |>
  ungroup()

# get subjects count of each OCMQNAM
adae_n <- adaeocmq |>
  filter(TRTEMFL == "Y", !!sym(popfl) == "Y", toupper(OCMQCLSS) == "NARROW", !is.na(OCMQNAM)) |>
  group_by(!!sym(trtvar), OCMQNAM) |>
  summarise(n = n_distinct(USUBJID)) |>
  ungroup() |>
  # fill in missing values for OCMQNAM which don't have incidence subjects
  complete(!!sym(trtvar), nesting(OCMQNAM), fill = list(n = 0))

# merge data to derive frequency count
adae_freq <- adae_n |>
  left_join(adsl_n, by = trtvar) |>
  mutate(
    prop = tidytlg::roundSAS((n / total) * 100, digits = 1),
    propC = tidytlg::roundSAS((n / total) * 100, digits = 1, as_char = TRUE)
  )

dt_freq <- function(trt) {
  df <- adae_freq |>
    # select compared active and placebo groups
    filter(!!sym(trtvar) %in% c(trt, "Placebo")) |>
    # flag higher/lower numbers to present on right/left side each
    group_by(OCMQNAM) |>
    mutate(
      flag = ifelse(prop == max(prop), 1, 0),
      TRT_ = !!sym(trtvar)
    ) |>
    # add total number to present in legend label
    mutate(
      !!sym(trtvar) := paste0(!!sym(trtvar), " (N=", total, ")")
    ) |>
    ungroup()

  # Restore levels
  df[[trtvar]] <- factor(
    df[[trtvar]],
    levels = c(
      unique(df[[trtvar]][df$TRT_ == "Xanomeline Low Dose"]),
      unique(df[[trtvar]][df$TRT_ == "Xanomeline High Dose"]),
      unique(df[[trtvar]][df$TRT_ == "Placebo"])
    )
  )

  df <- df |>
    mutate(trtmntvr = !!sym(trtvar))

  # create sorting order by decreasing incidence of OCMQNAM
  ae_order <- adae_freq |>
    filter(!!sym(trtvar) == trt) |>
    arrange(desc(prop)) |>
    mutate(srt = row_number()) |>
    select(OCMQNAM, srt)

  df <- df |>
    left_join(ae_order, by = "OCMQNAM") |>
    select(-TRT_)
}

df_freq_1 <- dt_freq("Xanomeline High Dose")
df_freq_2 <- dt_freq("Xanomeline Low Dose")


# prepare data frame for risk difference and 95% CI -------------

adsl_pt <- adsl |>
  filter(!!sym(popfl) == "Y") |>
  select(USUBJID, !!sym(trtvar))

adae_pt <- adaeocmq |>
  filter(TRTEMFL == "Y", !!sym(popfl) == "Y", !is.na(OCMQNAM)) |>
  distinct(USUBJID, !!sym(trtvar), OCMQNAM) |>
  mutate(rsp = TRUE)

# create each row response for subjects and OCMQNAM
adae_rsp <- expand_grid(
  distinct(adsl_pt, USUBJID, !!sym(trtvar)),
  distinct(adae_pt, OCMQNAM)
) |>
  left_join(adae_pt, by = c("USUBJID", trtvar, "OCMQNAM")) |>
  mutate(rsp = ifelse(is.na(rsp), FALSE, rsp))


# get OCMQNAM unique terms (used below in for loop)
ae_term <- as.vector(unlist(distinct(adae_rsp, OCMQNAM)))


dt_diff <- function(trt) {
  # select compared active and placebo groups
  df <- adae_rsp |>
    filter(!!sym(trtvar) %in% c(trt, "Placebo"))

  df <- df |>
    mutate(
      !!sym(trtvar) := droplevels(!!sym(trtvar)),
      !!sym(trtvar) := factor(!!sym(trtvar), levels = c("Placebo", trt)),
      trtmntvr = !!sym(trtvar)
    )

  # apply proportion difference method to
  # calculate risk difference and 95% CI of each OCMQNAM

  datalist <- list()

  for (i in ae_term) {
    ae <- subset(df, OCMQNAM == i)
    # here call down to prop_diff_wald, check in ?prop_diff_wald to select desired method
    res <- prop_diff_wald(rsp = ae$rsp, grp = ae$trtmntvr)
    datalist[[i]] <- as.data.frame(res) |>
      mutate(
        ci = row_number(),
        ci = case_when(
          ci == 1 ~ "lower",
          ci == 2 ~ "upper"
        )
      ) |>
      pivot_wider(
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
  dt$OCMQNAM <- row.names(dt)

  # create sorting order by decreasing incidence of FMQ
  ae_order <- adae_freq |>
    filter(!!sym(trtvar) == trt) |>
    arrange(desc(prop)) |>
    mutate(srt = row_number()) |>
    select(OCMQNAM, srt)

  dt <- dt |>
    left_join(ae_order, by = "OCMQNAM")

  return(dt)
}

df_diff_1 <- dt_diff("Xanomeline High Dose")
df_diff_2 <- dt_diff("Xanomeline Low Dose")


################################################################################
# Filter table to only keep those that meet x% criteria on any treatment column
################################################################################

df_freq1 <- df_freq_1 |>
  group_by(OCMQNAM) |>
  filter(all(prop >= 5)) |>
  ungroup()
df1_ocmqnam <- df_freq1 |>
  select(OCMQNAM) |>
  distinct()
df_diff1 <- df_diff_1 |> inner_join(df1_ocmqnam, by = c("OCMQNAM"))

df_freq2 <- df_freq_2 |>
  group_by(OCMQNAM) |>
  filter(all(prop >= 5)) |>
  ungroup()
df2_ocmqnam <- df_freq2 |>
  select(OCMQNAM) |>
  distinct()
df_diff2 <- df_diff_2 |> inner_join(df2_ocmqnam, by = c("OCMQNAM"))


################################################################################
# Generate plot:
################################################################################

# define function for creating the body of graphic:
g_diff <- function(df_freq, df_diff, brk_frq, lmt_frq, brk_diff, lmt_diff) {
  # frequency plot ----------------------------------
  plot_freq <- ggplot(
    df_freq,
    aes(x = reorder(.data$OCMQNAM, -.data$srt), y = .data$prop)
  ) +
    geom_point(aes(shape = .data$trtmntvr, color = .data$trtmntvr)) +
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
    scale_y_continuous(breaks = brk_frq, limits = lmt_frq) +

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
    aes(x = reorder(.data$OCMQNAM, -.data$srt), y = .data$diff)
  ) +
    geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), width = .4) +
    geom_line(aes()) +
    geom_point(aes(), size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip() +
    scale_y_continuous(breaks = brk_diff, limits = lmt_diff) +
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


# Create same Y axis values for both graph
min_val_frq_1 <- min(df_freq1$prop)
low_frq_1 <- floor(min_val_frq_1 / 10) * 10
max_val_frq_1 <- max(df_freq1$prop)
high_frq_1 <- ceiling(max_val_frq_1 / 10) * 10

min_val_diff_1 <- min(df_diff1$lower)
low_diff_1 <- floor(min_val_diff_1 / 10) * 10
max_val_diff_1 <- max(df_diff1$upper)
high_diff_1 <- ceiling(max_val_diff_1 / 10) * 10

###################################################################
min_val_frq_2 <- min(df_freq2$prop)
low_frq_2 <- floor(min_val_frq_2 / 10) * 10
max_val_frq_2 <- max(df_freq2$prop)
high_frq_2 <- ceiling(max_val_frq_2 / 10) * 10

min_val_diff_2 <- min(df_diff2$lower)
low_diff_2 <- floor(min_val_diff_2 / 10) * 10
max_val_diff_2 <- max(df_diff2$upper)
high_diff_2 <- ceiling(max_val_diff_2 / 10) * 10


low_frq <- min(low_frq_1, low_frq_2)
high_frq <- max(high_frq_1, high_frq_2)

low_diff <- min(low_diff_1, low_diff_2)
high_diff <- max(high_diff_1, high_diff_2)

# Generate breaks and limits for Xan High vs PBO

# y_breaks_freq <- sort(unique(c(seq(low_frq, high_frq, by = 20), high_frq)))
by_freq <- 20
y_breaks_freq <- seq(low_frq, ceiling(high_frq / by_freq) * by_freq, by = by_freq)
y_limits_freq <- c(low_frq, ceiling(high_frq / by_freq) * by_freq)

by_diff <- 10

y_breaks_diff <- seq(low_diff, ceiling(high_diff / by_diff) * by_diff, by = by_diff)
y_limits_diff <- c(low_diff, ceiling(high_diff / by_diff) * by_diff)

# call for Xan High vs PBO
pt_xanhigh <- g_diff(
  df_freq = df_freq1,
  brk_frq = y_breaks_freq,
  lmt_frq = y_limits_freq,
  df_diff = df_diff1,
  brk_diff = y_breaks_diff,
  lmt_diff = y_limits_diff
)

# Generate breaks and limits for Xan Low vs PBO
by_freq <- 20
# y_breaks_freq <- seq(low_frq, high_frq, by = 20)
y_breaks_freq <- seq(low_frq, ceiling(high_frq / by_freq) * by_freq, by = by_freq)
y_limits_freq <- c(low_frq, ceiling(high_frq / by_freq) * by_freq)


y_breaks_diff <- seq(low_diff, ceiling(high_diff / by_diff) * by_diff, by = by_diff)
y_limits_diff <- c(low_diff, ceiling(high_diff / by_diff) * by_diff)

# call for Xan Low vs PBO
pt_xan <- g_diff(
  df_freq = df_freq2,
  brk_frq = y_breaks_freq,
  lmt_frq = y_limits_freq,
  df_diff = df_diff2,
  brk_diff = y_breaks_diff,
  lmt_diff = y_limits_diff
)


################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- tab_titles

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
  res = 400,
  type = "cairo"
)
print(write_path(opath, pname_1)) ### print png path and name in log
print(pt_xan)
dev.off()

pname_2 <- paste0(tolower(tblid), "_2", ".png")

png(
  write_path(opath, pname_2),
  width = 22,
  height = 14,
  units = "cm",
  res = 400,
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
