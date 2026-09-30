###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              gefresp01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Forest Plot of Subgroup Analysis - HR
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADBDC, ADRESP
## Output:                    GEFRESP01.png
## Remarks:                   Template R script version using rtables framework
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


library(patchwork)

################################################################################
# Define output ID:
################################################################################

tblid <- "GEFRESP01"

################################################################################
# Get titles and footnotes:
################################################################################

# title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(ITTFL == "Y") |>
  mutate(
    trt = factor(
      case_when(
        TRT01P == "Dummy A - Tec-Dara" ~ "Dummy A",
        TRT01P == "Dummy B - DPd" ~ "Dummy B",
        TRT01P == "Dummy B - DVd" ~ "Dummy B"
      ),
      levels = c("Dummy B", "Dummy A")
    ),
    AGEGR2 = factor(AGEGR2, levels = c("<65 years", ">=65 years")),
    SEX = factor(SEX, levels = c("M", "F"), labels = c("Male", "Female")),
    RACEGR2 = factor(RACEGR2, levels = c("White", "Asian", "Other")),
    RENGRP1 = factor(RENGRP1, levels = c("<=60 mL/min", ">60 mL/min")),
    ECOGGR1 = factor(
      case_when(
        ECOGBL == 0 ~ "0",
        ECOGBL == 1 ~ ">=1",
        ECOGBL == 2 ~ ">=1"
      ),
      levels = c("0", ">=1")
    ),
    STRAT01 = factor(STRAT01, levels = c("DPd", "DVd")),
    PRLNGRP = factor(
      case_when(
        PRLINES == 1 ~ "1",
        PRLINES == 2 ~ "2 or 3",
        PRLINES == 3 ~ "2 or 3"
      ),
      levels = c("1", "2 or 3")
    )
  ) |>
  select(USUBJID, trt, AGEGR2, SEX, RACEGR2, RENGRP1, ECOGGR1, STRAT01, PRLNGRP)


adbdc <- haven::read_sas(read_path(a_in, "adbdc.sas7bdat")) |>
  filter(
    ITTFL == "Y" &
      PARAMCD %in% c("BSLISS", "BSOFPLAS", "STDRISK", "HRISKABN", "BSLPLSCE")
  )

adbdc_1 <- adbdc |>
  filter(PARAMCD != "BSLPLSCE") |>
  pivot_wider(
    id_cols = USUBJID,
    names_from = PARAMCD,
    values_from = AVALC
  )

adbdc_2 <- adbdc |>
  filter(PARAMCD == "BSLPLSCE") |>
  pivot_wider(
    id_cols = USUBJID,
    names_from = PARAMCD,
    values_from = AVAL
  )

adbdc_3 <- adbdc_1 |>
  full_join(adbdc_2, by = "USUBJID") |>
  mutate(
    ISSGRP = factor(BSLISS, levels = c("I", "II", "III")),
    PLASGRP = factor(
      case_when(
        BSOFPLAS == "0" ~ "No",
        BSOFPLAS == ">=1" ~ "Yes"
      ),
      levels = c("No", "Yes")
    ),
    BSLRISK = factor(
      case_when(
        HRISKABN == "Y" ~ "High-risk",
        STDRISK == "Y" ~ "Standard-risk"
      ),
      levels = c("High-risk", "Standard-risk")
    ),
    PLSCGRP = factor(
      case_when(
        BSLPLSCE <= 30 ~ "<=30",
        BSLPLSCE < 60 ~ ">30 to <60",
        BSLPLSCE >= 60 ~ ">=60"
      ),
      levels = c("<=30", ">30 to <60", ">=60")
    )
  ) |>
  select(USUBJID, ISSGRP, PLASGRP, BSLRISK, PLSCGRP)

adeff <- haven::read_sas(read_path(a_in, "adeff.sas7bdat")) |>
  filter(ITTFL == "Y" & PARAMCD == "IRCBRESP") |>
  select(USUBJID, PARAMCD, AVALC)

adrsp <- adsl |>
  left_join(adbdc_3, by = "USUBJID") |>
  left_join(adeff, by = "USUBJID") |>
  mutate(
    is_rsp = AVALC %in%
      c("Complete Response (CR)", "Stringent Complete Response (sCR)")
  )

# prepare data frame for survival subgroup analysis ---------------

## define subgroup variables
subgrp <- c(
  "AGEGR2",
  "SEX",
  "RACEGR2",
  "RENGRP1",
  "ECOGGR1",
  "STRAT01",
  "PRLNGRP",
  "ISSGRP",
  "PLASGRP",
  "BSLRISK",
  "PLSCGRP"
)

## generate events, odd ratio and 95%CI from tern
df_grouped <- extract_rsp_subgroups(
  variables = list(
    rsp = "is_rsp",
    arm = "trt",
    subgroups = subgrp
    # strata = "STRATA2"
  ),
  data = adrsp,
  conf_level = 0.95
)

## format columns in output presentation
prop_df <- df_grouped$prop
or_df <- df_grouped$or

prop_df2 <- prop_df |>
  mutate(
    arm = case_when(
      arm == "Dummy A" ~ "col1",
      arm == "Dummy B" ~ "col2"
    )
  ) |>
  pivot_wider(
    names_from = arm,
    names_glue = "{arm}_{.value}",
    values_from = c(n, n_rsp, prop)
  )

df <- prop_df2 |>
  inner_join(or_df, by = c("subgroup", "var", "var_label", "row_type")) |>
  select(-arm) |>
  mutate(
    # rounding in SAS rules
    across(
      c(or, lcl, ucl),
      ~ tidytlg::roundSAS(.x, digits = 2, as_char = TRUE, na_char = "NE")
    ),
    across(c(col1_prop, col2_prop), ~ .x * 100),
    across(
      c(col1_prop, col2_prop),
      ~ tidytlg::roundSAS(.x, digits = 1, as_char = TRUE, na_char = "NE")
    ),
    # formatting as table presentation
    or_ci = paste0(or, " (", lcl, ", ", ucl, ")"),
    col1_rsp = paste0(col1_n_rsp, "/", col1_n, " (", col1_prop, "%)"),
    col2_rsp = paste0(col2_n_rsp, "/", col2_n, " (", col2_prop, "%)")
  ) |>
  select(subgroup, var, var_label, col1_rsp, col2_rsp, or_ci, or, lcl, ucl)

df$or <- as.numeric(df$or)
df$lcl <- as.numeric(df$lcl)
df$ucl <- as.numeric(df$ucl)

## insert rows for each subgroup variable with blank values
subgrp_rows <- data.frame(
  subgroup = "",
  var = subgrp,
  var_label = subgrp,
  col1_rsp = NA,
  col2_rsp = NA,
  or_ci = NA,
  or = NA,
  lcl = NA,
  ucl = NA,
  stringsAsFactors = FALSE
)

tbl_df <- rbind(df, subgrp_rows)
tbl_df$ord1 <- match(tbl_df$var, subgrp) # create order for each subgroup variable

## sort table in the defined order
tbl_df <- tbl_df |>
  mutate(
    ord1 = ifelse(is.na(ord1), 0, ord1),
    ord2 = ifelse(subgroup == "", 1, 2)
  ) |>
  arrange(ord1, ord2) |>
  mutate(
    subgroup = ifelse(subgroup == "All Patients", "ALL", subgroup),
    var_label = case_when(
      var == "ALL" ~ "All subjects",
      var == "AGEGR2" ~ "Age",
      var == "SEX" ~ "Sex",
      var == "RACEGR2" ~ "Race",
      var == "RENGRP1" ~ "Baseline renal function",
      var == "ECOGGR1" ~ "ECOG performance score",
      var == "STRAT01" ~ "Investigator's choice of DPd or DVd",
      var == "PRLNGRP" ~ "Number of lines of prior therapy",
      var == "ISSGRP" ~ "Baseline ISS",
      var == "PLASGRP" ~ "Baseline soft-tissue plasmacytomas",
      var == "BSLRISK" ~ "Cytogenetic risk groups",
      var == "PLSCGRP" ~ "Bone marrow % plasma cells"
    )
  ) |>
  # create ord variable for row ordering
  tibble::rownames_to_column(var = "ord") |>
  mutate(ord = as.numeric(ord))

################################################################################
# Generate plot:
################################################################################

# create sub-group label vector for y-axis
ytxt <- tbl_df |>
  mutate(
    row_lbl = ifelse(
      (subgroup == "" | subgroup == "ALL"),
      var_label,
      paste0("  ", subgroup)
    )
  ) |>
  pull(row_lbl)

## generate forest plot ---------------

plot <- ggplot(tbl_df, aes(x = or, y = ord)) +
  # plot 95%CI
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), linewidth = 0.5, height = 0.5) +
  # plot odd ratio points
  geom_point(shape = 16, size = 2) +
  # plot the vertical line of or = 1.0
  geom_vline(xintercept = 1, color = "black", cex = 0.25) +

  # define title and axis labels
  scale_x_continuous(
    limits = c(0.01, 100),
    breaks = c(0.01, 0.1, 1, 10, 100),
    trans = "log10",
    labels = scales::label_number()
  ) +
  scale_y_continuous(
    breaks = c(1:nrow(tbl_df)),
    labels = ytxt,
    trans = "reverse"
  ) +
  theme_bw() +
  labs(
    title = "Odd Ratio and 95% CI",
    x = "\u2190 Favor Dummy B             Favor Dummy A \u2192",
    y = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      size = 8,
      margin = margin(b = 0, unit = "cm")
    ),
    axis.title.x = element_text(hjust = 0.3, size = 8),
    axis.text.y = element_text(hjust = 0),
    panel.grid.minor = element_blank()
  )


## generate table plot --------------

tblpt <- ggplot(tbl_df, aes(y = ord)) +
  geom_text(aes(x = 1, label = col1_rsp), size = 2.5) +
  geom_text(aes(x = 2, label = col2_rsp), size = 2.5) +
  geom_text(aes(x = 3, label = or_ci), size = 2.5) +
  scale_x_continuous(
    position = "top",
    breaks = c(1, 2, 3), # Specify the positions where you want ticks
    labels = c("Dummy A \n n/N(%)", "Dummy B \n n/N(%)", "Odd Ratio (95% CI)")
  ) +
  scale_y_continuous(trans = "reverse") +
  labs(y = "", x = "") +
  theme_minimal() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 7),
    # To expand the right side margin
    plot.margin = unit(c(0, 1, 0, 0), "cm"),
    panel.grid = element_blank()
  ) +
  # To display the complete results
  coord_cartesian(clip = "off")


# compose final object by putting plot, table, legend together --------

final <- (plot | tblpt) +
  plot_layout(widths = c(5, 5))


################################################################################
# Create png and output file:
################################################################################

# create png file and output figure
pname <- paste0(tolower(tblid), ".png")

# png(write_path(opath, pname)
png(
  write_path(opath, pname),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
# print(write_path(opath, pname)) ### print png path and name in log
print(write_path(opath, pname))
print(final)
dev.off()

# if (length( title_footer$main_footer) == 0) {
#   title_footer$main_footer <- NULL
# }

title <- "Forest Plot of Subgroup Analyses on CR or Better Response Based on Independent Review Committee (IRC) Assessment; Intent-to-Treat Analysis Set (Study 64007957MMY3001)"
main_footer <- ""

tidytlg::gentlg(
  tlf = "g",
  plotnames = write_path(opath, pname),
  plotwidth = 8,
  orientation = "landscape",
  opath = write_path(opath),
  file = tblid,
  title = title, # title_footer$title
  footers = main_footer
) # title_footer$main_footer
