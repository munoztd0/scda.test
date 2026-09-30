###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefmadsrp04.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefmadsrp04: Time to [Event] - Log-Rank
##                            Test[ - Double-Blind Phase]; Full Analysis Set (Study mdd)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADTTE
## Output:                    TEFMADSRP04.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:                     1
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
## Description:
################################################################################

# Environment ----

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)


library(dplyr)
library(rtables)
library(junco)
library(haven)

# Parameters ----

# Define output ID and file location.
tblid <- "TEFMADSRP04"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
tab_titles$main_footer <- tab_titles$main_footer[c(1, 2, 3, 6)]

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT03P"

# Define control group label used in the treatment variable.
ctrlab <- "Placebo"

# Define population flag used.
popfl <- "FASMA1FL"

# Define the potential strata variables for the log-rank test.
strata <- NULL

# Define time-to-event parameter to be used.
ttepar <- "TTRELAPS"

# Define label to be used for this time-to-event endpoint.
ttelbl <- "relapse"

# Define time unit.
timeunit <- "day"

# Decide between one and two sided p-values.
pval_sided <- "2"
# one of:
# "2" - two sided
# "-1" - one sided, less (treatment has lower hazard than control)
# "1" - one sided, greater (treatment has higher hazard than control)
# If one sided, the direction should be detailed in the footnote.

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived formats.
formats <- list(
  km_event_rate = jjcsformat_xx("xx.x (xx.x, xx.x)"),
  lr_stat_df = jjcsformat_xx("xx.x (xx.)"),
  pval = jjcsformat_pval_fct(alpha)
)
options(tern_default_na_str = rep("NE", 3))

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl), all_of(strata))

adsl$colspan_trt <- factor(
  ifelse(adsl[[trtvar]] == ctrlab, " ", "Active Study Agent"),
  levels = c("Active Study Agent", " ")
)

colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = ctrlab,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

ref_path <- c("colspan_trt", " ", trtvar, ctrlab)

## ADTTE ----

adtte <- haven::read_sas(read_path(a_in, "adtte.sas7bdat")) |>
  select(STUDYID, USUBJID, PARAMCD, CNSR, AVAL, all_of(trtvar))

adtte <- adtte |>
  select(USUBJID, PARAMCD, CNSR, AVAL) |>
  filter(PARAMCD == ttepar) |>
  select(-PARAMCD) |>
  mutate(
    USUBJID = factor(USUBJID),
    EVENT = factor(
      ifelse(CNSR == 1, "Number censored (%)", "Number of events (%)"),
      levels = c("Number censored (%)", "Number of events (%)")
    ),
    is_event = (CNSR == 0)
  )

## Analysis ----

ana <- adtte |>
  right_join(adsl, by = "USUBJID")

# Layout ----

lyt <- rtables::basic_table(
  top_level_section_div = " ",
  show_colcounts = TRUE,
  colcount_format = "N=xx"
) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  ## Event Counts ----
  analyze(
    vars = "EVENT",
    afun = s_proportion_factor,
    show_labels = "visible",
    var_labels = paste0("Time to ", ttelbl, "~[super a,b] (", timeunit, "s)"),
    extra_args = list(
      show_total = "top",
      total_label = "Number of assessments",
      use_alt_counts = FALSE
    )
  ) |>
  insert_blank_line() |>
  ## Kaplan-Meier ----
  tern::surv_time(
    vars = "AVAL",
    is_event = "is_event",
    show_labels = "hidden",
    na_str = default_na_str(),
    control = tern::control_surv_time(
      conf_level = conflvl,
      conf_type = "log-log"
    ),
    .stats = c("quantiles_lower", "median_ci_3d", "quantiles_upper"),
    .labels = c(
      quantiles_lower = paste0("25% quartile (", tern::f_conf_level(conflvl), ")"),
      median_ci_3d = paste0("Median (", tern::f_conf_level(conflvl), ")"),
      quantiles_upper = paste0("75% quartile (", tern::f_conf_level(conflvl), ")")
    ),
    .formats = c(
      quantiles_lower = formats$km_event_rate,
      median_ci_3d = formats$km_event_rate,
      quantiles_upper = formats$km_event_rate
    ),
    .indent_mods = c(
      quantiles_lower = 1,
      median_ci_3d = 1,
      quantiles_upper = 1
    )
  ) |>
  insert_blank_line() |>
  ## Hazard Ratio ----
  analyze(
    vars = "AVAL",
    afun = a_coxph_pairwise,
    na_str = default_na_str(),
    var_labels = "Statistical test",
    show_labels = "visible",
    table_names = "coxph_stratified",
    indent_mod = 1,
    extra_args = list(
      is_event = "is_event",
      strata = strata,
      control = tern::control_coxph(
        conf_level = conflvl,
        ties = "breslow"
      ),
      alternative = switch(
        pval_sided,
        `2` = "two.sided",
        `-1` = "less",
        `1` = "greater"
      ),
      .stats = c("lr_stat_df", "pvalue"),
      .formats = c(pvalue = formats$pval, lr_stat_df = formats$lr_stat_df),
      .labels = c(
        pvalue = paste0(
          abs(as.numeric(pval_sided)),
          "-sided p-value~[super c]"
        ),
        lr_stat_df = paste0(
          "Chisq (Overall DF) - ",
          ifelse(length(strata), "stratified", "unstratified")
        )
      ),
      ref_path = ref_path
    )
  )

# Output ----

result <- build_table(lyt, ana, alt_counts_df = adsl)


# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
