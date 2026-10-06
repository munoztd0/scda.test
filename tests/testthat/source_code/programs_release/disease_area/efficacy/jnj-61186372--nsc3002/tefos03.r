###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefos03.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefos03: Overall Survival by [Subgroup]–
##                            [Stratified/Unstratified] Analysis; Full Analysis Set (Study mmy, bc, lc)
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADTTEEF
## Output:                    TEFOS03.rtf
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
tblid <- "TEFOS03"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
# Warning: Title file should contain exactly one title record per Table ID
# Therefore need to make sure we only have one title record:
tab_titles$title <- tab_titles$title[1]
tab_titles$main_footer <- tab_titles$main_footer[c(2, 3, 4)]

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "CP"

# Define population flag used.
popfl <- "FASFL"

# Define subgroup variable used, alongside the label to be used.
subgroup <- "AGEGR1"
subgrlbl <- "Age Group"

# Define the strata variables for hazard ratio analysis.
strata <- c("STRATA1", "STRATA2", "STRATA3")

# Define time-to-event parameter to be used.
ttepar <- "RPFS"

# Define time unit.
timeunit <- "month"

# Define time points for event free rates reporting.
timepts <- c(6, 9, 12)

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived formats.
formats <- list(
  pval = jjcsformat_pval_fct(alpha),
  est_ci = jjcsformat_xx("xx.xx (xx.xx, xx.xx)"),
  range = jjcsformat_range_fct("xx.x"),
  km_event_rate = jjcsformat_xx("xx.x (xx.x, xx.x)")
)
set_default_na_str("NE")

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(
    !!trtvar := as.factor(.data[[trtvar]]),
    !!subgroup := as.factor(.data[[subgroup]])
  ) |>
  select(
    STUDYID,
    USUBJID,
    all_of(trtvar),
    all_of(popfl),
    all_of(strata),
    all_of(subgroup)
  )

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

## ADTTEEF ----

adtteef <- haven::read_sas(read_path(a_in, "adtteef.sas7bdat")) |>
  select(STUDYID, USUBJID, PARAMCD, CNSR, AVAL, all_of(trtvar))

adtteef <- adtteef |>
  select(USUBJID, PARAMCD, CNSR, AVAL) |>
  filter(PARAMCD == ttepar) |>
  select(-PARAMCD) |>
  mutate(
    USUBJID = factor(USUBJID),
    EVENT = factor(
      ifelse(CNSR == 1, "Censored", "Event"),
      levels = c("Event", "Censored")
    ),
    is_event = (EVENT == "Event")
  )

## Analysis ----

ana <- adtteef |>
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
  split_rows_by(
    subgroup,
    page_by = TRUE,
    split_fun = drop_split_levels
  ) |>
  summarize_row_counts(label_fstr = paste0(subgrlbl, ": %s")) |>
  insert_blank_line() |>
  ## Event Counts ----
  analyze(
    vars = "EVENT",
    afun = s_proportion_factor,
    show_labels = "hidden"
  ) |>
  insert_blank_line() |>
  ## Kaplan-Meier ----
  tern::surv_time(
    vars = "AVAL",
    is_event = "is_event",
    var_labels = paste0(
      "Kaplan-Meier estimate of time to event (",
      timeunit,
      "s)"
    ),
    show_labels = "visible",
    na_str = default_na_str(),
    control = tern::control_surv_time(
      conf_level = conflvl,
      conf_type = "log-log"
    ),
    .formats = c(
      quantiles_lower = formats$est_ci,
      median_ci_3d = formats$est_ci,
      quantiles_upper = formats$est_ci,
      range_with_cens_info = formats$range
    )
  ) |>
  insert_blank_line()
## Event-Free Rates ----
for (time_point in timepts) {
  lyt <- lyt |>
    analyze(
      vars = "AVAL",
      afun = a_event_free,
      show_labels = "hidden",
      na_str = default_na_str(),
      table_names = paste0("AVAL_event_free_", time_point),
      extra_args = list(
        time_unit = timeunit,
        is_event = "is_event",
        time_point = time_point,
        control = tern::control_surv_time(
          conf_level = conflvl,
          conf_type = "log-log"
        ),
        .formats = c(event_free_ci = formats$km_event_rate),
        .stats = "event_free_ci"
      )
    )
}
lyt <- lyt |>
  insert_blank_line() |>
  ## Hazard Ratio ----
  analyze(
    vars = "AVAL",
    afun = a_coxph_pairwise,
    show_labels = "hidden",
    na_str = default_na_str(),
    table_names = "coxph_stratified",
    extra_args = list(
      is_event = "is_event",
      strata = strata,
      control = tern::control_coxph(
        conf_level = conflvl,
        ties = "breslow"
      ),
      .stats = c("hr_ci_3d", "pvalue"),
      .formats = c(hr_ci_3d = formats$est_ci, pvalue = formats$pval),
      .labels = c(
        hr_ci_3d = paste0(
          "Hazard ratio (",
          tern::f_conf_level(conflvl),
          ")~[super a]"
        ),
        pvalue = "p-value~[super b]"
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
