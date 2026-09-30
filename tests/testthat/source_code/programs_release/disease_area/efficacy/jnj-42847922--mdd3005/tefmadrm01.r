###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefmadrm01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefmadrm01: Subjects Who Achieved
##                            Remission at [Time Point] - [Observed Case][ - Double-blind Phase];
##                            Full Analysis Set Analysis Set (Study mdd)
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADMADRSI
## Output:                    TEFMADRM01.rtf
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
tblid <- "TEFMADRM01"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Quetiapine XR"

# Define population flag used.
popfl <- "FAS1CON"

# Define the strata variables to use below.
strata <- c("RRSIWRS", "REGION1", "AGEGR2")
# For unstratified analyses, please use:
# strata <- NULL

# Define response parameter to be used.
resppar <- "MADR0117"

# Define response label to be used
# (could be "remission" or something else too).
resplbl <- "remission"

# Define between one and two sided p-values.
pval_sided <- "2"
# one of:
# "2" - two sided
# "-1" - one sided, less
# "1" - one sided, greater

# Define test method.
pval_method <- "cmh"
# one of:
# "chisq" - Chi-square test
# "fisher" - Fisher's exact test
# "cmh" - CMH test

# Define odds ratio estimation method.
or_method <- "exact"
# one of:
# "exact": odds ratio by conditional logistic regression
# "cmh": odds ratio by Cochran-Mantel-Haenszel (CMH)

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived formats specifications.
formats <- list(
  n_prop = jjcsformat_count_fraction,
  diff_est_ci = jjcsformat_xx("xx.x (xx.x, xx.x)"),
  or_rr_est_ci = jjcsformat_xx("xx.xx (xx.xx, xx.xx)"),
  pval = jjcsformat_pval_fct(alpha)
)

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

## ADMADRSI ----

admadrsi <- haven::read_sas(read_path(a_in, "admadrsi.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  filter(PARAMCD == resppar) |>
  mutate(AVISIT = factor(AVISIT, levels = get_visit_levels(AVISIT, AVISITN))) |>
  select(USUBJID, AVISIT, AVALC) |>
  mutate(is_resp = (AVALC == "Y"))

## Analysis ----

ana <- admadrsi |>
  inner_join(adsl, by = "USUBJID")

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
  ## For each visit ----
  split_rows_by("AVISIT", section_div = "") |>
  ## Counts ----
  summarize_row_groups(
    label_fstr = paste0("Subjects evaluable for ", resplbl, " at %s~[super a]"),
    format = "xx"
  ) |>
  ## Proportions ----
  estimate_proportion(
    "is_resp",
    conf_level = conflvl,
    .stats = "n_prop",
    .formats = c(n_prop = formats$n_prop),
    .labels = c(n_prop = paste0("Subjects with ", resplbl))
  ) |>
  ## Proportion difference ----
  analyze(
    vars = "is_resp",
    afun = a_proportion_diff_j,
    show_labels = "hidden",
    na_str = default_na_str(),
    table_names = "prop_diff",
    extra_args = list(
      conf_level = conflvl,
      method = "wald",
      .stats = "diff_est_ci",
      .labels = c(
        diff_est_ci = paste0(
          "% Difference (",
          tern::f_conf_level(conflvl),
          ")~[super b]"
        )
      ),
      .formats = c(diff_est_ci = formats$diff_est_ci),
      ref_path = ref_path
    )
  ) |>
  ## Odds ratio ----
  analyze(
    vars = "is_resp",
    afun = a_odds_ratio_j,
    table_names = "odds_ratio",
    na_str = default_na_str(),
    show_labels = "hidden",
    extra_args = list(
      variables = list(arm = trtvar, strata = strata),
      conf_level = conflvl,
      method = or_method,
      na_if_no_events = TRUE,
      .stats = "or_ci",
      .labels = c(
        or_ci = paste0(
          "Odds ratio (",
          tern::f_conf_level(conflvl),
          ")~[super c,d]"
        )
      ),
      .formats = c(or_ci = formats$or_rr_est_ci),
      .indent_mods = c(or_ci = 0L),
      ref_path = ref_path
    )
  ) |>
  ## Relative risk ----
  analyze(
    "is_resp",
    afun = a_relative_risk,
    table_names = "rel_risk",
    show_labels = "hidden",
    na_str = default_na_str(),
    extra_args = list(
      variables = list(strata = strata),
      conf_level = conflvl,
      .stats = "rel_risk_ci",
      .labels = c(
        rel_risk_ci = paste0(
          "Relative risk of ",
          resplbl,
          " (",
          tern::f_conf_level(conflvl),
          ")~[super c]"
        )
      ),
      .formats = c(rel_risk_ci = formats$or_rr_est_ci),
      .indent_mods = c(rel_risk_ci = 0L),
      ref_path = ref_path
    )
  ) |>
  ## P-value ----
  analyze(
    "is_resp",
    afun = a_test_proportion_diff,
    table_names = "pval",
    show_labels = "hidden",
    na_str = default_na_str(),
    extra_args = list(
      variables = list(strata = strata),
      method = pval_method,
      alternative = switch(
        pval_sided,
        "2" = "two.sided",
        "-1" = "less",
        "1" = "greater"
      ),
      .labels = c(pval = paste0(abs(as.numeric(pval_sided)), "-sided p-value")),
      .formats = c(pval = formats$pval),
      .indent_mods = c(pval = 0L),
      ref_path = ref_path
    )
  )

# Output ----

result <- build_table(lyt, df = ana, alt_counts_df = adsl)


# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
