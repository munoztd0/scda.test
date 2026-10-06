###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefpromis01a.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefpromis01a:
##                            Change From [Induction] Baseline in PROMIS-29 Domain
##                            T-scores and Pain Intensity Over Time
## Author:                    Technology Solutions
## Date:                      2026-09-30
## Input:                     adsl, adpromsi
## Output:                    tefpromis01a.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:
## Modified By:
## Reporting effort:
## Date:
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
tblid <- "TEFPROMIS01a"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TR03PG1"

# Define control group label used in the treatment variable.
ctrlab <- "Placebo SC q4w"

# Define the order of treatment group labels in the treatment variable.
trtlab <- c(
  "Guselkumab 100 mg SC q8w",
  "Guselkumab 200 mg SC q4w"
)

# Define population flags used.
popfl <- expr(FAS2FL == "Y" & RAND2FL == "Y" & .data[[trtvar]] != "" & !is.na(.data[[trtvar]]))

# Domain level flags used.
domfl <- expr(!is.na(AVAL) & BASETYPE == "Induction")

# Define the covariates to use in the MMRM models.
# Note: The main effects and interaction of AVISIT and trtvar are always
# included already.
covariates <- c("BASE", "RASTRAM1", "RASTRAM2")

# Additional interaction terms based on AVISIT, trtvar and covariates
# to be included in the MMRM models.
interactions <- c("BASE:AVISIT", "BASE:RASTRAM1", "BASE:RASTRAM2")

# Define visits in the right order incl. labels to use. The first visit
# is used as the baseline visit.
visits <- c(
  "INDUCTION BASELINE" ~ "Induction baseline",
  "MAINTENANCE BASELINE" ~ "Maintenance baseline",
  "WEEK M-28" ~ "Week M-28~[super c]",
  "WEEK M-44" ~ "Week M-44~[super c]"
)

# Define response parameters to be used and analyzed separately, together
# with the labels to use.
resppar <- c(
  "PDEPRESP" ~ "Depression (T-score)",
  "PANXP" ~ "Anxiety (T-score)",
  "PPHYSP" ~ "Physical function (T-score)",
  "PPNINFP" ~ "Pain interference (T-score)",
  "PFATIGP" ~ "Fatigue (T-score)",
  "PSLEEPP" ~ "Sleep disturbance (T-score)",
  "PSOCIALP" ~ "Ability to participate in social roles and activities (T-score)",
  "PPNINTP" ~ "Pain intensity (0-10)"
)

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived formats specifications.
formats <- list(
  mean_sd = jjcsformat_xx("xx.xx (xx.xxx)"),
  median = jjcsformat_xx("xx.xx"),
  range = jjcsformat_xx("xx.x, xx.x"),
  iqr = jjcsformat_xx("xx.xx, xx.xx"),
  lsmean_estci = jjcsformat_xx("xx.xx (xx.xx, xx.xx)"),
  lsmean_diffci = jjcsformat_xx("xx.xx (xx.xx, xx.xx)"),
  pval = jjcsformat_pval_fct(alpha)
)

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!popfl) |>
  mutate(!!trtvar := factor(.data[[trtvar]], levels = c(ctrlab, trtlab))) |>
  select(STUDYID, USUBJID, all_of(trtvar), any_of(covariates))

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

## ADPROMSI ----

adpromsi <- haven::read_sas(read_path(a_in, "adpromsi.sas7bdat")) |>
  filter(!!popfl) |>
  filter(!!domfl) |>
  filter(
    PARAMCD %in% sapply(resppar, junco::leftside),
    AVISIT %in% sapply(visits, junco::leftside)
  ) |>
  select(USUBJID, AVISIT, PARAMCD, AVAL, CHG, any_of(covariates)) |>
  mutate(
    USUBJID = factor(USUBJID),
    PARAMCD = factor(PARAMCD, levels = sapply(resppar, junco::leftside)),
    AVISIT = factor(AVISIT, levels = sapply(visits, junco::leftside)),
    PARAMLBL = recode_values(
      PARAMCD,
      !!!resppar
    ),
    AVISITLBL = recode_values(
      AVISIT,
      !!!visits
    )
  )

## Analysis ----

ana <- adpromsi |>
  inner_join(adsl, by = "USUBJID")

# Layout ----

# First determine baseline visit label and level.
baseline_visit_label <- rightside(visits[[1]])
baseline_visit_level <- leftside(visits[[1]])

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
    "PARAMCD",
    labels_var = "PARAMLBL",
    page_by = TRUE,
    section_div = " "
  ) |>
  split_rows_by(
    "AVISIT",
    labels_var = "AVISITLBL",
    section_div = " "
  ) |>
  analyze(
    "AVAL",
    afun = tern::a_summary,
    show_labels = "hidden",
    extra_args = list(
      .stats = c("n", "mean_sd", "median", "range", "quantiles"),
      .indent_mods = c(
        n = 0,
        mean_sd = 1,
        median = 1,
        range = 1,
        quantiles = 1
      ),
      .labels = c(
        n = "N",
        mean_sd = "Mean (SD)",
        median = "Median",
        range = "Min, max",
        quantiles = "Interquartile range"
      ),
      .formats = c(
        n = "xx",
        mean_sd = formats$mean_sd,
        median = formats$median,
        range = formats$range,
        quantiles = formats$iqr
      ),
      quantiles = c(0.25, 0.75),
      quantile_type = 2,
      conf_level = conflvl
    )
  ) |>
  analyze(
    "CHG",
    afun = a_summary_j_with_exclude,
    var_labels = paste("Change from", baseline_visit_label),
    extra_args = list(
      exclude_levels = list(AVISIT = baseline_visit_level),
      .stats = c("n", "mean_sd", "median", "range", "quantiles"),
      .indent_mods = c(
        n = 0,
        mean_sd = 1,
        median = 1,
        range = 1,
        quantiles = 1
      ),
      .labels = c(
        n = "N",
        mean_sd = "Mean (SD)",
        median = "Median",
        range = "Min, max",
        quantiles = "Interquartile range"
      ),
      .formats = c(
        n = "xx",
        mean_sd = formats$mean_sd,
        median = formats$median,
        range = formats$range,
        quantiles = formats$iqr
      ),
      quantiles = c(0.25, 0.75),
      quantile_type = 2,
      conf_level = conflvl
    )
  ) |>
  analyze(
    "CHG",
    afun = a_summarize_mmrm_with_exclude,
    na_str = default_na_str(),
    show_labels = "hidden",
    table_names = "chg_mmrm",
    extra_args = list(
      exclude_levels = list(AVISIT = baseline_visit_level),
      variables = list(
        arm = trtvar,
        covariates = c(covariates, interactions),
        id = "USUBJID",
        visit = "AVISIT"
      ),
      conf_level = conflvl,
      cor_struct = "unstructured",
      ref_levels = setNames(
        list(baseline_visit_level, ctrlab),
        c("AVISIT", trtvar)
      ),
      weights_emmeans = "equal",
      method = "Kenward-Roger",
      vcov = "Kenward-Roger-Linear",
      alternative = "two.sided",
      .stats = c("n", "adj_mean_est_ci", "diff_mean_est_ci", "p_value"),
      .indent_mods = c(
        n = 1,
        adj_mean_est_ci = 2,
        diff_mean_est_ci = 2,
        p_value = 3
      ),
      .labels = c(
        n = "N",
        adj_mean_est_ci = paste0("LS mean (", tern::f_conf_level(conflvl), ")"),
        diff_mean_est_ci = paste0(
          "LS mean difference (",
          tern::f_conf_level(conflvl),
          ")"
        ),
        p_value = "p-value"
      ),
      .formats = c(
        adj_mean_est_ci = formats$lsmean_estci,
        diff_mean_est_ci = formats$lsmean_diffci,
        p_value = formats$pval
      )
    )
  )

# Output ----

# Build table, and prune NULL rows.
result <- build_table(lyt, ana, alt_counts_df = adsl, round_type = "sas") |>
  prune_table(prune_func = keep_rows(keep_non_null_rows))

# We add section dividers only now (and not via insert_blank_line in layout)
# because we want to avoid multiple consecutive blank lines in the result.
section_div_at_path(result, c("PARAMCD", "*", "AVISIT", "*", "AVAL", "quantiles")) <- " "
section_div_at_path(result, c("PARAMCD", "*", "AVISIT", "*", "CHG", "quantiles")) <- " "

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = 
  result,
  file = fileid,
  fontspec = formatters::font_spec("Times", 8L, 1.2),
  orientation = "landscape"
)
