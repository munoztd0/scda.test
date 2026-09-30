###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefprom01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefprom01:
##                            TEFPSSD04: [Primary/Secondary Endpoint Analysis] ([Primary Estimand], [Composite Strategy]):
##                            Change From Baseline in PROMIS-29 Domain T-scores and Pain Intensity
##                            at [Time Point]; Full Analysis Set (Study psoriasis)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADPROMI
## Output:                    TEFPROM01.rtf
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
tblid <- "TEFPROM01"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "PLACEBO"

# Define the order of treatment group labels in the treatment variable.
trtlab <- c(
  "JNJ-77242113 25 MG QD",
  "JNJ-77242113 50 MG QD",
  "JNJ-77242113 25 MG BID",
  "JNJ-77242113 100 MG QD",
  "JNJ-77242113 100 MG BID"
)

# Define population flags used.
popfl <- "FASFL"

# Define the covariates to use in the MMRM models.
# Note: The main effects and interaction of AVISIT and trtvar are always
# included already.
covariates <- c("BASE", "STRATWTG")

# Additional interaction terms based on AVISIT, trtvar and covariates
# to be included in the MMRM models.
interactions <- c("BASE:AVISIT", "STRATWTG:AVISIT")

# Define visits in the right order. The first visit
# is used as the baseline visit.
visits <- c("Week 0", "Week 8", "Week 16")

# Define visit at which to analyze the results, this must be included in visits above.
anavisit <- "Week 16"

# Define response parameters to be used and analyzed separately.
resppar <- c(
  "PRPHYSP",
  "PRANXP",
  "PRDEPRSP",
  "PRFATIGP",
  "PRSLEEPP",
  "PRSOCIAP",
  "PRPNINFP",
  "PRPNINTP"
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
  lsmean_estci = jjcsformat_xx("xx.x (xx.xx, xx.xx)"),
  lsmean_diffci = jjcsformat_xx("xx.x (xx.xx, xx.xx)"),
  pval = jjcsformat_pval_fct(alpha)
)

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
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

## ADPROMI ----

adpromi <- haven::read_sas(read_path(a_in, "adpromi.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(PARAMCD %in% resppar, AVISIT %in% visits) |>
  select(USUBJID, AVISIT, PARAMCD, PARAM, AVAL, CHG, any_of(covariates)) |>
  mutate(
    USUBJID = factor(USUBJID),
    AVISIT = factor(AVISIT, levels = visits)
  )

## Analysis ----

ana <- adpromi |>
  inner_join(adsl, by = "USUBJID") |>
  mutate(
    AVISITLBL = case_match(
      AVISIT,
      visits[1] ~ "Baseline",
      .default = paste0("Change from baseline at ", AVISIT)
    ),
    # Important change here: For the descriptive table part, switch the variable
    # analyzed based on the visit, baseline for baseline visit and afterwards
    # change from baseline.
    AVAL = ifelse(AVISIT == visits[1], AVAL, CHG)
  )

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
  split_rows_by("PARAM", page_by = FALSE, section_div = " ") |>
  split_rows_by(
    "AVISIT",
    split_fun = keep_split_levels(only = c(visits[1], anavisit)),
    labels_var = "AVISITLBL"
  ) |>
  analyze_vars(
    "AVAL",
    show_labels = "hidden",
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
    control = control_analyze_vars(
      quantiles = c(0.25, 0.75),
      quantile_type = 2
    )
  ) |>
  analyze(
    "CHG",
    afun = a_summarize_mmrm,
    na_str = default_na_str(),
    show_labels = "hidden",
    extra_args = list(
      variables = list(
        arm = trtvar,
        covariates = covariates,
        id = "USUBJID",
        visit = "AVISIT"
      ),
      conf_level = conflvl,
      cor_struct = "unstructured",
      ref_levels = setNames(
        list(visits[1], ctrlab),
        c("AVISIT", trtvar)
      ),
      weights_emmeans = "equal",
      method = "Kenward-Roger",
      vcov = "Kenward-Roger-Linear",
      alternative = "two.sided",
      .stats = c("adj_mean_est_ci", "diff_mean_est_ci", "p_value"),
      .indent_mods = c(
        adj_mean_est_ci = 1,
        diff_mean_est_ci = 1,
        p_value = 2
      ),
      .labels = c(
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

result <- build_table(lyt, ana, alt_counts_df = adsl) |>
  prune_table(prune_func = keep_rows(keep_non_null_rows))

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = 
  result,
  file = fileid,
  fontspec = formatters::font_spec("Times", 8L, 1.2),
  orientation = "landscape"
)
