###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefpasi05.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefpasi05: PASI Score Analysis
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     adsl.sas7bdat, adparspi.sas7bdat
## Output:                    tefpasi05.rtf
## Remarks:
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
library(junco)
library(haven)

################################################################################
# Define script level parameters:
################################################################################

tblid <- "TEFPASI05"
fileid <- write_path(opath, tblid)

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

popfl <- "FASFL"
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

# Subset for vists of interest.
visits <- c(
  "Week 0",
  "Week 1",
  "Week 2",
  "Week 4",
  "Week 8",
  "Week 12",
  "Week 16"
)

# Define visit at which to analyze the results, this must be included in visits above.
anavisit <- "Week 16"

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Define covariates and interactions.
covariates <- c("TRT01P", "AVISIT", "BASE", "STRATWTG")
interactions <- c("TRT01P:AVISIT", "BASE:AVISIT", "STRATWTG:AVISIT")

# Define formats specifications.
formats <- list(
  mean_sd = jjcsformat_xx("xx.xx (xx.xxx)"),
  median = jjcsformat_xx("xx.xx"),
  range = jjcsformat_xx("xx.x, xx.x"),
  iqr = jjcsformat_xx("xx.xx, xx.xx"),
  lsmean_estci = jjcsformat_xx("xx.x (xx.xx, xx.xx)"),
  lsmean_diffci = jjcsformat_xx("xx.x (xx.xx, xx.xx)"),
  pval = jjcsformat_pval_fct(alpha)
)

################################################################################
# Process data:
################################################################################

# Read in SAS dataset and convert to R dataframe.
adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(
    !!popfl := factor(.data[[popfl]]),
    !!trtvar := factor(.data[[trtvar]], levels = c(ctrlab, trtlab))
  ) |>
  create_colspan_var(
    non_active_grp = ctrlab,
    non_active_grp_span_lbl = " ",
    active_grp_span_lbl = "Active Study Agent",
    colspan_var = "colspan_trt",
    trt_var = trtvar
  ) |>
  select(
    USUBJID,
    !!rlang::sym(popfl),
    !!rlang::sym(trtvar),
    "STRATWTG",
    colspan_trt
  )

# Read in SAS dataset and convert to R dataframe.
adpasi <- haven::read_sas(read_path(a_in, "adpasii.sas7bdat")) |>
  filter(
    !!rlang::sym(popfl) == "Y" & PARAMCD == "PATOTP" & AVISIT %in% visits
  ) |>
  mutate(
    AVISIT = forcats::fct_reorder(factor(AVISIT), AVISITN),
    AVISITLBL = case_match(
      AVISIT,
      visits[1] ~ "Baseline",
      .default = paste0("Change from baseline at ", AVISIT)
    ),
    APCHG = APCHG * -1,
    AVAL = ifelse(AVISIT == "Week 0", BASE, APCHG)
  ) |>
  select(USUBJID, AVISIT, AVISITLBL, AVAL, BASE, APCHG)

adpasi <- inner_join(x = adsl, y = adpasi, by = "USUBJID")

################################################################################
# Define layout and build table:
################################################################################

# Map each treatment group to appropriate columns spanning header.
colspan_trt_map <- create_colspan_map(
  df = adsl,
  non_active_grp = ctrlab,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

lyt <- basic_table(
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
    "AVISIT",
    split_fun = keep_split_levels(only = c(visits[1], anavisit)),
    labels_var = "AVISITLBL"
  ) |>
  analyze_vars(
    vars = "AVAL",
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
    vars = "APCHG",
    afun = a_summarize_mmrm,
    na_str = default_na_str(),
    show_labels = "hidden",
    extra_args = list(
      variables = list(
        id = "USUBJID",
        arm = trtvar,
        visit = "AVISIT",
        covariates = c(covariates, interactions)
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

result <- build_table(lyt, adpasi, alt_counts_df = adsl)
result

################################################################################
# Post-Processing:
################################################################################

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

################################################################################
# Convert to tbl file and output table:
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
