###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefdlqi03.R
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefdlqi03: Change From Baseline in DLQI
##                            Component Scores at Week 8; Full Analysis Set (Study psoriasis)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADDLQII
## Output:                    TEFDLQI03.rtf
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

# Parameters ----

# Define output ID and file location.
tblid <- "TEFDLQI03"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)

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

# Define the covariates to use in the ANCOVA model.
# Note: The main effects and interaction of AVISIT and trtvar are always
# included already.
covariates <- c("BASE", "STRATWTG")

# Single timepoint only - only occurs in title not in body!
visits <- c("Week 8")

# Define response parameters to be used and analyzed separately.
# specify in the proper order here or later in levels of factor

resppar <- c(
  "DLSYFP",
  "DLDLYP",
  "DLLEIP",
  "DLWRKP",
  "DLPERP",
  "DLTRTP"
)

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived formats specifications.
formats <- list(
  mean_sd = jjcsformat_xx("xx.x (xx.xx)"),
  median = jjcsformat_xx("xx.x"),
  range = jjcsformat_xx("xx., xx."),
  iqr = jjcsformat_xx("xx.x, xx.x"),
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

## ADDLQII ----

addlqii <- haven::read_sas(read_path(a_in, "addlqii.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(PARAMCD %in% resppar, AVISIT %in% visits, !is.na(BASE)) |>
  select(USUBJID, AVISIT, PARAMCD, PARAM, AVAL, CHG, any_of(covariates)) |>
  mutate(
    USUBJID = factor(USUBJID),
    AVISIT = factor(AVISIT, levels = visits),
    PARAMCD = factor(PARAMCD, levels = resppar)
  )

## Analysis ----

ana <- addlqii |>
  inner_join(adsl, by = "USUBJID")

# Layout ----
ref_path <- c("colspan_trt", " ", trtvar, ctrlab)

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
    labels_var = "PARAM",
    page_by = FALSE,
    section_div = " "
  ) |>
  analyze_vars(
    "CHG",
    table_names = "descriptives",
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
    table_names = "model",
    afun = a_ancova,
    na_str = default_na_str(),
    show_labels = "hidden",
    var_labels = "bla",
    extra_args = list(
      variables = list(
        arm = trtvar,
        covariates = covariates
      ),
      conf_level = conflvl,
      weights_emmeans = "equal",
      ref_path = ref_path,
      .stats = c("lsmean_diff_with_ci", "pval"),
      .indent_mods = c(lsmean_diff_with_ci = 1, pval = 2),
      .labels = c(
        lsmean_diff_with_ci = paste0(
          "LS mean difference (",
          tern::f_conf_level(conflvl),
          ")"
        )
      ),
      .formats = c(
        lsmean_diff_with_ci = formats$lsmean_diffci,
        pval = formats$pval
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
