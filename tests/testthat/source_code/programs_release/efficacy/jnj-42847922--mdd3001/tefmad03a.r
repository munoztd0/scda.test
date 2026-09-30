###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefmad03a.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefmad03a: Change From Baseline[(DB)] to [Time Point]:
##                            MMRM [Observed Case] Analysis[ - Double-blind Phase];
##                            Full Analysis Set Analysis Set (Study mdd)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADMADRSI
## Output:                    TEFMAD03a.rtf
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
tblid <- "TEFMAD03a"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Placebo"

# Define population flags used.
popfl <- "FAS1FL"

# Define analysis set label.
poplbl <- "FAS1"

# Define the covariates to use in the MMRM.
# Note: The main effects and interaction of AVISIT and trtvar are always
# included already.
covariates <- c("BASE", "RRSIWRS", "COUNTRY", "AGEGR2")

# Define analysis domain flags to be used (for ADMADRS).
domfl <- quote(ANL03FL == "Y" & DTYPE != "ENDPOINT")

# Define response parameter to be used (for ADMADRS).
resppar <- "MADES1"

# Define time point to be used.
timept <- "Day 43"

# Define baseline time point label.
basetpt <- "Baseline (DB)"

# Define between one and two sided p-values.
pval_sided <- "2"
# one of:
# "2" - two sided
# "-1" - one sided, less
# "1" - one sided, greater

# Double blind labeling in the output or not.
doubleblind <- TRUE

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived formats specifications.
formats <- list(
  mean_sd = jjcsformat_xx("xx.x (xx.xx)"),
  median = jjcsformat_xx("xx.x"),
  range = jjcsformat_xx("xx, xx"),
  pval = jjcsformat_pval_fct(alpha),
  diff_mean_se = jjcsformat_xx("xx.x (xx.xx)"),
  diff_mean_ci = jjcsformat_xx("(xx.xx, xx.xx)")
)
dblabel <- if (doubleblind) " (DB)" else ""

set_default_na_str("NE")

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  mutate(!!trtvar := forcats::fct_relevel(!!rlang::sym(trtvar), ctrlab)) |>
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
  trt_var = trtvar,
  active_first = FALSE
)

ref_path <- c("colspan_trt", " ", trtvar, ctrlab)

## ADMADRSI ----

admadrsi <- haven::read_sas(read_path(a_in, "admadrsi.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(!!domfl) |>
  filter(PARAMCD == resppar) |>
  select(USUBJID, AVISIT, AVAL, BASE, CHG) |>
  mutate(
    USUBJID = factor(USUBJID),
    AVISIT = factor(AVISIT)
  )

## Analysis ----

# For descriptive statistics.
ana_desc <- admadrsi |>
  filter(AVISIT %in% c(basetpt, timept)) |>
  mutate(
    AVISIT = stringr::str_replace(AVISIT, timept, paste0(timept, dblabel))
  ) |>
  inner_join(adsl, by = "USUBJID")

# For MMRM.
ana_mmrm <- admadrsi |>
  inner_join(adsl, by = "USUBJID") |>
  # Drop baseline visit (because there is no change to model yet).
  filter(AVISIT != basetpt) |>
  mutate(AVISIT = droplevels(AVISIT))

# MMRM analysis ----

model_obj <- fit_mmrm_j(
  vars = list(
    response = "CHG",
    covariates = covariates,
    id = "USUBJID",
    arm = trtvar,
    visit = "AVISIT"
  ),
  data = ana_mmrm,
  cor_struct = "unstructured",
  conf_level = conflvl,
  method = "Kenward-Roger",
  vcov = "Kenward-Roger-Linear",
  weights_emmeans = "equal"
)
ana_mmrm_results <- broom::tidy(model_obj) |>
  filter(AVISIT == timept) |>
  mutate(
    colspan_trt = factor(
      ifelse(!!rlang::sym(trtvar) == ctrlab, " ", "Active Study Agent"),
      levels = c("Active Study Agent", " ")
    )
  )

# Layout ----

## Descriptive analysis ----

lyt_desc <- rtables::basic_table(
  top_level_section_div = " ",
  show_colcounts = FALSE
) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  analyze(
    "STUDYID",
    var_labels = paste0("Analysis Set: ", poplbl),
    afun = function(x, .alt_df) nrow(.alt_df),
    format = "xx"
  ) |>
  split_rows_by("AVISIT", section_div = "") |>
  analyze_values("AVAL", formats = formats) |>
  analyze_values(
    "CHG",
    var_labels = paste0("Change from baseline", dblabel),
    formats = formats,
    show_labels = "visible",
    nested = FALSE
  )

## MMRM results ----

lyt_mmrm <- rtables::basic_table(
  top_level_section_div = " ",
  show_colcounts = FALSE
) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  analyze(
    trtvar,
    afun = a_lsmeans,
    show_labels = "hidden",
    na_str = default_na_str(),
    extra_args = list(
      alternative = switch(
        pval_sided,
        `2` = "two.sided",
        `-1` = "less",
        `1` = "greater"
      ),
      .stats = c(
        "p_value",
        "diff_mean_se",
        "diff_mean_ci"
      ),
      .labels = c(
        p_value = paste0(
          abs(as.numeric(pval_sided)),
          "-sided p-value (minus ",
          ctrlab,
          ")~[super a,b]"
        ),
        diff_mean_se = "Diff. of LS Means (SE)",
        diff_mean_ci = f_conf_level(conflvl)
      ),
      .formats = c(
        p_value = formats$pval,
        diff_mean_se = formats$diff_mean_se,
        diff_mean_ci = formats$diff_mean_ci
      ),
      .indent_mods = c(
        p_value = 1,
        diff_mean_se = 1,
        diff_mean_ci = 1
      ),
      ref_path = ref_path
    )
  )

# Output ----

result_desc <- build_table(lyt_desc, df = ana_desc, alt_counts_df = adsl)
result_mmrm <- build_table(
  lyt_mmrm,
  df = ana_mmrm_results,
  alt_counts_df = adsl
)

# This workaround makes sure we can rbind the two table parts:
col_info(result_mmrm) <- col_info(result_desc)

result <- rbind(result_desc, result_mmrm)


# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
