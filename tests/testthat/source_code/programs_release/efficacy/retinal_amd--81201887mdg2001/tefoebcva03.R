###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefoebcva03.r
## R Version:                 4.5.2
## junco Version:             0.1.6 (tbc)
## Short Description:         Program to create tefoebcva03:
##                            [Efficacy Variable]: Change From Baseline in
##                            [Treated vs. Untreated Eye] Over Time [by Subgroup]:
##                            MMRM [ - Double-blind Phase];
##                            [Analysis Set] Analysis Set (Study specialty ophthalmology)
## Author:                    C&SP Methodology
## Date:                      04 May 2026
## Input:                     adsl, adeff
## Output:                    tefoebcva03.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:                     1
## Modified By:
## Reporting effort:
## Date:
## Description:
################################################################################

# Environment ----

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))


library(dplyr)
library(rtables)
library(junco)
library(haven)
library(tern)

# Parameters ----

# Define output ID and file location.
tblid <- "TEFOEBCVA03"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Sham Procedure"

# Define population flags used.
popfl <- "FASFL"

# Define analysis set label.
poplbl <- "Full"

# Define subgroup variable used, alongside the label to be used and the order of the
# categories.
subgroup <- "SEX"
subgrlbl <- "Sex"
subgrorder <- c("F", "M")

# Define the covariates to use in the MMRM.
# Note: The main effects and interaction of AVISIT and trtvar are always
# included already.
covariates <- c("BASE", "STRAT1")

# Define analysis domain flags to be used (for ADEFF).
domfl <- quote(EST01RFL == "Y" & ANLGRPP == "Study Eye")

# Define response parameter to be used (for ADEFF).
resppar <- "SQGAAREA"

# Define post-baseline visits be used for analysis.
visits <- c("Month 3", "Month 6", "Month 9", "Month 12", "Month 15", "Month 18")

# Define between one and two sided p-values.
pval_sided <- "2"
# one of:
# "2" - two sided
# "-1" - one sided, less
# "1" - one sided, greater

# Define which multiplicity adjustment to use for multiple arm
# comparisons, in each visit of interest, for p-values and CIs.
mult_adjust <- "step-down-dunnett"
# one of:
# "none" - no adjustment
# "dunnett" - one-step Dunnett adjustment for p-values, simultaneous CIs.
# "step-down-dunnett" - step-down Dunnett adjustment for p-values.
# Note that in this case still the same simultaneous CIs are calculated
# with Dunnett. It is currently not possible to get step-down Dunnett
# adjusted CIs.

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived formats specifications.
formats <- list(
  adj_mean_se = jjcsformat_xx("xx.x (xx.xx)"),
  adj_mean_ci = jjcsformat_xx("(xx.xx, xx.xx)"),
  diff_mean_se = jjcsformat_xx("xx.x (xx.xx)"),
  diff_mean_ci = jjcsformat_xx("(xx.xx, xx.xx)"),
  diff_pval = jjcsformat_pval_fct(alpha)
)

adjusted_lbl <- ifelse(mult_adjust == "none", "", "Adjusted ")
set_default_na_str("NE")

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  mutate(!!trtvar := forcats::fct_relevel(!!rlang::sym(trtvar), ctrlab)) |>
  mutate(!!subgroup := factor(.data[[subgroup]], levels = subgrorder)) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(subgroup), any_of(covariates))

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

## ADEFF ----

adeff <- haven::read_sas(read_path(a_in, "adeff.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(!!domfl) |>
  filter(PARAMCD == resppar) |>
  filter(AVISIT %in% visits) |>
  select(USUBJID, AVISIT, AVAL, BASE, CHG) |>
  mutate(
    USUBJID = factor(USUBJID),
    AVISIT = factor(AVISIT, levels = visits)
  ) |>
  arrange(USUBJID, AVISIT)

## Analysis ----

ana <- adeff |>
  inner_join(adsl, by = "USUBJID")

# MMRM analysis ----

variables <- list(
  response = "CHG",
  covariates = covariates,
  id = "USUBJID",
  arm = trtvar,
  visit = "AVISIT"
)

ana_by_subgr <- split(ana, ana[[subgroup]], drop = TRUE)
model_by_subgr <- lapply(
  ana_by_subgr,
  FUN = fit_mmrm_j,
  vars = variables,
  cor_struct = "unstructured",
  conf_level = conflvl,
  method = "Kenward-Roger",
  vcov = "Kenward-Roger-Linear",
  weights_emmeans = "equal",
  mult_adj_emmeans = mult_adjust
)
results_by_subgr <- lapply(model_by_subgr, broom::tidy)
ana_results <- bind_rows(results_by_subgr, .id = subgroup)

results <- ana_results |>
  mutate(
    colspan_trt = factor(
      ifelse(!!rlang::sym(trtvar) == ctrlab, " ", "Active Study Agent"),
      levels = c("Active Study Agent", " ")
    )
  )

# Layout ----

lyt <- rtables::basic_table(
  top_level_section_div = " ",
  show_colcounts = FALSE
) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  split_rows_by(
    subgroup,
    page_by = TRUE,
    split_fun = drop_split_levels,
    split_label = subgrlbl,
    section_div = " "
  ) |>
  summarize_row_counts() |>
  split_rows_by(variables$visit, section_div = "") |>
  ## LS Means part ----
  analyze(
    trtvar,
    afun = a_lsmeans,
    show_labels = "hidden",
    na_str = default_na_str(),
    extra_args = list(
      .stats = c(
        "n",
        "adj_mean_se",
        "adj_mean_ci"
      ),
      .labels = c(
        n = "N",
        adj_mean_se = "LS Mean (SE)",
        adj_mean_ci = "95% CI"
      ),
      .formats = c(
        adj_mean_se = formats$adj_mean_se,
        adj_mean_ci = formats$adj_mean_ci
      ),
      .indent_mods = c(
        adj_mean_se = 1,
        adj_mean_ci = 1
      ),
      ref_path = ref_path
    )
  ) |>
  insert_blank_line() |>
  ## Difference in LS means part ----
  analyze(
    trtvar,
    table_names = "diff_ls_means",
    afun = a_lsmeans,
    show_labels = "visible",
    var_labels = "Difference in change from baseline between treatments (active-sham)",
    na_str = default_na_str(),
    extra_args = list(
      alternative = switch(
        pval_sided,
        `2` = "two.sided",
        `-1` = "less",
        `1` = "greater"
      ),
      .stats = c(
        "n",
        "diff_mean_se",
        "diff_mean_ci",
        "p_value"
      ),
      .labels = c(
        n = "N",
        diff_mean_se = "LS Mean (SE)",
        diff_mean_ci = paste0(adjusted_lbl, f_conf_level(conflvl)),
        p_value = paste0(
          adjusted_lbl,
          abs(as.numeric(pval_sided)),
          "-sided p-value~[super a,b]"
        )
      ),
      .formats = c(
        diff_mean_se = formats$diff_mean_se,
        diff_mean_ci = formats$diff_mean_ci,
        p_value = formats$pval
      ),
      .indent_mods = c(
        n = 1,
        diff_mean_se = 2,
        diff_mean_ci = 2,
        p_value = 2
      ),
      ref_path = ref_path
    )
  )

# Output ----

result <- build_table(
  lyt,
  df = results,
  alt_counts_df = adsl,
  round_type = "sas"
)

# Add blank lines after each subgroup page split label.
section_div_at_path(result, c(subgroup, "*", "@content", "*")) <- " "

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
