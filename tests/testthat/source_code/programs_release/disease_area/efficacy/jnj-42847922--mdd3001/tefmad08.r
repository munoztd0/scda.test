###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefmad08.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefmad08: Change From Baseline[(DB)] to
##                            [Time Point]: Tipping Point Multiple Imputation Analysis With Delta
##                            Adjustments[ - Double-blind Phase]; Full Analysis Set Analysis Set (Study mdd)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADMADRSI, ADICE, ADDISP
## Output:                    TEFMAD08.rtf
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

# Test run or real production run? Please set to FALSE for use in the study analysis.
testrun <- TRUE

# Define output ID and file location.
tblid <- "TEFMAD08"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Placebo"

# Define treatment group label to look at here (only one treatment group allowed).
trtlab <- "Seltorexant 20 mg"

# Define population flags used.
popfl <- "FAS1FL"

# Define the covariates to use in the MMRM (imputation and analysis).
# Note: AVISIT and trtvar are always included already.
covariates <- c("BASE", "RRSIWRS", "COUNTRY", "AGEGR2")

# Interaction terms based on covariates, AVISIT and trtvar
# to be used in the MMRM (imputation and analysis).
interactions <- c(paste(trtvar, "* AVISIT"))

# Define analysis domain flags to be used (for ADMADRSI).
domfl <- quote(ANL03FL == "Y" & DTYPE != "ENDPOINT")

# Define response parameter to be used (for ADMADRSI).
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

# Settings for the imputation method.
library(rbmi)
set.seed(3469) # Needed because method_bayes generates a random seed.
imputation_method <- rbmi::method_bayes(
  control = rbmi::control_bayes(
    warmup = 200,
    thin = ifelse(testrun, 1, 5)
  ),
  n_samples = ifelse(testrun, 10, 500) # Number of imputed datasets to be used.
)

# Reference based imputation method to be used.
reference_method <- "MAR"

# References to be used during imputation.
reference_spec <- c(
  "Placebo" = "Placebo",
  "Seltorexant 20 mg" = "Placebo"
)

# Define significance threshold to use (important for p-value formatting).
# Note that here a non-zero alpha should be used, because this is the threshold
# used to find the tipping point.
alpha <- 0.05

# Increasing sequence of delta values which are added to the defined subset
# of the control arm in turn, until the p-value is larger than alpha.
deltas <- seq(from = 0, to = 20, by = ifelse(testrun, 20, 1))

# Definition of the patients in the treatment arm where delta adjustment will be applied:
# - Patients are in ADICE (i.e. had an intercurrent event, ICE)
# - Patients do not have any records in ADDISP fulfilling the conditions in `delta_addisp_exclude_flags`:
delta_addisp_exclude_flags <- list(
  # No additional filters on top of being contained in ADICE:
  "All Subjects Who Discontinued and/or Switched Study Drug and/or Underlying Antidepressant" = quote(
    FALSE
  ),
  # Exclude patients from delta adjustment with the following disposition information:
  # PARAM = 'Reason for Treatment Discontinuation' and APHASEN = 1 and FAS1FL = 'Y';
  # & (AVALC in ('LOST TO FOLLOW-UP','WITHDRAWAL BY SUBJECT','OTHER')  or CRIT9FL = 'Y');
  "All Subjects Who Discontinued and/or Switched Study Drug and/or Underlying Antidepressant for Study Agent Reasons Only" = quote(
    PARAM == "Reason for Treatment Discontinuation" &
      APHASEN == 1 &
      FAS1FL == "Y" &
      (AVALC %in%
        c("LOST TO FOLLOW-UP", "WITHDRAWAL BY SUBJECT", "OTHER") |
        CRIT9FL == "Y")
  ),
  # Exclude patients from delta adjustment with the following disposition information:
  # PARAM = 'Reason for Treatment Discontinuation' and APHASEN = 1 and FAS1FL = 'Y';
  # & AVALC in ('LOST TO FOLLOW-UP','WITHDRAWAL BY SUBJECT','OTHER');
  "All Subjects Who Discontinued and/or Switched Study Drug and/or Underlying Antidepressant for Study Agent Reasons Only Including Reasons Due to COVID-19" = quote(
    PARAM == "Reason for Treatment Discontinuation" &
      APHASEN == 1 &
      FAS1FL == "Y" &
      (AVALC %in% c("LOST TO FOLLOW-UP", "WITHDRAWAL BY SUBJECT", "OTHER"))
  )
)

# Number of cores to parallelize the computations on, because they take a lot of time.
ncores <- 4

# Double blind labeling in the output or not.
doubleblind <- TRUE

# Define confidence level to use.
conflvl <- 0.95

# Derived formats specifications.
formats <- list(
  pval = jjcsformat_pval_fct(alpha),
  diff_mean_se = jjcsformat_xx("xx.x (xx.xx)"),
  diff_mean_ci = jjcsformat_xx("(xx.xx, xx.xx)")
)
dblabel <- if (doubleblind) " (DB)" else ""

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  mutate(!!trtvar := forcats::fct_relevel(!!rlang::sym(trtvar), ctrlab)) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl), any_of(covariates))

## ADMADRSI ----

admadrsi <- haven::read_sas(read_path(a_in, "admadrsi.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  filter(PARAMCD == resppar) |>
  filter(!!domfl) |>
  select(USUBJID, AVISIT, AVAL, BASE, CHG) |>
  mutate(
    USUBJID = factor(USUBJID),
    AVISIT = factor(AVISIT)
  )

## ADICE ----

adice <- haven::read_sas(read_path(a_in, "adice.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y"))

## ADDISP ----

addisp <- haven::read_sas(read_path(a_in, "addisp.sas7bdat"))

## Analysis ----

# For multiple imputation.
ana_mi <- admadrsi |>
  inner_join(adsl, by = "USUBJID") |>
  # Complete patient/visits grid.
  expand_locf(
    USUBJID = levels(admadrsi$USUBJID),
    AVISIT = levels(admadrsi$AVISIT),
    vars = c("STUDYID", popfl, trtvar, covariates),
    group = c("USUBJID"),
    order = c("USUBJID", "AVISIT")
  ) |>
  # Drop baseline visit (because there is no change to model yet).
  filter(AVISIT != basetpt) |>
  mutate(AVISIT = droplevels(AVISIT))

# Multiple Imputation ----

## Intercurrent events ----

# This could look more complicated, when different strategies are used
# for different events. Here we follow the approach that
# - intermittent missing values are imputed using hypothetical strategy (MAR)
# - final missing values are imputed using the copy-reference (CR) strategy
#   and therefore in `ana_ice` the first such missing value is defined.
# All complete patients, or patients who only have intermittent missing values
# are therefore not included in `ana_ice`.
ana_ice <- ana_mi |>
  select(USUBJID, AVISIT, CHG) |>
  arrange(USUBJID, AVISIT) |>
  group_by(USUBJID) |>
  summarize(AVISIT = find_missing_chg_after_avisit(pick(AVISIT, CHG))) |>
  ungroup() |>
  filter(!is.na(AVISIT)) |>
  mutate(AVISIT = factor(AVISIT, levels = levels(ana_mi$AVISIT))) |>
  mutate(strategy = reference_method)

## Bayesian imputation ----

vars <- rbmi::set_vars(
  outcome = "CHG",
  visit = "AVISIT",
  subjid = "USUBJID",
  group = trtvar,
  covariates = c(covariates, interactions)
)

set.seed(12345)
wrapper <- if (testrun) suppressWarnings else identity
draw_obj <- wrapper(rbmi::draws(
  data = ana_mi,
  data_ice = ana_ice, # Here we pass the intercurrent events information already.
  vars = vars,
  method = imputation_method,
  quiet = TRUE
))
impute_obj <- rbmi::impute(
  draw_obj,
  references = reference_spec
)

## Tipping point analyses ----

# Obtain the subset of subjects where the delta adjustment will be applied,
# cf. the the definition in Parameters section.
flag_to_subset <- function(delta_addisp_exclude_flag) {
  adice_ids <- adice |>
    filter(!!rlang::sym(trtvar) == trtlab) |>
    pull(USUBJID) |>
    unique()
  addisp_exclude_ids <- addisp |>
    filter(!!delta_addisp_exclude_flag) |>
    pull(USUBJID) |>
    unique()
  setdiff(adice_ids, addisp_exclude_ids)
}

# Do this for all subsets.
delta_addisp_subsets <- lapply(delta_addisp_exclude_flags, flag_to_subset)

# Find the unique subsets, because only there the results can differ.
delta_addisp_unique_subsets <- unique(delta_addisp_subsets)

# Match the original subsets to the unique subsets.
delta_addisp_subset_map <- match(
  delta_addisp_subsets,
  delta_addisp_unique_subsets
)

# Function to perform tipping point analysis with one delta subset definition.
analyze_delta_subset <- function(delta_subset_ids) {
  # Get initial delta data frame template.
  delta_df_init <- delta_template(impute_obj) |>
    mutate(apply_delta = is_missing & (USUBJID %in% delta_subset_ids))

  # Initialize loop over delta values.
  pval <- 0
  delta_index <- 0

  while (pval < alpha & delta_index < length(deltas)) {
    ### Delta definition ----

    delta_index <- delta_index + 1
    delta_used <- deltas[delta_index]
    delta_df <- delta_df_init |>
      mutate(delta = delta_used * apply_delta)

    ### Analysis step ----

    ana_obj <- analyse(
      impute_obj,
      fun = rbmi_mmrm,
      delta = delta_df,
      vars = set_vars(
        subjid = "USUBJID",
        outcome = "CHG",
        visit = "AVISIT",
        group = trtvar,
        covariates = c(covariates, interactions)
      ),
      cov_struct = "us",
      weights = "equal",
      reml = TRUE,
      method = "Kenward-Roger",
      vcov = "Kenward-Roger-Linear",
      ncores = ifelse(testrun, 1, cl)
    )

    ### Pooling ----

    pool_obj <- junco:::rbmi_pool(
      ana_obj,
      conf.level = conflvl,
      alternative = switch(
        pval_sided,
        "2" = "two.sided",
        "-1" = "less",
        "1" = "greater"
      )
    )

    ### p-value extraction ----

    ana_mi_results <- broom::tidy(pool_obj, visits = levels(ana_mi$AVISIT)) |>
      filter(visit == timept) |>
      rename(!!trtvar := group) |>
      filter(!!rlang::sym(trtvar) == trtlab) |>
      mutate(!!trtvar := paste0(trtlab, " vs. ", ctrlab))

    pval <- ana_mi_results |>
      pull(p_value)
  }

  list(
    pval = pval,
    delta_used = delta_used,
    ana_mi_results = ana_mi_results
  )
}

## Execute tipping point analyses ----

# Start cluster to speed up computations.
if (!testrun) {
  cl <- rbmi::make_rbmi_cluster(
    ncores = ncores,
    objects = list(
      rbmi_mmrm = rbmi_mmrm,
      rbmi_mmrm_single_info = rbmi_mmrm_single_info
    ),
    packages = c("mmrm", "emmeans", "checkmate")
  )
}

# Perform computations for unique subsets.
delta_unique_subset_results <- lapply(
  delta_addisp_unique_subsets,
  analyze_delta_subset
)

# Now we don't need the cluster R processes any longer.
if (!testrun) {
  parallel::stopCluster(cl)
}

# Map back to original (potentially duplicated) subsets.
delta_subset_results <- delta_unique_subset_results[delta_addisp_subset_map]

# Layout ----

## Multiple imputation results ----

lyt_base <- rtables::basic_table(
  top_level_section_div = " ",
  show_colcounts = FALSE
) |>
  split_cols_by(trtvar)

get_table <- function(delta_subset_result, delta_subset_label) {
  # Because the layout here includes result information (namely the delta_used)
  # we need to include the RBMI tabulation part of the layout in here.
  lyt_base |>
    analyze(
      trtvar,
      afun = a_rbmi_lsmeans,
      na_str = default_na_str(),
      var_labels = paste0(
        "Subset of subjects with delta adjusted imputed values: ",
        delta_subset_label,
        " ",
        ctrlab,
        "-adjustment=0; ",
        trtlab,
        "-adjustment=0 to Δ*~[super a]"
      ),
      show_labels = "visible",
      extra_args = list(
        .stats = c(
          "additional_title_row",
          "diff_mean_se",
          "diff_mean_ci",
          "p_value"
        ),
        .labels = c(
          additional_title_row = paste0(
            ctrlab,
            "-adjustment=0; ",
            trtlab,
            "-adjustment=",
            delta_subset_result$delta_used
          ),
          diff_mean_se = paste0(
            "Diff. of LS means (SE) (",
            trtlab,
            " minus ",
            ctrlab,
            ")~[super b]"
          ),
          diff_mean_ci = paste0(
            tern::f_conf_level(conflvl),
            " confidence interval on diff."
          ),
          p_value = paste0(pval_sided, "-sided p-value")
        ),
        .formats = c(
          p_value = formats$pval,
          diff_mean_se = formats$diff_mean_se,
          diff_mean_ci = formats$diff_mean_ci
        ),
        .indent_mods = c(diff_mean_se = 1L),
        ref_path = c(trtvar, ctrlab)
      )
    ) |>
    build_table(delta_subset_result$ana_mi_results)
}

# Output ----

subset_tables <- mapply(
  get_table,
  delta_subset_result = delta_subset_results,
  delta_subset_label = names(delta_addisp_exclude_flags)
)
result <- do.call(rbind, subset_tables)


# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
