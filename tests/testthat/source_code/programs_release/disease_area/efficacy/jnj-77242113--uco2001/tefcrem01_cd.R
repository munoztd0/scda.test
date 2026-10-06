###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefcrem01_cd.r
## R Version:                 4.5.2
## junco Version:             0.1.6
## Short Description:         Program to create tefcrem01_cd:
##                            Subjects in Clinical Remission
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adibdqi
## Output:                    tefcrem01_cd.rtf
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
tblid <- "TEFCREM01_CD"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "Placebo"

# Define the order of treatment group labels in the treatment variable.
trtlab <- c("JNJ-77242113 100 mg QD", "JNJ-77242113 200 mg QD", "JNJ-77242113 400 mg QD")

# Define population flags used.
popfl <- expr(FASFL == "Y" & .data[[trtvar]] != "" & !is.na(.data[[trtvar]]))

# Define the stratification variables to use for CMH.
stratvar <- c("STRAT01", "STRAT02")
# For unstratified analyses, please use:
# strata <- NULL

# Define domain level flags.
domfl <- quote(EST01RFL == "Y")

# Define visits from which to define sustained response.
timepoints_response <- c("Week 12", "Week 28")
last_timepoint_response <- tail(timepoints_response, 1L)

# Define response parameter to be used.
# along with its label to use in the table.
resppar <- "IBDQREM"
resplbl <- paste0(
  "Subjects with sustained IBDQ remission at ",
  last_timepoint_response,
  "~[super a,b,c]"
)

# Proportion estimation method.
prop_method <- "wald / clopper-pearson"
# based on one of the following methods:
# "clopper-pearson" - Clopper-Pearson exact confidence intervals
# "wald" - normal approximation, i.e. Wald statistic
# "wald / clopper-pearson" - conditionally either use Wald or Clopper-Pearson.
# See ?tern::estimate_proportion for details.
#
# The Clopper-Pearson method is used when any of these are true:
# (a) the number of responders is `num_limit` (or less),
# (b) all subjects except `num_limit` (or less) have observed response,
# or (c) the observed group size is less than `denom_limit`.
# See ?junco::cond_proportion_j for details.

# Default choices, only needed for "wald / clopper-pearson" as described above:
num_limit <- 0
denom_limit <- 10

# Define proportion difference estimation method.
prop_diff_method <- "cmh_sato"
# based on one of the following statistics:
# "cmh" - CMH test
# "cmh_sato" - CMH test with Sato variance estimator
# "wald" - Wald
# see ?tern::prop_diff for all possible options.

# Define test method.
pval_method <- "cmh_sato"
# one of:
# "cmh" - CMH test
# "cmh_sato" - CMH test with Sato variance estimator
# "chisq" - Chi-square test
# see ?tern::test_proportion_diff for all possible options.

# Define between one and two sided p-values.
pval_sided <- "2"
# one of:
# "2" - two sided
# "-1" - one sided, less
# "1" - one sided, greater

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
# If a value larger than 0 is chosen, then no rounding is applied for values just below:
# For example, 0.0048 is not rounded to 0.005 but stays at 0.0048 if alpha = 0.005 is set.
# Please see ?jjcsformat_pval_fct for more details and examples.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived formats specifications.
formats <- list(
  est_prop = jjcsformat_count_fraction,
  prop_ci = jjcsformat_xx("(xx.x%, xx.x%)"),
  est_prop_diff = jjcsformat_xx("xx.x% (xx.x%, xx.x%)"),
  pval = jjcsformat_pval_fct(alpha)
)

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!popfl) |>
  mutate(!!trtvar := factor(.data[[trtvar]], levels = c(ctrlab, trtlab))) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(stratvar)) |>
  mutate(across(all_of(stratvar), ~ factor(.x)))

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

## ADIBDQI ----

adibdqi <- haven::read_sas(read_path(a_in, "adibdqi.sas7bdat")) |>
  filter(!!popfl) |>
  filter(!!domfl) |>
  filter(
    PARAMCD == resppar,
    AVISIT %in% timepoints_response
  ) |>
  select(STUDYID, USUBJID, AVISIT, AVALC)

## Analysis ----

ana <- adibdqi |>
  # Define sustained remission when all the timepoints considered have remission.
  summarize(AVALC = if (all(AVALC == "Y")) "Y" else "N", .by = c(STUDYID, USUBJID)) |>
  # Right join here because missing patients need to be counted as non-responders.
  right_join(adsl, by = c("STUDYID", "USUBJID")) |>
  mutate(
    # Translate to logical response variable.
    response = !is.na(AVALC) & AVALC == "Y",
    visit = last_timepoint_response
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
  split_rows_by("visit") |>
  # Note: We currently need to use the tern layout function here, just using
  # analyze() with a_proportion does not work.
  analyze(
    "STUDYID",
    var_labels = "N",
    show_labels = "hidden",
    afun = function(df) nrow(df),
    format = "xx"
  )

lyt <- if (prop_method == "wald / clopper-pearson") {
  lyt |>
    analyze(
      vars = "response",
      afun = a_cond_proportion_j,
      table_names = "est_prop",
      show_labels = "hidden",
      extra_args = list(
        conf_level = conflvl,
        num_limit = num_limit,
        denom_limit = denom_limit,
        .stats = c("n_prop", "prop_ci"),
        .labels = c(
          "n_prop" = resplbl,
          "prop_ci" = paste0(tern::f_conf_level(conflvl), "~[super d]")
        ),
        .formats = c(
          "n_prop" = formats$est_prop,
          "prop_ci" = formats$prop_ci
        ),
        .indent_mods = 1
      )
    )
} else {
  lyt |>
    estimate_proportion(
      "response",
      table_names = "est_prop",
      show_labels = "hidden",
      method = prop_method,
      .stats = c("n_prop", "prop_ci"),
      .labels = c(
        "n_prop" = resplbl,
        "prop_ci" = paste0(tern::f_conf_level(conflvl), "~[super d]")
      ),
      .formats = c(
        "n_prop" = formats$est_prop,
        "prop_ci" = formats$prop_ci
      ),
      .indent_mods = 1
    )
}

lyt <- lyt |>
  analyze(
    vars = "response",
    afun = a_proportion_diff_j,
    na_str = default_na_str(),
    table_names = "est_prop_diff",
    show_labels = "hidden",
    extra_args = list(
      .stats = "diff_est_ci",
      .labels = c(
        diff_est_ci = paste0(
          ifelse(!is.null(stratvar), "Adjusted treatment", "Treatment"),
          " difference (",
          tern::f_conf_level(conflvl),
          ")~[super e]"
        )
      ),
      .formats = c(diff_est_ci = formats$est_prop_diff),
      .indent_mods = 1,
      method = prop_diff_method,
      variables = list(strata = stratvar),
      conf_level = conflvl,
      ref_path = ref_path
    )
  ) |>
  analyze(
    vars = "response",
    afun = junco::a_test_proportion_diff,
    table_names = "pval",
    show_labels = "hidden",
    na_str = default_na_str(),
    extra_args = list(
      method = pval_method,
      alternative = switch(
        pval_sided,
        "2" = "two.sided",
        "-1" = "less",
        "1" = "greater"
      ),
      variables = list(strata = stratvar),
      .formats = c("pval" = formats$pval),
      .indent_mods = 2,
      ref_path = ref_path,
      .labels = c(pval = "p-value~[super e]")
    )
  )

# Output ----

result <- build_table(lyt, ana, alt_counts_df = adsl, round_type = "sas")

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
