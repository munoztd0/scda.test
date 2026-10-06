###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefcyto04a.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefcyto04a: Overall Best Confirmed
##                            Response Based on Computerized Algorithm by Cytogenetic Risk Group at
##                            Baseline for Subjects With [DRd]; Full Analysis Set (Study mmy, bc, lc)
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADEFF
## Output:                    TEFCYTO04a.rtf
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
tblid <- "TEFCYTO04a"
fileid <- write_path(opath, tblid)

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
tab_titles$main_footer <- tab_titles$main_footer[c(1, 4, 5, 6)]

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define treatment variable grouping to be used.
grouping <- list(
  "Dummy A" = c("Dummy A - Tec-Dara"),
  "Dummy B" = c("Dummy B - DPd", "Dummy B - DVd")
)

# Define control among groups defined above.
ctrlab <- "Dummy B"

# Define population flag used (here take measurable disease at baseline flag).
popfl <- "ITTFL"

# Define subgroup variable to be used.
# Note: This should be chosen as the cytogenic risk variable.
subgrvar <- "RACEGR1"

# Define subgroup label to be used.
subgrlbl <- "Race"

# Define the strata variables to use below.
strata <- c("STRAT02A", "STRAT03")
# For unstratified analyses, please use:
# strata <- NULL

# Define response parameter to be used.
resppar <- "IRCBRESP"

# Define method for response rate confidence intervals.
# See ?tern::s_proportion for possible options.
cimethod <- "clopper-pearson"

# Define method for comparing response.
method <- "or_cmh"
# one of:
# "rr": relative risk
# "or_logistic": odds ratio by logistic regression
# "or_cmh": odds ratio by Cochran-Mantel-Haenszel (CMH)

# If unstratified, then this choice of the test matters:
testmethod <- "fisher"
# one of:
# "fisher": Fisher's exact test
# "chisq": Chi-Square test

# Define significance threshold to use (important for p-value formatting).
# 0 means no formal testing is applied, therefore standard p-value rounding applies.
alpha <- 0

# Define confidence level to use.
conflvl <- 0.95

# Derived formats and methods specifications.
formats <- list(
  comp_stat_ci = jjcsformat_xx("xx.xx (xx.xx, xx.xx)"),
  pval = jjcsformat_pval_fct(alpha),
  prop_ci = jjcsformat_xx("(xx.x%, xx.x%)")
)

methods <- list(
  comp_stat_ci = method,
  pval = testmethod,
  prop_ci = cimethod
)

# Data ----

## ADSL ----

grpvar <- "GRPVAR"

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  mutate(
    !!grpvar := do.call(
      forcats::fct_collapse,
      c(list(.data[[trtvar]]), grouping)
    )
  ) |>
  mutate(
    !!grpvar := forcats::fct_relevel(.data[[grpvar]], ctrlab, after = Inf)
  ) |>
  mutate(!!subgrvar := as.factor(.data[[subgrvar]])) |>
  select(
    STUDYID,
    USUBJID,
    all_of(grpvar),
    all_of(subgrvar),
    all_of(popfl),
    all_of(strata)
  )

## ADEFF ----

adeff <- haven::read_sas(read_path(a_in, "adeff.sas7bdat")) |>
  select(USUBJID, PARAMCD, AVALC) |>
  filter(PARAMCD == resppar) |>
  select(-PARAMCD) |>
  mutate(
    USUBJID = factor(USUBJID),
    AVALC = forcats::fct_recode(
      AVALC,
      "sCR" = "Stringent Complete Response (sCR)",
      "CR" = "Complete Response (CR)",
      "VGPR" = "Very Good Partial Response (VGPR)",
      "PR" = "Partial Response (PR)",
      "MR" = "Minimal Response (MR)",
      "SD" = "Stable Disease (SD)",
      "PD" = "Progressive Disease (PD)",
      "NE" = "Not Evaluable (NE)"
    )
  )

## Analysis ----
ana <- adeff |>
  mutate(
    rsp_lab = jjcs_mmy_rsp_label(AVALC),
    is_rsp_any = AVALC %in% c("sCR", "CR", "VGPR", "PR"),
    is_rsp_cbr = AVALC %in% c("sCR", "CR", "VGPR", "PR", "MR"),
    is_rsp_vgpr_or_better = AVALC %in% c("sCR", "CR", "VGPR"),
    is_rsp_cr_or_better = AVALC %in% c("sCR", "CR")
  ) |>
  right_join(adsl, by = "USUBJID")

# Functions ----

analyze_response_logical <- function(lyt, vars, label, include_comp) {
  analyze(
    lyt = lyt,
    vars = vars,
    afun = resp01_acfun,
    show_labels = "hidden",
    extra_args = list(
      arm = grpvar,
      include_comp = include_comp,
      conf_level = conflvl,
      strata = strata,
      label = label,
      methods = methods,
      formats = formats
    )
  )
}

# Layout ----

lyt <- basic_table() |>
  split_cols_by(
    grpvar,
    show_colcounts = TRUE,
    # Note: It needs to stay "Overall" here (this won't be shown anyway).
    split_fun = add_overall_level("Overall", label = "", first = FALSE)
  ) |>
  split_cols_by(
    "STUDYID", # The particular choice of this is irrelevant.
    split_fun = resp01_split_fun_fct(method = method, conf_level = conflvl)
  ) |>
  ## Subgroups ----
  split_rows_by(
    subgrvar,
    page_by = TRUE,
    split_fun = drop_split_levels
  ) |>
  summarize_row_groups(
    cfun = resp01_counts_cfun,
    extra_args = list(label_fstr = paste0(subgrlbl, ": %s"))
  ) |>
  ## Response category ----
  analyze(
    vars = "rsp_lab",
    afun = resp01_acfun,
    show_labels = "visible",
    var_labels = "Response category",
    extra_args = list(
      arm = grpvar,
      include_comp = c(),
      conf_level = conflvl,
      strata = strata,
      formats = formats,
      methods = methods
    )
  ) |>
  insert_blank_line() |>
  ## Overall response ----
  analyze_response_logical(
    vars = "is_rsp_any",
    label = "Overall response (sCR+CR+VGPR+PR)",
    include_comp = FALSE
  ) |>
  ## Clinical benefit rate ----
  analyze_response_logical(
    vars = "is_rsp_cbr",
    label = "Clinical benefit rate (overall response+MR)",
    include_comp = FALSE
  ) |>
  ## VGPR or better ----
  analyze_response_logical(
    vars = "is_rsp_vgpr_or_better",
    label = "VGPR or better (sCR+CR+VGPR)",
    include_comp = TRUE
  ) |>
  ## CR or better ----
  analyze_response_logical(
    vars = "is_rsp_cr_or_better",
    label = "CR or better (sCR+CR)",
    include_comp = TRUE
  ) |>
  append_topleft("Response")

# Output ----

result <- build_table(lyt, ana, alt_counts_df = adsl)

# Post-hoc we can suppress the column count for the Overall column.
colcount_visible(result, c(grpvar, "Overall")) <- FALSE

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
