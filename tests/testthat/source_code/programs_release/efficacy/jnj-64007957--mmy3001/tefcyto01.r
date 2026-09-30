###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefcyto01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefcyto01: Overall Best Confirmed
##                            Response Based on Computerized Algorithm by Cytogenetic High Risk Marker;
##                            [Analysis Set] Analysis Set (Study mmy, bc, lc)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADEFF
## Output:                    TEFCYTO01.rtf
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
tblid <- "TEFCYTO01"
fileid <- write_path(opath, tblid)

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

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
  mutate(is_rsp = AVALC %in% c("sCR", "CR", "VGPR", "PR")) |>
  right_join(adsl, by = "USUBJID")

# Functions ----

analyze_response_logical <- function(lyt, vars, label, include_comp = TRUE) {
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
  ## Overall response ----
  analyze_response_logical(
    vars = "is_rsp",
    label = paste0("Overall response rate by ", subgrlbl)
  ) |>
  split_rows_by(subgrvar, indent_mod = 1) |>
  summarize_row_groups(
    var = "is_rsp",
    cfun = resp01_acfun,
    extra_args = list(
      arm = grpvar,
      include_comp = TRUE,
      conf_level = conflvl,
      strata = strata,
      methods = methods,
      formats = formats
    )
  )

# Output ----

result <- build_table(lyt, ana, alt_counts_df = adsl)

# Post-hoc we can suppress the column count for the Overall column.
colcount_visible(result, c(grpvar, "Overall")) <- FALSE

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
