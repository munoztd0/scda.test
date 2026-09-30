###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefdlqi01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefdlqi01: DLQI Response Analysis
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     adsl.sas7bdat, addlrspi.sas7bdat
## Output:                    tefdlqi01.rtf
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

tblid <- "TEFDLQI01"
fileid <- write_path(opath, tblid)
popfl <- "FASFL"
trtvar <- "TRT01P"

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

timepoint <- "Week 16"
resp_paramcd <- "DLLE1P"
resp_txt <- "Subjects with DLQI of 0 or 1"
bsl_var <- "BLGT1FL"
bsl_txt <- "Subjects with DLQI > 1 at baseline"
stratvar <- "STRATWTG"

################################################################################
# Process data:
################################################################################

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(
    !!popfl := factor(!!rlang::sym(popfl)),
    !!trtvar := factor(
      !!rlang::sym(trtvar),
      levels = c(
        "JNJ-77242113 25 MG QD",
        "JNJ-77242113 50 MG QD",
        "JNJ-77242113 25 MG BID",
        "JNJ-77242113 100 MG QD",
        "JNJ-77242113 100 MG BID",
        "PLACEBO"
      )
    )
  ) |>
  create_colspan_var(
    non_active_grp = "PLACEBO",
    non_active_grp_span_lbl = " ",
    active_grp_span_lbl = "Active Study Agent",
    colspan_var = "colspan_trt",
    trt_var = trtvar
  ) |>
  select(
    USUBJID,
    !!rlang::sym(popfl),
    !!rlang::sym(trtvar),
    colspan_trt,
    !!rlang::sym(stratvar)
  )

adrsp <- haven::read_sas(read_path(a_in, "addlrspi.sas7bdat")) |>
  filter(
    !!rlang::sym(popfl) == "Y" &
      AVISIT == timepoint &
      PARAMCD == resp_paramcd &
      !!rlang::sym(bsl_var) == "Y"
  ) |>
  mutate(
    !!bsl_var := !!rlang::sym(bsl_var) == "Y",
    AVISIT = forcats::fct_reorder(factor(AVISIT), AVISITN),
    response = case_when(
      AVALC == "Y" ~ TRUE,
      TRUE ~ FALSE
    )
  ) |>
  select(USUBJID, AVISIT, PARAMCD, AVALC, response, !!rlang::sym(bsl_var))

adrsp <- inner_join(x = adsl, y = adrsp, by = "USUBJID")

################################################################################
# Define layout and build table:
################################################################################

# Map each treatment group to appropriate columns spanning header.
colspan_trt_map <- create_colspan_map(
  df = adsl,
  trt_var = trtvar,
  colspan_var = "colspan_trt",
  non_active_grp = "PLACEBO",
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent"
)

# Specification for reference column
ref_path <- c("colspan_trt", " ", trtvar, "PLACEBO")

lyt <- basic_table(show_colcounts = TRUE, colcount_format = "N=xx") |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  analyze(
    vars = bsl_var,
    show_labels = "hidden",
    afun = function(df) {
      in_rows(
        sum(df[[bsl_var]]),
        .formats = "xx.",
        .labels = bsl_txt
      )
    }
  ) |>
  insert_blank_line() |>
  estimate_proportion(
    vars = "response",
    table_names = "est_prop",
    .stats = c("n_prop"),
    .labels = c("n_prop" = resp_txt),
    .formats = c("n_prop" = jjcsformat_count_fraction)
  ) |>
  analyze(
    vars = "response",
    afun = a_test_proportion_diff,
    table_names = "pval",
    show_labels = "hidden",
    na_str = default_na_str(),
    extra_args = list(
      method = "cmh",
      variables = list(strata = stratvar),
      ref_path = ref_path,
      .labels = c(pval = "p-value")
    )
  )

result <- build_table(lyt, adrsp, alt_counts_df = adsl)

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
