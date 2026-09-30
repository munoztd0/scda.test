###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefiga02.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefiga02: IGA Response Analysis
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     adsl.sas7bdat, adigrspi.sas7bdat
## Output:                    tefiga02.rtf
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

tblid <- "TEFIGA02"
fileid <- write_path(opath, tblid)
popfl <- "FASFL"
trtvar <- "TRT01P"

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

timepoint <- "Week 16"
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

adrsp <- haven::read_sas(read_path(a_in, "adigrspi.sas7bdat")) |>
  filter(
    !!rlang::sym(popfl) == "Y" &
      AVISIT == timepoint &
      PARAMCD %in% c("IGMLDP", "IGMINP", "IGCLRP")
  ) |>
  mutate(
    AVISIT = forcats::fct_reorder(factor(AVISIT), AVISITN),
    response = case_when(
      AVALC == "Y" ~ TRUE,
      TRUE ~ FALSE
    )
  ) |>
  select(USUBJID, AVISIT, PARAMCD, AVALC, response)

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
  split_rows_by("PARAMCD", child_labels = "hidden", section_div = c(" ")) |>
  estimate_proportion(
    vars = "response",
    table_names = "est_prop",
    .stats = c("n_prop"),
    .labels = c("n_prop" = "n_prop"),
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

obj_label(tt_at_path(result, c("PARAMCD", "IGCLRP", "est_prop", "n_prop"))) <-
  "IGA of cleared (0)"

obj_label(tt_at_path(result, c("PARAMCD", "IGMINP", "est_prop", "n_prop"))) <-
  "IGA of cleared (0), or minimal (1)"

obj_label(tt_at_path(result, c("PARAMCD", "IGMLDP", "est_prop", "n_prop"))) <-
  "IGA of cleared (0), minimal (1), or mild (2)"

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)
################################################################################
# Convert to tbl file and output table:
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
