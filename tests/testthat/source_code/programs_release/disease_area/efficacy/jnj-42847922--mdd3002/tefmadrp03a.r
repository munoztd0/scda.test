###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefmadrm03a.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefmadrm03a: [MADRS] Subjects Who Achieved Remission Over Time; Full Analysis Set
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADMADRSI.
## Output:                    TEFMADRM03a.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:                     1
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
## Description:
################################################################################

################################################################################
# Prep Environment
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

################################################################################
# - Define output ID and file location
# - Define treatment variable used (default=TRT01P)
# - Define population flag used (default=SDBFL)
# - Define control treatment arm
# - Define parameter code PARAMCD needed for output
################################################################################

tblid <- "TEFMADRP03a"
fileid <- write_path(opath, tblid)

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

trtvar <- "TRT01P"
popfl <- "FAS2FL"

ctrl_grp <- "Placebo"
# param_code <- "MADRESP"
param_code <- "MADR0113"

################################################################################
# Read in ADSL and ADMADRS datasets
################################################################################
adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y" & !!rlang::sym(trtvar) != "") |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl))

adsl$colspan_trt <- factor(
  ifelse(adsl[[trtvar]] == ctrl_grp, " ", "Active Study Agent"),
  levels = c("Active Study Agent", " ")
)

# Observed data
admadrs_o <- haven::read_sas(read_path(a_in, "admadrsi.sas7bdat")) |>
  # remove/update ANLxxFL='Y' as required
  filter(ANL02FL == "Y" & PARAMCD == param_code & APHASEN == 1 & DTYPE == "") |>
  mutate(
    PARAMCD = as.factor(PARAMCD),
    AVISIT = as.factor(AVISIT)
  ) |>
  mutate(DTYPE = "OBSERVED") |>
  select(USUBJID, PARAMCD, ANL02FL, AVISIT, APHASEN, AVALC, DTYPE)

# Observed data
admadrs_o <- haven::read_sas(read_path(a_in, "admadrsi.sas7bdat")) |>
  filter(ANL02FL == "Y" & PARAMCD == param_code & APHASEN == 1 & DTYPE == "") |>
  mutate(
    PARAMCD = as.factor(PARAMCD),
    AVISIT = as.factor(AVISIT),
    DTYPE = "OBSERVED",
    DTYPE = factor(DTYPE, levels = c("OBSERVED")),
    dtypelbl = "Observed case"
  ) |>
  # adapt below as needed
  mutate(AVISIT = factor(AVISIT, levels = c("Day 15", "Day 29", "Day 43"))) |>
  select(USUBJID, PARAMCD, ANL02FL, AVISIT, APHASEN, AVALC, DTYPE, dtypelbl)

# Imputed data
admadrs_i <- haven::read_sas(read_path(a_in, "admadrsi.sas7bdat")) |>
  filter(ANL02FL == "Y" & PARAMCD == param_code & APHASEN == 1) |>
  mutate(
    PARAMCD = as.factor(PARAMCD),
    AVISIT = as.factor(AVISIT),
    DTYPE = "IMPUTED",
    DTYPE = factor(DTYPE, levels = c("IMPUTED")),
    dtypelbl = "Imputed~[super a]"
  ) |>
  # adapt below as needed
  mutate(AVISIT = factor(AVISIT, levels = c("Day 15", "Day 29", "Day 43"))) |>
  select(USUBJID, PARAMCD, ANL02FL, AVISIT, APHASEN, AVALC, DTYPE, dtypelbl)

# join data together
madrs1 <- admadrs_o |> inner_join(adsl, by = c("USUBJID"))
madrs2 <- admadrs_i |> inner_join(adsl, by = c("USUBJID"))

colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = ctrl_grp,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

################################################################################
# Define layout and build table:
################################################################################

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
  split_rows_by("DTYPE", labels_var = "dtypelbl") |>
  split_rows_by("AVISIT") |>
  analyze(
    "AVISIT",
    table_names = "AVAL1x",
    show_labels = "hidden",
    indent_mod = 0L,
    afun = function(x) {
      list(
        "N" = rcell(length(x), format = jjcsformat_xx("xx"))
      )
    }
  ) |>
  analyze(
    "AVALC",
    afun = a_freq_j,
    extra_args = list(
      denom = "n_df",
      val = "Y",
      label = "Subjects with response",
      .stats = c("n_df", "count_unique_fraction")
    ),
    show_labels = "hidden",
    indent_mod = 0L
  )

result1 <- build_table(lyt, madrs1, alt_counts_df = adsl)
result2 <- build_table(lyt, madrs2, alt_counts_df = adsl)

result <- rbind(result1, result2)

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

################################################################################
# Convert to tbl file and output table
################################################################################
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
