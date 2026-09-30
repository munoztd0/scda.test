###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefeq5d01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefeq5d01:
##                              [Primary/Secondary Endpoint Analysis] ([Composite/Main Estimand])
##                              Change From [Induction] Baseline in Health State EQ-5D VAS Score
##                              and EQ-5D Dimensions at [Time Point];  Analysis Set (Study gi)
## Author:                    Technology Solutions
## Date:                      2026-09-30
## Input:                     ADSL, ADEFF
## Output:                    TEFEQ5D01.rtf
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
tblid <- "TEFEQ5D01"
fileid <- write_path(opath, tblid)

# Please select titles and footnotes as appropriate.
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define treatment arms to use.
trtlvls <- c(
  "Guselkumab 200 mg IV -> Guselkumab 100 mg SC q8w",
  "Guselkumab 200 mg IV -> Guselkumab 200 mg SC q4w",
  "Ustekinumab 6 mg/kg"
)

# Define control group label used in the treatment variable.
ctrlab <- "Placebo"

# Define population flag used.
popfl <- "PASFL"

# Define EQ5D dimension parameters to use.
paramcds <- c(
  "Mobility" = "EQ5DMOBP",
  "Self-Care" = "EQ5DSCP",
  "Usual Activities" = "EQ5DUAP",
  "Pain/Discomfort" = "EQ5DPDP",
  "Anxiety/Depression" = "EQ5DADP"
)

# Define EQ5D VAS parameter to use.
vasparam <- "EQ5DVASP"

# Define visits to use. The first one needs to be the baseline visit.
# Use "Induction" here e.g. to use it consistently throughout the table.
visits <- c("Baseline", "Week 12", "Week 24", "Week 48")

# Define the strata variables to use below.
strata <- c("RASTRA1", "RASTRA2", "RASTRA3", "RASTRA4")
# For unstratified analyses, please use:
# strata <- NULL

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
  iqr = jjcsformat_xx("xx.x, xx.x"),
  lsmean_se = jjcsformat_xx("xx.x (xx.xx)"),
  lsmean_diffci = jjcsformat_xx("xx.x (xx.xx, xx.xx)"),
  pval = jjcsformat_pval_fct(alpha)
)

# Data ----

## ADSL ----

adsl <- read_sas(file.path(a_in[[environ]], "adsl.sas7bdat")) |>
  filter(
    !!rlang::sym(popfl) == "Y",
    !!rlang::sym(trtvar) %in% c(ctrlab, trtlvls)
  ) |>
  mutate(!!trtvar := factor(.data[[trtvar]], levels = c(ctrlab, trtlvls))) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(strata))

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

## ADEQ5D ----

adeq5d_full <- read_sas(file.path(a_in[[environ]], "adeq5d.sas7bdat"))
adeq5d <- adeq5d_full |>
  filter(!!rlang::sym(popfl) == "Y") |>
  filter(
    toupper(PARAMCD) %in%
      c(paramcds, vasparam) &
      (toupper(AVISIT) %in% toupper(visits) | ABLFL == "Y") &
      round(AVISITN) == AVISITN &
      !is.na(AVAL) &
      ANL03FL == "Y"
  ) |>
  select(
    USUBJID,
    PARAMCD,
    PARAM,
    AVAL,
    AVALC,
    CHG,
    CHGCAT1,
    ANL03FL,
    AVISIT,
    AVISITN,
    ABLFL
  ) |>
  mutate(
    AVISIT = case_when(
      ABLFL == "Y" ~ visits[1],
      TRUE ~ stringr::str_to_sentence(AVISIT)
    )
  )
stopifnot(all(adeq5d$AVISIT %in% visits))

### Label map ----

# Construct a mapping table with all levels of the parameters.
# We use the entire dataset to maximize the chance that all possible levels are included.
# (Note that an alternative approach could be to work with an Excel file that covers all
# levels for all parameters.)
vals <- adeq5d_full |>
  filter(toupper(PARAMCD) %in% paramcds) |>
  filter(!is.na(AVAL)) |>
  select(PARAMCD, PARAM, AVAL, AVALC) |>
  arrange(PARAMCD, AVAL) |>
  unique()

valmax <- vals |>
  group_by(PARAMCD) |>
  summarize(max = max(AVAL))
# Ensure that all parameters have 5 levels.
# Otherwise this has to be handled on study level.
stopifnot(all(valmax$max == 5))

# Convert dataframe into a `label_map` tibble that can be used with the `a_freq_j`
# function.
xlabel_map <- vals |>
  mutate(value = AVALC, label = tolower(AVALC), var = "AVALC") |>
  select(PARAMCD, value, label, var) |>
  mutate(label = sub("i have", "have", label)) |>
  mutate(label = sub("i am ", "", label)) |>
  mutate(label = stringr::str_to_sentence(label))

### VAS ----

adeq5d_vas <- adeq5d |>
  filter(PARAMCD == vasparam)

# Create baseline variable.
adeq5d_vas_bl <- adeq5d_vas |>
  filter(AVISIT == visits[1]) |>
  mutate(BASE_VAS = AVAL) |>
  select(USUBJID, BASE_VAS)

# Compute change from baseline to double check.
adeq5d_vas <- adeq5d_vas |>
  left_join(adeq5d_vas_bl, by = "USUBJID") |>
  mutate(CHG_VAS = AVAL - BASE_VAS)
stopifnot(all.equal(
  adeq5d_vas$CHG_VAS,
  ifelse(adeq5d_vas$AVISIT == visits[1], 0, adeq5d_vas$CHG),
  check.attributes = FALSE
))

### Dimensions ----

adeq5d_dims <- adeq5d |>
  filter(PARAMCD %in% paramcds) |>
  mutate(PARAMCD = factor(PARAMCD, levels = paramcds)) |>
  mutate(PARAM = forcats::fct_recode(PARAMCD, !!!paramcds)) |>
  mutate(AVALC = factor(AVALC, levels = vals$AVALC))

# Check which subjects have baseline, if not make CHGCAT1 missing.
adeq5d_dims_bl <- adeq5d_dims |>
  filter(AVISIT == visits[1]) |>
  mutate(has_bl = "Y", BASE_AVALC = AVALC) |>
  select(USUBJID, PARAMCD, has_bl, BASE_AVALC)

adeq5d_dims <- adeq5d_dims |>
  left_join(adeq5d_dims_bl, by = c("USUBJID", "PARAMCD")) |>
  mutate(
    CHGCAT1 = case_when(
      is.na(has_bl) ~ "",
      TRUE ~ CHGCAT1
    )
  ) |>
  mutate(
    CHGCAT1 = factor(CHGCAT1, levels = c("Improved", "No change", "Worsened"))
  )

## Analysis ----

### VAS ----

ana_vas <- adeq5d_vas |>
  inner_join(adsl, by = "USUBJID") |>
  mutate(
    VASLABEL = "EQ-5D VAS",
    AVISIT_superscript = paste0(AVISIT, "~[super sup1]")
  )

### Dimensions ----

ana_dims <- adeq5d_dims |>
  inner_join(adsl, by = "USUBJID") |>
  mutate(
    EQ5DLABEL = "EQ-5D dimensions"
  )

# Layout ----

## VAS ----

lyt_vas <- basic_table(show_colcounts = TRUE) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  split_rows_by("VASLABEL") |>
  split_rows_by(
    "AVISIT",
    split_fun = keep_split_levels(only = visits[1])
  ) |>
  analyze_vars(
    "BASE_VAS",
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
  insert_blank_line() |>
  split_rows_by(
    "AVISIT",
    labels_var = "AVISIT_superscript",
    split_fun = drop_and_remove_levels(excl = visits[1]),
    section_div = " ",
    indent = 1
  ) |>
  analyze_vars(
    "AVAL",
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
    ),
    table_names = "aval_desc",
    show_labels = "hidden"
  ) |>
  insert_blank_line() |>
  analyze_vars(
    "CHG",
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
    ),
    table_names = "chg_desc",
    show_labels = "visible",
    var_labels = paste("Change from", tolower(visits[1]))
  ) |>
  insert_blank_line() |>
  analyze(
    "CHG",
    afun = a_summarize_ancova_j,
    na_str = default_na_str(),
    show_labels = "hidden",
    indent = 1,
    extra_args = list(
      variables = list(
        arm = trtvar,
        covariates = strata
      ),
      conf_level = conflvl,
      weights_emmeans = "equal",
      ref_path = ref_path,
      .stats = c("lsmean_se", "lsmean_diffci", "pval"),
      .indent_mods = c(lsmean_se = 1, lsmean_diffci = 1, pval = 2),
      .labels = c(
        lsmean_se = "LS Means (SE)~[super sup_LS]",
        lsmean_diffci = paste0(
          "LS mean difference (",
          tern::f_conf_level(conflvl),
          ")~[super sup_LSdiff]"
        ),
        pval = "p-value~[super sup_p]"
      ),
      .formats = c(
        lsmean_se = formats$lsmean_se,
        lsmean_diffci = formats$lsmean_diffci,
        pval = formats$pval
      )
    )
  )

## Dimensions ----

lyt_dims <- basic_table(show_colcounts = TRUE) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  split_rows_by("EQ5DLABEL") |>
  split_rows_by("PARAMCD", labels_var = "PARAM", section_div = " ") |>
  split_rows_by("AVISIT", split_fun = drop_split_levels, section_div = " ") |>
  analyze(
    "AVALC",
    show_labels = "hidden",
    afun = a_freq_j,
    section_div = " ",
    extra_args = list(
      label_map = xlabel_map,
      riskdiff = FALSE,
      .stats = c("n_df", "count_unique_fraction"),
      denom = "n_df"
    )
  ) |>
  analyze(
    "CHGCAT1",
    var_labels = paste0("Change from ", tolower(visits[1]), "~[super sup1]"),
    afun = a_freq_j_with_exclude,
    extra_args = list(
      exclude_levels = list(AVISIT = visits[1]),
      riskdiff = FALSE,
      .stats = c("n_df", "count_unique_fraction"),
      denom = "n_df"
    )
  ) |>
  analyze(
    "CHGCAT1",
    table_names = "CHG_cmh_rms",
    afun = a_cmhrms_j_with_exclude,
    show_labels = "hidden",
    extra_args = list(
      exclude_levels = list(AVISIT = visits[1]),
      ref_path = ref_path,
      variables = list(
        arm = trtvar,
        strata = strata
      ),
      .formats = c(pval = formats$pval),
      .labels = c(pval = "p-value~[super sup_pCMH]")
    ),
    indent_mod = 2
  )

# Output ----

result_vas <- build_table(lyt_vas, df = ana_vas, alt_counts_df = adsl)
result_dims <- build_table(lyt_dims, df = ana_dims, alt_counts_df = adsl) |>
  prune_table(prune_func = keep_rows(keep_non_null_rows))
result <- rbind(result_vas, result_dims)

# Add title and main footnotes.
result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
