###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefiga05.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefiga05: IGA Response Analysis by Region
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     adsl.sas7bdat, adigrspi.sas7bdat
## Output:                    tefiga05.rtf
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

tblid <- "TEFIGA05"
fileid <- write_path(opath, tblid)
popfl <- "FASFL"
trtvar <- "TRT01P"

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

timepoint <- "Week 16"
resp_paramcd <- "IGMINP"

################################################################################
# Process data:
################################################################################

# Read in SAS dataset and convert to R dataframe.
adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  mutate(
    !!popfl := factor(!!rlang::sym(popfl)),
    GEORGR1 = factor(GEOGR1),
    COUNTRY = factor(COUNTRY),
    SITEID = factor(SITEID),
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
  filter(!!rlang::sym(popfl) == "Y") |>
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
    GEOGR1,
    COUNTRY,
    SITEID
  )

adrsp <- haven::read_sas(read_path(a_in, "adigrspi.sas7bdat")) |>
  filter(
    !!rlang::sym(popfl) == "Y" & PARAMCD == resp_paramcd & AVISIT == timepoint
  ) |>
  mutate(AVALC = factor(AVALC)) |>
  select(USUBJID, PARAMCD, AVALC)

adrsp <- inner_join(x = adsl, y = adrsp, by = "USUBJID")

################################################################################
# Define layout and build table:
################################################################################

# Map each treatment group to appropriate columns spanning header.
colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = "PLACEBO",
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

# Map each country to the appropriate region.
region_map <- tribble(
  ~GEOGR1         ,
  ~COUNTRY        ,
  "Asia Pacific"  ,
  "JPN"           ,
  "Asia Pacific"  ,
  "KOR"           ,
  "Asia Pacific"  ,
  "TWN"           ,
  "Europe"        ,
  "FRA"           ,
  "Europe"        ,
  "DEU"           ,
  "Europe"        ,
  "POL"           ,
  "Europe"        ,
  "ESP"           ,
  "Europe"        ,
  "GBR"           ,
  "North America" ,
  "CAN"           ,
  "North America" ,
  "USA"
)

lyt <- basic_table(show_colcounts = TRUE, colcount_format = "N=xx") |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  split_rows_by(
    var = "GEOGR1",
    section_div = " ",
    split_fun = trim_levels_to_map(map = region_map)
  ) |>
  summarize_row_groups(
    var = "GEOGR1",
    cfun = response_by_var,
    extra_args = list(
      resp_var = "AVALC",
      .format = jjcsformat_fraction_count_denom
    )
  ) |>
  split_rows_by(
    var = "COUNTRY",
    split_fun = trim_levels_in_group(innervar = "SITEID")
  ) |>
  summarize_row_groups(
    var = "COUNTRY",
    cfun = response_by_var,
    extra_args = list(
      resp_var = "AVALC",
      .format = jjcsformat_fraction_count_denom
    )
  ) |>
  analyze(
    vars = "SITEID",
    afun = response_by_var,
    extra_args = list(
      resp_var = "AVALC",
      .format = jjcsformat_fraction_count_denom
    )
  )

result <- build_table(lyt, df = adrsp, alt_counts_df = adsl)

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
