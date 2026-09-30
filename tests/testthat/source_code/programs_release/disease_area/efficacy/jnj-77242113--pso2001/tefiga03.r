###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefiga03.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefiga03: IGA Response Analysis
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     adsl.sas7bdat, adigrspi.sas7bdat
## Output:                    tefiga03.rtf
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

tblid <- "TEFIGA03"
fileid <- write_path(opath, tblid)
popfl <- "FASFL"
trtvar <- "TRT01P"

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

timepoints <- c("Week 1", "Week 2", "Week 4", "Week 8", "Week 12", "Week 16")

################################################################################
# Process data:
################################################################################

# Read in SAS dataset and convert to R dataframe.
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
  select(USUBJID, !!rlang::sym(popfl), !!rlang::sym(trtvar), colspan_trt)

# Read in SAS dataset and convert to R dataframe.
adrsp <- haven::read_sas(read_path(a_in, "adigrspi.sas7bdat")) |>
  filter(
    !!rlang::sym(popfl) == "Y" &
      PARAMCD %in% c("IGMLDP", "IGMINP", "IGCLRP") &
      AVISIT %in% timepoints
  ) |>
  mutate(
    AVISIT = forcats::fct_reorder(factor(AVISIT), AVISITN),
    response = factor(
      case_when(
        PARAMCD == "IGCLRP" & AVALC == "Y" ~ "IGA of cleared (0)",
        PARAMCD == "IGMINP" &
          AVALC == "Y" ~
          "IGA of cleared (0) or minimal (1)",
        PARAMCD == "IGMLDP" &
          AVALC == "Y" ~
          "IGA of cleared (0), minimal (1), or mild (2)",
        TRUE ~ "Nonresponder"
      ),
      levels = c(
        "Nonresponder",
        "IGA of cleared (0)",
        "IGA of cleared (0) or minimal (1)",
        "IGA of cleared (0), minimal (1), or mild (2)"
      )
    )
  ) |>
  select(USUBJID, AVISIT, response)

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

lyt <- basic_table(show_colcounts = TRUE, colcount_format = "N=xx") |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  split_rows_by(
    "AVISIT",
    split_fun = keep_split_levels(only = timepoints),
    section_div = " "
  ) |>
  analyze(
    vars = "response",
    afun = a_freq_j,
    extra_args = list(
      denom = "n_df",
      .stats = c("n_df", "count_unique_fraction"),
      val = c(
        "IGA of cleared (0)",
        "IGA of cleared (0) or minimal (1)",
        "IGA of cleared (0), minimal (1), or mild (2)"
      )
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
