###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsidev02.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create tsidev02: Subjects Who Did Not Meet
##                            Inclusion Criteria or Who Met Exclusion
##                            Criteria by Region and Center
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, ie
## Output:                    tsidev02.rtf
## Remarks:
## R-functions:
## R-function Sample Call:
##
## Modification History:
##  Rev #:
##  Modified By:
##  Reporting Effort:
##  Date:
##  Description:
###############################################################################

################################################################################
# Prep environment:
################################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)
library(dplyr)
library(rtables)
library(junco)
library(forcats)

################################################################################
# Define script level parameters:
################################################################################

tblid <- "TSIDEV02"
fileid <- write_path(opath, tblid)
titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
popfl <- "FASFL"
trtvar <- "TRT01A"
ctrl_grp <- "Placebo"


################################################################################
# Initial Read in of adsl dataset
################################################################################

adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) |>
  df_na() |>
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c(
        "Xanomeline Low Dose",
        "Xanomeline High Dose",
        "Placebo"
      )
    )
  )

adsl <- adsl |>
  mutate(SITEID = fct_drop(SITEID))


ie <- haven::read_sas(envsetup::read_path(d_in, "ie.sas7bdat")) |>
  df_na()


################################################################################
# Further script level parameters, after having read in main data
################################################################################

demog_vars <- c("REGION1", "SITEID")
demog_displ_vars <- demog_vars

## make it named vars so that demog_vars[xx] with xx subset of vars still works
names(demog_vars) <- demog_vars
## retrieve labels
demog_labels <- formatters::var_labels(adsl)[demog_vars]

################################################################################
# Process data:
################################################################################

# filter and restrict to population of interest
adsl <- adsl |>
  filter(.data[[popfl]] == "Y") |>
  select(
    USUBJID,
    STUDYID,
    starts_with("TRT01"),
    all_of(popfl),
    all_of(unique(c(demog_vars, demog_displ_vars)))
  )

adsl$colspan_trt <- factor(
  ifelse(adsl[[trtvar]] == ctrl_grp, " ", "Active Study Agent"),
  levels = c("Active Study Agent", " ")
)

# to ensure alphabetical ordering according to COUNTRY variable
adsl$REGION1 <- factor(
  as.character(adsl$REGION1),
  levels = sort(unique(as.character(adsl$REGION1)))
)

ie <- ie |>
  select(STUDYID, USUBJID, IECAT, IETEST, IETESTCD) |>
  mutate(
    IECAT1 = factor(case_when(
      IECAT == "INCLUSION" ~ "Subjects who did not meet inclusion criteria",
      TRUE ~ NA_character_
    )),
    IETSTCD1 = factor(case_when(IECAT == "INCLUSION" ~ IETESTCD, TRUE ~ NA_character_)),
    IECAT2 = factor(case_when(
      IECAT == "EXCLUSION" ~ "Subjects who met exclusion criteria",
      TRUE ~ NA_character_
    )),
    IETSTCD2 = factor(case_when(IECAT == "EXCLUSION" ~ IETESTCD, TRUE ~ NA_character_))
  ) |>
  filter(!is.na(IECAT))

adsl_ie <- adsl |>
  inner_join(ie, by = c("STUDYID", "USUBJID")) |>
  mutate(SITEID = fct_drop(SITEID))

colspan_trt_map <- create_colspan_map(
  adsl_ie,
  non_active_grp = ctrl_grp,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

################################################################################
# Define layout and build table:
################################################################################
extra_args <- list(
  .stats = "count_unique_fraction"
)

lyt <- basic_table(
  top_level_section_div = " ",
  show_colcounts = TRUE,
  colcount_format = "N=xx"
) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  add_overall_col("Total") |>
  # Subjects who did not meet inclusion criteria
  split_rows_by(
    "IECAT1",
    section_div = " "
  ) |>
  summarize_row_groups(
    var = "IECAT1",
    cfun = a_freq_j,
    extra_args = append(
      extra_args,
      list(
        denom = "n_altdf",
        .labels_n = ""
      )
    )
  ) |>
  analyze(
    "IETSTCD1",
    afun = a_freq_j,
    extra_args = append(
      extra_args,
      list(
        denom = "n_altdf",
        .labels_n = ""
      )
    )
  ) |>
  # Subjects who met exclusion criteria
  split_rows_by(
    "IECAT2",
    section_div = " "
  ) |>
  summarize_row_groups(
    var = "IECAT2",
    cfun = a_freq_j,
    extra_args = append(
      extra_args,
      list(
        denom = "n_altdf",
        .labels_n = ""
      )
    )
  ) |>
  analyze(
    "IETSTCD2",
    afun = a_freq_j,
    extra_args = append(
      extra_args,
      list(
        denom = "n_altdf",
        .labels_n = ""
      )
    )
  )

result <- build_table(lyt, adsl_ie, alt_counts_df = adsl, round_type = "sas")

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, titles)

################################################################################
# Convert to tbl file and output table:
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid)
