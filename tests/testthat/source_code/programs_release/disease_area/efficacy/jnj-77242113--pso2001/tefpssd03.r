###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefpssd03.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefpssd03:
##                            [Primary/Secondary Endpoint Analysis] ([Primary Estimand], [Composite Strategy]):
##                            Subjects Achieving PSSD Individual Sign Scale Score of 0 at [Time Point]
##                            Among Subjects With Baseline PSSD Individual Scale Score ≥1;
##                            Full Analysis Set (Study psoriasis)
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     adsl.sas7bdat, adparspi.sas7bdat
## Output:                    tefpssd03.rtf
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

################################################################################
# Define script level parameters:
################################################################################

tblid <- "TEFPSSD03"
fileid <- write_path(opath, tblid)
popfl <- "FASFL"
trtvar <- "TRT01P"

tab_titles <- get_titles_internal(tblid)

timepoint <- "Week 16"
stratvar <- "STRATWTG"


### PSSD parameters setup
pssdparamcds <- tribble(
  ~PARAMCD                                                                  ,
  ~PARAMN                                                                   ,
  ~PARCAT1                                                                  ,
  ~PARAM                                                                    ,
  "PSSGN0P"                                                                 ,
                                                                          1 ,
  "SIGN SCORE"                                                              ,
  "PSSD Sign Score of 0 (PE)"                                               ,
  "PSDRY0P"                                                                 ,
                                                                          2 ,
  "SIGN SCORE"                                                              ,
  "PSSD Dryness Score of 0 (PE)"                                            ,
  "PSCRK0P"                                                                 ,
                                                                          3 ,
  "SIGN SCORE"                                                              ,
  "PSSD Cracking Score of 0 (PE)"                                           ,
  "PSSCL0P"                                                                 ,
                                                                          4 ,
  "SIGN SCORE"                                                              ,
  "PSSD Scaling (Build-Up of Skin) Score of 0 (PE)"                         ,
  "PSSHD0P"                                                                 ,
                                                                          5 ,
  "SIGN SCORE"                                                              ,
  "PSSD Shedding or Flaking Score of 0 (PE)"                                ,
  "PSRED0P"                                                                 ,
                                                                          6 ,
  "SIGN SCORE"                                                              ,
  "PSSD Redness Score of 0 (PE)"                                            ,
  "PSBLD0P"                                                                 ,
                                                                          7 ,
  "SIGN SCORE"                                                              ,
  "PSSD Bleeding Score of 0 (PE)"                                           ,
  "PSSYM0P"                                                                 ,
                                                                         11 ,
  "SYMPTOM SCORE"                                                           ,
  "PSSD Symptom Score of 0 (PE)"                                            ,
  "PSITC0P"                                                                 ,
                                                                         12 ,
  "SYMPTOM SCORE"                                                           ,
  "PSSD Itch Score of 0 (PE)"                                               ,
  "PSTGT0P"                                                                 ,
                                                                         13 ,
  "SYMPTOM SCORE"                                                           ,
  "PSSD Skin Tightness Score of 0 (PE)"                                     ,
  "PSBRN0P"                                                                 ,
                                                                         14 ,
  "SYMPTOM SCORE"                                                           ,
  "PSSD Burning Score of 0 (PE)"                                            ,
  "PSSTG0P"                                                                 ,
                                                                         15 ,
  "SYMPTOM SCORE"                                                           ,
  "PSSD Stinging Score of 0 (PE)"                                           ,
  "PSPAI0P"                                                                 ,
                                                                         16 ,
  "SYMPTOM SCORE"                                                           ,
  "PSSD Pain from Your Psoriasis Lesions Score of 0 (PE)"                   ,
  "PSSGNIP"                                                                 ,
                                                                         21 ,
  "SIGN SCORE IMPROVE"                                                      ,
  "PSSD Sign Clinically Meaningful Improvement (PE)"                        ,
  "PSDRYIP"                                                                 ,
                                                                         22 ,
  "SIGN SCORE IMPROVE"                                                      ,
  "PSSD Dryness Clinically Meaningful Improvement (PE)"                     ,
  "PSCRKIP"                                                                 ,
                                                                         23 ,
  "SIGN SCORE IMPROVE"                                                      ,
  "PSSD Cracking Clinically Meaningful Improvement (PE)"                    ,
  "PSSCLIP"                                                                 ,
                                                                         24 ,
  "SIGN SCORE IMPROVE"                                                      ,
  "PSSD Scaling (Build-Up of Skin) Clinically Meaningful Improvement (PE)"  ,
  "PSSHDIP"                                                                 ,
                                                                         25 ,
  "SIGN SCORE IMPROVE"                                                      ,
  "PSSD Shedding or Flaking Clinically Meaningful Improvement (PE)"         ,
  "PSREDIP"                                                                 ,
                                                                         26 ,
  "SIGN SCORE IMPROVE"                                                      ,
  "PSSD Redness Clinically Meaningful Improvement (PE)"                     ,
  "PSBLDIP"                                                                 ,
                                                                         27 ,
  "SIGN SCORE IMPROVE"                                                      ,
  "PSSD Bleeding Clinically Meaningful Improvement (PE)"                    ,
  "PSSYMIP"                                                                 ,
                                                                         31 ,
  "SYMPTOM SCORE IMPROVE"                                                   ,
  "PSSD Symptom Clinically Meaningful Improvement (PE)"                     ,
  "PSITCIP"                                                                 ,
                                                                         32 ,
  "SYMPTOM SCORE IMPROVE"                                                   ,
  "PSSD Itch Clinically Meaningful Improvement (PE)"                        ,
  "PSTGTIP"                                                                 ,
                                                                         33 ,
  "SYMPTOM SCORE IMPROVE"                                                   ,
  "PSSD Skin Tightness Clinically Meaningful Improvement (PE)"              ,
  "PSBRNIP"                                                                 ,
                                                                         34 ,
  "SYMPTOM SCORE IMPROVE"                                                   ,
  "PSSD Burning Clinically Meaningful Improvement (PE)"                     ,
  "PSSTGIP"                                                                 ,
                                                                         35 ,
  "SYMPTOM SCORE IMPROVE"                                                   ,
  "PSSD Stinging Clinically Meaningful Improvement (PE)"                    ,
  "PSPAIIP"                                                                 ,
                                                                         36 ,
  "SYMPTOM SCORE IMPROVE"                                                   ,
  "PSSD Pain from Psoriasis Lesions Clinically Meaningful Improvement (PE)" ,
) %>%
  arrange(PARAMN)


# specify in the proper order you want on the table
paramcds <- pssdparamcds %>%
  filter(PARCAT1 == "SYMPTOM SCORE") %>%
  filter(!(PARAMCD %in% c("PSSGN0P", "PSSYM0P"))) %>%
  arrange(PARAMN) %>%
  pull(PARAMCD)


# mapping to show proper response label - use in a_freq_j as label_map
mapdf <- pssdparamcds %>%
  mutate(
    label = tolower(stringr::str_replace(PARAM, stringr::fixed("(PE)"), "at"))
  ) %>%
  mutate(
    label = stringr::str_replace(label, stringr::fixed("pssd"), "PSSD")
  ) %>%
  mutate(label = paste("Subjects achieving", label, timepoint)) %>%
  mutate(value = "Y") %>%
  filter(PARAMCD %in% paramcds)


# proper labels for filtering subjects on baseline condition
bltext <- pssdparamcds %>%
  mutate(
    label = tolower(stringr::str_replace(
      PARAM,
      stringr::fixed("of 0 (PE)"),
      ">= 1"
    ))
  ) %>%
  mutate(
    label = tolower(stringr::str_replace(
      label,
      stringr::fixed("clinically meaningful improvement (pe)"),
      "score >= 1"
    ))
  ) %>%
  mutate(
    label = stringr::str_replace(label, stringr::fixed("pssd"), "PSSD")
  ) %>%
  mutate(label = paste("Subjects with baseline", label)) %>%
  rename(BLTEXT = label) %>%
  filter(PARAMCD %in% paramcds)

# prepare for proper baseline condition, might depend on paramcd
# for signs and symptoms : >= 1 - BLGE1FL

blcrit <- pssdparamcds %>%
  mutate(BLSELVAR = "BLGE1FL") %>%
  filter(PARAMCD %in% paramcds)

# use case_when if different vars are needed for different paramcds, eg BLGE4FL for some or BLGE3FL for others
# ensure to update bltext dataframe accordingly
# ensure to look at mapdf as well
# eg
blcrit <- blcrit %>%
  mutate(
    BLSELVAR = case_when(
      PARAMCD == "PSSYMIP" ~ "BLGE40FL",
      PARAMCD == "PSITCIP" ~ "BLGE4FL",
      PARAMCD == "PSTGTIP" ~ "BLGE4FL",
      PARAMCD == "PSBRNIP" ~ "BLGE4FL",
      PARAMCD == "PSSTGIP" ~ "BLGE3FL",
      PARAMCD == "PSPAIIP" ~ "BLGE4FL",
      TRUE ~ BLSELVAR
    )
  )
bltext <- bltext %>%
  mutate(
    BLTEXT = case_when(
      PARAMCD == "PSSYMIP" ~ "Subjects with baseline symptom score  >= 40",
      PARAMCD == "PSITCIP" ~ "Subjects with baseline itch score >= 4",
      PARAMCD == "PSTGTIP" ~ "Subjects with baseline skin score  >= 4",
      PARAMCD == "PSBRNIP" ~ "Subjects with baseline burning score  >= 4",
      PARAMCD == "PSSTGIP" ~ "Subjects with baseline stinging score  >= 3",
      PARAMCD == "PSPAIIP" ~ "Subjects with baseline pain score  >= 4",
      TRUE ~ BLTEXT
    )
  )
mapdf <- mapdf %>%
  mutate(
    label = case_when(
      PARAMCD == "PSSYMIP" ~ "Subjects achieving >= 40 points improvement in symptom score",
      PARAMCD == "PSITCIP" ~ "Subjects achieving >= 4 points improvement in itch score",
      PARAMCD == "PSTGTIP" ~ "Subjects achieving >= 4 points improvement in skin score",
      PARAMCD == "PSBRNIP" ~ "Subjects achieving >= 4 points improvement in burning score",
      PARAMCD == "PSSTGIP" ~ "Subjects achieving >= 3 points improvement in stinging score",
      PARAMCD == "PSPAIIP" ~ "Subjects achieving >= 4 points improvement in pain score",
      TRUE ~ label
    )
  )

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
  select(
    USUBJID,
    !!rlang::sym(popfl),
    !!rlang::sym(trtvar),
    colspan_trt,
    !!rlang::sym(stratvar)
  )

# Read in SAS dataset and convert to R dataframe.
adresp <- haven::read_sas(read_path(a_in, "adpsrspi.sas7bdat")) |>
  filter(
    !!rlang::sym(popfl) == "Y" & AVISIT == timepoint & PARAMCD %in% paramcds
  ) |>
  select(
    USUBJID,
    AVISIT,
    AVISITN,
    PARAMCD,
    PARAM,
    starts_with("BLGE"),
    AVALC
  ) |>
  inner_join(bltext) |>
  inner_join(blcrit) %>%
  #### setup to filter on baseline selection criteria
  rowwise() %>%
  mutate(BLSEL = get(BLSELVAR) == "Y") %>%
  ungroup() %>%
  #### actual filter on baseline selection criteria
  filter(BLSEL == TRUE)

adresp <- adresp %>%
  mutate(AVISIT = forcats::fct_reorder(factor(AVISIT), AVISITN)) |>
  mutate(PARAMCD = factor(PARAMCD, levels = paramcds)) |>
  mutate(
    response = case_when(
      AVALC == "Y" ~ TRUE,
      TRUE ~ FALSE
    )
  )

adresp <- inner_join(x = adsl, y = adresp, by = "USUBJID")

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

ref_path <- c("colspan_trt", " ", trtvar, "PLACEBO")

lyt <- basic_table(show_colcounts = TRUE, colcount_format = "N=xx") |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  split_rows_by("PARAMCD", labels_var = "BLTEXT", section_div = " ") |>
  summarize_row_groups(
    "PARAMCD",
    cfun = a_freq_j,
    extra_args = list(.stats = c("n_df"))
  ) |>
  analyze(
    vars = "AVALC",
    afun = a_freq_j,
    na_str = default_na_str(),
    table_names = "est_prop",
    show_labels = "hidden",
    extra_args = list(
      val = "Y",
      .stats = c("count_unique_fraction"),
      denom = "n_df",
      label_map = mapdf
    )
  ) |>
  # optional - remove/disable if not needed
  analyze(
    vars = "response",
    afun = a_proportion_diff_j,
    na_str = default_na_str(),
    table_names = "est_prop_diff",
    show_labels = "hidden",
    extra_args = list(
      .stats = c("diff_est_ci"),
      .labels = c(diff_est_ci = "% Difference (95% CI)"),
      .formats = c(diff_est_ci = jjcsformat_xx("xx.x (xx.x, xx.x)")),
      .indent_mods = 1,
      method = "cmh",
      variables = list(strata = stratvar),
      ref_path = ref_path
    )
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

result <- build_table(lyt, adresp, alt_counts_df = adsl)

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)
################################################################################
# Convert to tbl file and output table:
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
