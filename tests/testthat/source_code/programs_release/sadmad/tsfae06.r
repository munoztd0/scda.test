###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsfae06.r
## R version:                 4.5.2
## junco Version:             0.1.6.9000
## Short Description:         Program to create tsfae06: Subjects With Treatment
##                            -emergent Adverse Events of [Special] Interest by
##                            AE Grouping and Preferred Term – [SAD/MAD] [Part 1]
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adae
## Output:                    tsfae06.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
##  Rev #:
##  Modified By:
##  Reporting Effort:
##  Date:
##  Description:
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

################################################################################
# Define script level parameters:
################################################################################

################################################################################
# - Define output ID, file location and titles
# - Define the customized query variable (CQzzNAM) or Standardized MedDRA Query(SMQzzNAM) for the AEs of special interest
# - Choose whether or not you want to present a combined active treatment column (default=TRUE)
# - Define what the control treatment group is for your study (e.g Placebo)
################################################################################

tblid <- "TSFAE06"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
active_study_lbl <- "Active Study Agent"
trtvar <- "TRT01A"
popfl <- "SAFFL"
studypart_var <- "PARTC"
studyprt_val <- "PART 1"
actarm_str <- "SAD"
special_interest_var <- c("CQ07NAM")
combined_colspan_trt <- TRUE
ctrl_grp <- "Pooled Placebo"



################################################################################
# Process Data:
################################################################################

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) %>%
  df_na() %>%
  filter(
    !!rlang::sym(studypart_var) == studyprt_val,
    grepl(actarm_str, ACTARM, ignore.case = TRUE),
    !!sym(popfl) == "Y"
  )

adsl <- adsl %>%
  mutate(
    !!trtvar := as.character(.data[[trtvar]])
  ) %>%
  mutate(
    !!trtvar := ifelse(
      .data[[trtvar]] == "Placebo",
      ctrl_grp,
      .data[[trtvar]]
    ),
    !!trtvar := factor(.data[[trtvar]])
  ) %>%
  select(
    USUBJID,
    all_of(c(popfl, trtvar))
  )

adsl <- adsl %>%
  create_colspan_var(
    non_active_grp = ctrl_grp,
    non_active_grp_span_lbl = " ",
    active_grp_span_lbl = active_study_lbl,
    colspan_var = "colspan_trt",
    trt_var = trtvar
  )

active_levels <- setdiff(
  levels(factor(adsl[[trtvar]])),
  ctrl_grp
)

if (combined_colspan_trt) {
  add_combo <- add_combo_facet(
    "Combined",
    label = "Combined",
    levels = active_levels
  )

  rm_combo_from_placebo <- cond_rm_facets(
    facets = "Combined",
    ancestor_pos = NA,
    value = " ",
    split = "colspan_trt"
  )

  mysplit <- make_split_fun(
    post = list(
      add_combo,
      rm_combo_from_placebo
    )
  )
}

totdf <- tribble(
  ~valname               ,
  ~label                 ,
  ~levelcombo            ,
  ~exargs                ,
  "Total"                ,
  "Total"                ,
  levels(adsl[[trtvar]]) ,
  list()
)

adae <- haven::read_sas(envsetup::read_path(a_in, "adae.sas7bdat")) %>%
  mutate(
    AEDECOD = case_when(
      AEDECOD == "" ~ paste0("Uncoded: ", AETERM),
      .default = AEDECOD
    )
  ) %>%
  df_na() %>%
  filter(TRTEMFL == "Y", !!sym(popfl) == "Y", if_any(all_of(special_interest_var), ~ !is.na(.))) %>%
  select(
    USUBJID,
    TRTEMFL,
    AEDECOD,
    all_of(special_interest_var)
  ) %>%
  tidyr::pivot_longer(
    cols = all_of(special_interest_var),
    names_to = "special_interest_var",
    values_to = "special_interest_val"
  ) %>%
  select(
    -any_of(trtvar)
  )


# join data together
ae <- adae %>% inner_join(., adsl, by = c("USUBJID"))

colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = ctrl_grp,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = active_study_lbl,
  colspan_var = "colspan_trt",
  trt_var = trtvar
)


################################################################################
# Define layout and build table:
################################################################################
extra_args <- list(
  denom = "n_altdf",
  .stats = c("count_unique_fraction")
)

lyt <- basic_table(
  top_level_section_div = " ",
  show_colcounts = TRUE,
  colcount_format = "N=xx"
) %>%
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  )

if (combined_colspan_trt == TRUE) {
  lyt <- lyt %>%
    split_cols_by(trtvar, split_fun = mysplit)
} else {
  lyt <- lyt %>%
    split_cols_by(trtvar)
}

lyt <- lyt %>%
  split_cols_by(
    trtvar,
    split_fun = add_combo_levels(totdf, keep_levels = "Total"),
    nested = FALSE
  )

lyt <- lyt %>%
  split_rows_by(
    "special_interest_val",
    split_label = "AE Grouping",
    split_fun = trim_levels_in_group("AEDECOD"),
    label_pos = "topleft",
    indent_mod = 0,
    section_div = c(" ")
  ) %>%
  summarize_row_groups(
    "special_interest_val",
    cfun = a_freq_j,
    extra_args = append(extra_args, NULL)
  ) %>%
  analyze(
    "AEDECOD",
    afun = a_freq_j,
    extra_args = extra_args,
    show_labels = "hidden",
  ) %>%
  append_topleft(c(" ", "  Preferred Term, n (%)"))

result <- build_table(lyt, ae, alt_counts_df = adsl, round_type = "sas")

# If there is no data display "No data to display" text
if (nrow(ae) == 0) {
  result <- safe_prune_table(result, empty_msg = "No data to report")
}
################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

################################################################################
# Convert to tbl file and output table
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
