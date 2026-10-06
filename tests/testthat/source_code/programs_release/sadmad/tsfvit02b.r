###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsfvit02b.r
## R version:                 4.5.1
## junco version:             0.1.6
## Short Description:         Subjects With On-treatment Clinically Important
##                            Vital Signs Based on FDA Toxicity Grading Scale
##                            for Healthy Adult and Adolescent Volunteers -
##                            [SAD/MAD][Part 1]
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, advs
## Output:                    tsfvit02b.rtf
## Remarks:                   Should only be used when abnormalities are based
##                            upon vital science toxicity grading
##
## Modification History:
##  Rev #:
##  Modified By:
##  Reporting Effort:
##  Date:
##  Description:
################################################################################

is_template_pgm <- TRUE #### !!!! set to FALSE for actual study!!
test_stop <- FALSE #### !!!! set to FALSE for actual study!!

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
tblid <- "TSFVIT02b"
fileid <- write_path(opath, tblid)
titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
popfl <- "SAFFL"
trtvar <- "TRT01A"
ctrl_grp <- "Pooled Placebo"
combined_colspan_trt <- TRUE

# add 'Grade >= 2' row to the output, set to FALSE if not required
add_ge_2_row <- TRUE

## if the option TRTEMFL needs to be added to the TLF
trtemfl <- TRUE

## ANL flag variable for worst on-treatment low grade
anl_low_fl <- "ANL04FL"

## ANL flag variable for worst on-treatment high grade
anl_high_fl <- "ANL05FL"

studypart_var <- "PARTC"
studyprt_val <- "PART 1"

actarm_str <- "SAD"

active_study_lbl <- "Active Study Agent"


################################################################################
# Initial processing of data + check if table is valid for trial:
################################################################################
advs_complete <- haven::read_sas(read_path(a_in, "advs.sas7bdat")) %>%
  df_na()

################################################################################
# Check if table should be produced for study
################################################################################

check_fda_cols <- bds_check_abn_nci_fda_daids(
  bds_df = advs_complete,
  tblid = tblid,
  template = is_template_pgm,
  test_stop = test_stop
)

if (check_fda_cols != "OK") {
  stop("Inappropriate table for current study", call. = FALSE)
}

vstoxgrade_file <- read_path(dpspath, "vstoxgrade.xlsx")
vstoxgrade_sheets <- readxl::excel_sheets(path = vstoxgrade_file)

# Check if FDA toxicity grading sheet is available
if (!"VSTOXADULTFDA_2007" %in% vstoxgrade_sheets) {
  stop("FDA toxicity grading sheet missing")
}

vstoxgrade_defs <- readxl::read_excel(vstoxgrade_file, sheet = "VSTOXADULTFDA_2007")

vstoxgrade_defs <- vstoxgrade_defs %>%
  mutate(
    ATOXDSCLH = TOXTERM,
    ATOXGRLH = paste("Grade", TOXGRD)
  ) %>%
  rename(ATOXDIR = INDICATR) %>%
  select(ATOXDSCLH, ATOXGRLH, ATOXDIR) %>%
  distinct()

################################################################################
# Process Data:
################################################################################

adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) %>%
  df_na() %>%
  filter(
    !!rlang::sym(studypart_var) == studyprt_val,
    grepl(actarm_str, ACTARM, ignore.case = TRUE),
    .data[[popfl]] == "Y"
  ) %>%
  select(USUBJID, all_of(c(popfl, trtvar))) %>%
  mutate(
    !!trtvar := if_else(
      .data[[trtvar]] == "Placebo",
      ctrl_grp,
      .data[[trtvar]]
    ),
    !!trtvar := factor(.data[[trtvar]])
  )

adsl$colspan_trt <- factor(
  if_else(adsl[[trtvar]] == ctrl_grp, " ", active_study_lbl),
  levels = c(active_study_lbl, " ")
)

colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = ctrl_grp,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = active_study_lbl,
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

if (combined_colspan_trt) {
  add_combo <- add_combo_facet(
    "Combined",
    label = "Combined",
    levels = setdiff(unique(adsl[[trtvar]]), ctrl_grp)
  )

  rm_combo_from_placebo <- cond_rm_facets(
    facets = "Combined",
    ancestor_pos = NA,
    value = " ",
    split = "colspan_trt"
  )

  mysplit <- make_split_fun(post = list(add_combo, rm_combo_from_placebo))
}

advs00 <- advs_complete %>%
  select(
    USUBJID,
    starts_with("ATOX"),
    starts_with("ANL"),
    ONTRTFL,
    TRTEMFL,
    AVISIT,
    PARAMCD
  ) %>%
  inner_join(adsl) %>%
  mutate(
    ATOXGRL = as.character(ATOXGRL),
    ATOXGRH = as.character(ATOXGRH)
  ) %>%
  relocate(
    .,
    USUBJID,
    all_of(anl_low_fl),
    all_of(anl_high_fl),
    ONTRTFL,
    TRTEMFL,
    ATOXGRL,
    ATOXGRH,
    ATOXDSCL,
    ATOXDSCH
  )

# Filtered ADVS
filtered_advs <- advs00 %>%
  filter(ONTRTFL == 'Y')

### low grades : ATOXDSCL ATOXGRL ANL04FL
### Note on Worst On-treatment
### note: by filter ANL04FL/ANL05FL, this table is restricted to On-treatment
### values, per definition of ANL04FL/ANL05FL therefor, no need to add ONTRTFL
### in filter if derivation of ANL04FL/ANL05FL is not restricted to ONTRTFL
### records, adding ONTRTFL here will not give the correct answer either
### as mixing worst with other period is not giving the proper selection !!!

filtered_advs_low <- filtered_advs %>%
  filter(.data[[anl_low_fl]] == "Y" & !is.na(ATOXDSCL) & !is.na(ATOXGRL)) %>%
  mutate(
    ATOXDSCLH = ATOXDSCL,
    ATOXGRLH = ATOXGRL,
    ATOXDIR = "LOW"
  ) %>%
  select(USUBJID, starts_with("PAR"), starts_with("ATOX"), PARAMCD, AVISIT, TRTEMFL) %>%
  select(-c(ATOXGRL, ATOXGRH, ATOXDSCL, ATOXDSCH))

### high grades: ATOXDSCH ATOXGRH ANL05FL
filtered_advs_high <- filtered_advs %>%
  filter(.data[[anl_high_fl]] == "Y" & !is.na(ATOXDSCH) & !is.na(ATOXGRH)) %>%
  mutate(
    ATOXDSCLH = ATOXDSCH,
    ATOXGRLH = ATOXGRH,
    ATOXDIR = "HIGH"
  ) %>%
  select(USUBJID, starts_with("PAR"), starts_with("ATOX"), AVISIT, PARAMCD, TRTEMFL) %>%
  select(-c(ATOXGRL, ATOXGRH, ATOXDSCL, ATOXDSCH))

## combine Low and high into advs_tox
filtered_advs_tox <- bind_rows(
  filtered_advs_low,
  filtered_advs_high
) %>%
  select(-c(ATOXGR, ATOXGRN)) %>%
  inner_join(adsl, by = "USUBJID")

#### DO NOT USE TRTEMFL = Y in filter, as this will remove subjects from both
#### numerator and denominator, instead : set ATOXGRLH to a non-reportable value
#### (ie Grade 0) and keep in dataset
if (trtemfl) {
  filtered_advs_tox <- filtered_advs_tox %>%
    mutate(
      ATOXGRLH = case_when(
        is.na(TRTEMFL) | TRTEMFL != "Y" ~ "0",
        TRUE ~ ATOXGRLH
      )
    )
}

filtered_advs_tox <- filtered_advs_tox %>%
  mutate(
    ATOXGRLH = factor(paste("Grade", ATOXGRLH), levels = paste("Grade", 0:5)),
    ATOXDIR = factor(ATOXDIR, levels = c("LOW", "HIGH"))
  ) %>%
  distinct()

check_non_unique_subject <- filtered_advs_tox %>%
  group_by(USUBJID, PARAMCD, ATOXDSCLH) %>%
  summarize(n_subject = n(), .groups = "drop") %>%
  filter(n_subject > 1)

if (nrow(check_non_unique_subject)) {
  message(
    "Please review your data selection process, subject has multiple records"
  )
}

### add relevant extra vars to vstoxgrade_defs, only restrict to those actually in trial
vstoxgrade_defs <- vstoxgrade_defs %>%
  inner_join(
    .,
    unique(
      filtered_advs_tox %>%
        select(PARAMCD, ATOXDIR, ATOXDSCLH)
    ),
    relationship = "many-to-many"
  )

### Define param_map to be used in layout
param_map <- vstoxgrade_defs %>%
  select(ATOXDIR, PARAMCD, ATOXDSCLH, ATOXGRLH) %>%
  mutate(
    ATOXDIR = factor(ATOXDIR, levels = c("LOW", "HIGH"))
  ) %>%
  arrange(ATOXDSCLH, ATOXDIR) %>%
  mutate(
    PARAMCD = as.character(PARAMCD),
    ATOXDIR = as.character(ATOXDIR),
    ATOXDSCLH = as.character(ATOXDSCLH),
    ATOXGRLH = as.character(ATOXGRLH)
  )

if (add_ge_2_row) {
  filtered_advs_tox_ge2 <- filtered_advs_tox %>%
    filter(ATOXGRLH %in% paste("Grade", 2:5)) %>%
    mutate(
      ATOXGRLH = "Grade >= 2"
    )

  filtered_advs_tox <- bind_rows(
    filtered_advs_tox,
    filtered_advs_tox_ge2
  ) %>%
    mutate(
      ATOXGRLH = factor(
        ATOXGRLH,
        levels = c(paste("Grade", 0:5), "Grade >= 2")
      )
    )

  param_map_ge2 <- param_map %>%
    filter(ATOXGRLH %in% paste("Grade", 2:5)) %>%
    distinct(ATOXDIR, PARAMCD, ATOXDSCLH) %>%
    mutate(
      ATOXGRLH = "Grade >= 2"
    )

  param_map <- bind_rows(
    param_map,
    param_map_ge2
  ) %>%
    mutate(
      ATOXDIR = factor(ATOXDIR, levels = c("LOW", "HIGH"))
    ) %>%
    arrange(ATOXDSCLH, ATOXDIR) %>%
    mutate(ATOXDIR = as.character(ATOXDIR))
}

################################################################################
# Define layout and build table:
###############################################################################
ref_path <- c("colspan_trt", " ", trtvar, ctrl_grp)

extra_args_rr <- list(
  method = "wald",
  denom = "n_df",
  ref_path = ref_path,
  .stats = c("denom", "count_unique_fraction")
)

lyt <- basic_table(show_colcounts = TRUE, colcount_format = "N=xx") %>%
  append_topleft("Vital Sign") %>%
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) %>%
  {
    if (combined_colspan_trt) {
      split_cols_by(., trtvar, split_fun = mysplit)
    } else {
      split_cols_by(., trtvar)
    }
  } %>%
  split_rows_by(
    "ATOXDSCLH",
    label_pos = "topleft",
    child_labels = "visible",
    split_label = " ",
    ### trim_levels_to_map needs to be applied at ALL split_rows_by levels
    split_fun = trim_levels_to_map(param_map),
    section_div = " "
  ) %>%
  append_topleft("  Grade, n (%)") %>%
  analyze(
    "ATOXGRLH",
    a_freq_j,
    extra_args = extra_args_rr,
    show_labels = "hidden",
    indent_mod = 0L
  )

result <- build_table(lyt, filtered_advs_tox, alt_counts_df = adsl, round_type = 'sas')

################################################################################
# Post-Processing: Remove Grade 0 line
################################################################################

remove_grade0 <- function(tr) {
  if (is(tr, "DataRow") & (tr@label == "Grade 0")) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

result <- result %>% prune_table(prune_func = keep_rows(remove_grade0))

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, titles)

################################################################################
# Convert to tbl file and output table
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
