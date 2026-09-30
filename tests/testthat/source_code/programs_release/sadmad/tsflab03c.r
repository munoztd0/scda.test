###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsflab03c.r
## R version:                 4.5.2
## junco version:             0.1.3
## Short Description:         Program to create tsflab03c:
##                            Subjects With ≥1 Laboratory Values Meeting Specified Grades
##                            Based on Worst On-treatment Value Using DAIDS Criteria – [SAD/MAD] [Part 1]
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adlb or adlbc or adlc
## Output:                    tsflab03c.rtf
## Remarks:                   Should only be used when abnormalities are based upon lab toxicity grading
##
## Modification History:
##  Rev #:
##  Modified By:
##  Reporting Effort:
##  Date:
##  Description:
################################################################################

#### BEFORE YOU START USING THIS PROGRAM ENSURE THE FOLLOWING: For your trial you should EITHER use lab toxicity grading (lbtoxgrade file) or Abnormality criteria (markedly abnormal file)
#### As ANL04FL is a single variable to indicate worst, you CANNOT work with both classification methods
#### If your study uses markedly abnormal file for adlb, do not produce this table, instead use TSFLAB02 and TSFLAB04a,TSFLAB04b instead

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


tblid <- "TSFLAB03c"
fileid <- write_path(opath, tblid)
# Load SADMAD-specific titles from a dedicated CSV file.
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

## if the option TRTEMFL needs to be added to the TLF
trtemfl <- TRUE


popfl <- "SAFFL"
trtvar <- "TRT01A"

ctrl_grp <- "Placebo"

combined_colspan_trt <- TRUE

## ANL flag variable for worst on-treatment low grade
anl_low_fl <- "ANL04FL"

## ANL flag variable for worst on-treatment high grade
anl_high_fl <- "ANL05FL"

## For analysis on SI units: use adlb dataset
## For analysis on Conventional units: use adlbc dataset -- shell is in conventional units

ad_domain <- "adlb"


studypart_var <- "PARTC" # need to be changed to STUDYPRT variable
studyprt_val <- "PART 1"

actarm_str <- "SAD"


active_study_lbl <- "Active Study Agent"

################################################################################
# Initial processing of data + check if table is valid for trial:
################################################################################

adlb_complete <- haven::read_sas(read_path(a_in, paste0(tolower(ad_domain), ".sas7bdat"))) %>% df_na()

################################################################################
# Check if table should be produced for study
################################################################################

check_bds_toxicity_criteria <- bds_check_abn_nci_fda_daids(
  bds_df = adlb_complete,
  tblid = tblid,
  template = is_template_pgm,
  test_stop = test_stop
)

if (check_bds_toxicity_criteria != "OK") {
  stop("Inappropriate table for current study", call. = FALSE)
}


################################################################################
# Process lab toxicity file
################################################################################

lbtoxgrade_file <- read_path(dpspath, "lbtoxgrade.xlsx")
lbtoxgrade_sheets <- readxl::excel_sheets(path = lbtoxgrade_file)

#Check if DAIDS toxicity grading sheet is available
if (!"DAIDS21c" %in% lbtoxgrade_sheets) {
  stop("DAIDS toxicity grading sheet missing")
}


lbtoxgrade_defs <- readxl::read_excel(lbtoxgrade_file, sheet = "DAIDS21c")

lbtoxgrade_defs <- unique(
  lbtoxgrade_defs %>%
    select(
      TOXTERM,
      TOXGRD,
      INDICATR
    ) %>%
    mutate(
      ATOXDSCLH = TOXTERM,
      ATOXGRLH = paste("Grade", TOXGRD)
    ) %>%
    rename(
      ATOXDIR = INDICATR
    ) %>%
    select(
      ATOXDSCLH,
      ATOXGRLH,
      ATOXDIR
    )
)

################################################################################
# Process Data:
################################################################################
adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) %>%
  df_na() %>%
  filter(
    !!rlang::sym(studypart_var) == studyprt_val,
    grepl(actarm_str, ACTARM, ignore.case = TRUE),
  ) %>%
  filter(.data[[popfl]] == "Y")


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


colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = ctrl_grp,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = active_study_lbl,
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

adlb <- adlb_complete %>%
  select(
    USUBJID,
    AVISITN,
    AVISIT,
    starts_with("PAR"),
    starts_with("ATOX"),
    starts_with("ANL"),
    ONTRTFL,
    TRTEMFL,
    AVAL,
    APOBLFL,
    ABLFL,
    LVOTFL
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
    AVISIT,
    ATOXGRL,
    ATOXGRH,
    ATOXDSCL,
    ATOXDSCH,
    PARAMCD,
    AVISIT,
    AVAL,
    APOBLFL,
    ABLFL
  )

# Filtered ADLB
filtered_adlb <- adlb %>%
  filter(PARCAT4 == "Graded tests" & ONTRTFL == 'Y')

filtered_adlb <- filtered_adlb %>% mutate(PARCAT1 = tools::toTitleCase(tolower(PARCAT1)))


# Note on Worst On-treatment
# Note: by filter ANL04FL/ANL05FL, this table is restricted to On-treatment values, per definition of ANL04FL/ANL05FL therefor, no need to add ONTRTFL in filter
# if derivation of ANL04FL/ANL05FL is not restricted to ONTRTFL records, adding ONTRTFL here will not give the correct answer either
# as mixing worst with other period is not giving the proper selection !!!

#Low grades : ATOXDSCL ATOXGRL ANL04FL
filtered_adlb_low <- filtered_adlb %>%
  filter(.data[[anl_low_fl]] == "Y" & !is.na(ATOXDSCL) & !is.na(ATOXGRL)) %>%
  mutate(
    ATOXDSCLH = ATOXDSCL,
    ATOXGRLH = ATOXGRL,
    ATOXDIR = "LOW"
  ) %>%
  select(USUBJID, starts_with("PAR"), starts_with("ATOX"), TRTEMFL) %>%
  select(-c(ATOXGRL, ATOXGRH, ATOXDSCL, ATOXDSCH))

#High grades: ATOXDSCH ATOXGRH ANL05FL
filtered_adlb_high <- filtered_adlb %>%
  filter(.data[[anl_high_fl]] == "Y" & !is.na(ATOXDSCH) & !is.na(ATOXGRH)) %>%
  mutate(
    ATOXDSCLH = ATOXDSCH,
    ATOXGRLH = ATOXGRH,
    ATOXDIR = "HIGH"
  ) %>%
  select(USUBJID, starts_with("PAR"), starts_with("ATOX"), TRTEMFL) %>%
  select(-c(ATOXGRL, ATOXGRH, ATOXDSCL, ATOXDSCH))

#Combine Low and high into adlb_tox
filtered_adlb_tox <-
  bind_rows(
    filtered_adlb_low,
    filtered_adlb_high
  ) %>%
  select(-c(ATOXGR, ATOXGRN)) %>%
  inner_join(adsl)


#DO NOT USE TRTEMFL = Y in filter, as this will remove subjects from both numerator and denominator
#instead : set ATOXGRLH to a non-reportable value (ie Grade 0) and keep in dataset
if (trtemfl) {
  filtered_adlb_tox <- filtered_adlb_tox %>%
    mutate(
      ATOXGRLH = case_when(
        is.na(TRTEMFL) | TRTEMFL != "Y" ~ "0",
        TRUE ~ ATOXGRLH
      )
    )
}

#Convert some to factors - lyt(layout) will fail if these are not factors
filtered_adlb_tox <-
  filtered_adlb_tox %>%
  mutate(
    ATOXGRLH = factor(paste("Grade", ATOXGRLH), levels = paste("Grade", 0:5)),
    ATOXDIR = factor(ATOXDIR, levels = c("LOW", "HIGH"))
  )

filtered_adlb_tox <- unique(
  filtered_adlb_tox
)

check_non_unique_subject <- filtered_adlb_tox %>%
  group_by(USUBJID, PARAMCD, ATOXDSCLH) %>%
  summarize(n_subject = n()) %>%
  filter(n_subject > 1)

if (nrow(check_non_unique_subject)) {
  message(
    "Please review your data selection process, subject has multiple records"
  )
}


#Add relevant extra vars to lbtoxgrade_defs, only restrict to those actually in trial
lbtoxgrade_defs <- lbtoxgrade_defs %>%
  inner_join(
    .,
    unique(
      filtered_adlb_tox %>%
        select(PARAMCD, PARAM, ATOXDIR, ATOXDSCLH, PARCAT1)
    ),
    relationship = "many-to-many"
  )

#Define param_map to be used in layout
param_map <- lbtoxgrade_defs %>%
  select(PARAM, PARAMCD, ATOXDIR, ATOXDSCLH, ATOXGRLH, PARCAT1) %>%

  mutate(
    sort_test = case_when(
      grepl("Hypo", ATOXDSCLH, ignore.case = TRUE) ~
        gsub("Hypo", "", ATOXDSCLH, ignore.case = TRUE),

      grepl("Hyper", ATOXDSCLH, ignore.case = TRUE) ~
        gsub("Hyper", "", ATOXDSCLH, ignore.case = TRUE),

      grepl(" low$", ATOXDSCLH, ignore.case = TRUE) ~
        gsub(" low$", "", ATOXDSCLH, ignore.case = TRUE),

      grepl(" high$", ATOXDSCLH, ignore.case = TRUE) ~
        gsub(" high$", "", ATOXDSCLH, ignore.case = TRUE),

      grepl(" decrease$", ATOXDSCLH, ignore.case = TRUE) ~
        gsub(" decrease$", "", ATOXDSCLH, ignore.case = TRUE),

      grepl(" increase$", ATOXDSCLH, ignore.case = TRUE) ~
        gsub(" increase$", "", ATOXDSCLH, ignore.case = TRUE),

      TRUE ~ ATOXDSCLH
    ),
    #Actual sorting -- alphabetic by laboratory test, LOW before HIGH within same term
    sort_order = case_when(
      grepl("Hypo", ATOXDSCLH, ignore.case = TRUE) ~ 1,
      grepl(" low$", ATOXDSCLH, ignore.case = TRUE) ~ 1,
      grepl(" decrease$", ATOXDSCLH, ignore.case = TRUE) ~ 1,

      grepl("Hyper", ATOXDSCLH, ignore.case = TRUE) ~ 2,
      grepl(" high$", ATOXDSCLH, ignore.case = TRUE) ~ 2,
      grepl(" increase$", ATOXDSCLH, ignore.case = TRUE) ~ 2,

      TRUE ~ 1
    )
  ) %>%

  arrange(
    PARCAT1,
    sort_test,
    sort_order
  ) %>%

  select(-sort_test, -sort_order) %>%
  ### !!!! no factors are allowed in this split_fun map definition
  mutate(
    PARAMCD = as.character(PARAMCD),
    PARAM = as.character(PARAM),
    ATOXDIR = as.character(ATOXDIR),
    ATOXDSCLH = as.character(ATOXDSCLH),
    ATOXGRLH = as.character(ATOXGRLH),
    PARCAT1 = as.character(PARCAT1)
  )


################################################################################
# Define layout and build table:
################################################################################

extra_args_freq <- list(
  denom = "n_df",
  .stats = c("denom", "count_unique_fraction")
)

lyt <- basic_table(show_colcounts = TRUE, colcount_format = "N=xx") %>%
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
    "PARCAT1",
    split_fun = drop_split_levels,
    label_pos = "topleft",
    child_labels = "visible",
    split_label = "Laboratory Category"
  ) %>%
  split_rows_by(
    "ATOXDSCLH",
    label_pos = "topleft",
    child_labels = "visible",
    split_label = "Laboratory Test",
    ### trim_levels_to_map needs to be applied at ALL split_rows_by levels
    split_fun = trim_levels_to_map(param_map),
    section_div = " "
  ) %>%
  append_topleft("    Grade, n (%)") %>%
  analyze(
    "ATOXGRLH",
    a_freq_j,
    extra_args = extra_args_freq,
    show_labels = "hidden",
    indent_mod = 0L
  )

result <- build_table(lyt, filtered_adlb_tox, alt_counts_df = adsl, round_type = 'sas')

################################################################################
# Post-Processing:
# - Remove Grade 0 line
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

result <- set_titles(result, tab_titles)

################################################################################
# Convert to tbl file and output table
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
