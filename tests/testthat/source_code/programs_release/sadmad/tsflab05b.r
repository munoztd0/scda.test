###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsflab05b.r
## R version:                 4.5.2
## junco version:             0.1.3
## Short Description:         Program to create tsflab05b:
##                            Subjects With [Last/Any] On-treatment Laboratory
##                            Values [≥ Grade 2] Based on Worst On-treatment Value
##                            Using FDA Toxicity Criteria – [SAD/MAD] [Part 1]
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adlb or adlbc
## Output:                    tsflab05b.rtf
## Remarks:                   Should only be used when abnormalities are based upon lab toxicity file
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
#### If your study uses markedly abnormal file for adlb, do not produce this table -- instead use TSFLAB04a,TSFLAB04b instead

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


tblid <- "TSFLAB05b"
# Load SADMAD-specific titles from a dedicated CSV file.
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
fileid <- write_path(opath, tblid)

# Population flag (default=SAFFL).
popfl <- "SAFFL"

# Actual treatment variable (default=TRT01A).
trtvar <- "TRT01A"

ctrl_grp <- "Pooled Placebo"

combined_colspan_trt <- TRUE

grade_threshold <- "2"

## For toxicity grades, the units should not matter, work with adlb dataset as default

ad_domain <- "adlb"

#Table options:(LAST/ANY)
last_any <- "LAST"

#If ANY, then Subjects with Any on-treatment value >= Level 2 will be presented (ANL04FL/ANL05FL/ONTRTFL will be used here)
#if Last, then Subjects with Last on-treatment value >= Level 2 will be presented (LVOTFL will be used here)

## if the option TRTEMFL needs to be added to the TLF -- ensure the same setting as in tsflabxxxx
trtemfl <- TRUE


studypart_var <- "PARTC" # need to be changed to STUDYPRT variable
studyprt_val <- "PART 1"

actarm_str <- "SAD"


active_study_lbl <- "Active Study Agent"

lbvars <- c("LBTESTCD", "LBTEST", "LBSPEC", "LBMETHOD")

flagvars <- c("ONTRTFL", "TRTEMFL", "LVOTFL")


################################################################################
# Initial processing of data
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

#Check if FDA toxicity grading sheet is available
if (!"LBTOXFDA_2007" %in% lbtoxgrade_sheets) {
  stop("FDA toxicity grading sheet missing")
}

lbtoxgrade_defs <- readxl::read_excel(lbtoxgrade_file, sheet = "LBTOXFDA_2007")

lbtoxgrade_defs <- unique(
  unique(
    lbtoxgrade_defs %>%
      select(all_of(lbvars), TOXTERM, TOXGRD, INDICATR)
  ) %>%
    mutate(
      ATOXDSCLH = TOXTERM,
      ATOXGRLH = paste("Grade", TOXGRD)
    ) %>%
    rename(ATOXDIR = INDICATR) %>%
    select(all_of(lbvars), ATOXDSCLH, ATOXDIR)
)

################################################################################
# Initial processing of data
################################################################################

#Be aware, there are toxicity terms that are based upon diff tests (example "Neutrophil Count Decreased": NEUT and NEUTSG) !!!!!
#if both tests are included in ADaM dataset, review your derivations carefully
attention_terms <-
  unique(lbtoxgrade_defs %>% select(ATOXDSCLH, all_of(lbvars))) %>%
  group_by(ATOXDSCLH) %>%
  mutate(n = n_distinct(LBTESTCD)) %>%
  filter(n > 1)

ad_toxterms <- bind_rows(
  unique(
    adlb_complete %>%
      select(PARAMCD, ATOXDSCL) %>%
      mutate(TOXTERM = ATOXDSCL, TOXDIR = "LOW")
  ),
  unique(
    adlb_complete %>%
      select(PARAMCD, ATOXDSCH) %>%
      mutate(TOXTERM = ATOXDSCH, TOXDIR = "HIGH")
  )
) %>%
  select(PARAMCD, TOXTERM, TOXDIR) %>%
  filter(!is.na(TOXTERM)) %>%
  arrange(TOXTERM, TOXDIR, PARAMCD)

attention_ad_toxterms <- ad_toxterms %>%
  group_by(TOXTERM) %>%
  mutate(n = n_distinct(PARAMCD)) %>%
  filter(n > 1)

#From here onwards: avoid using PARAMCD, to ensure toxterm is combined
#could also work with param_lookup and lbtoxgrade_defs
#ALERT: do not include PARAMCD here as some toxicity terms are based on more than one PARAMCD: here : NEUT and NEUTSG both have TOXTERM = Neutrophil Count Decreased
toxterms <- unique(
  adlb_complete %>%
    filter(!(is.na(ATOXDSCL) & is.na(ATOXDSCH))) %>%
    ### do not include PARAMCD here!!!!
    select(ATOXDSCL, ATOXDSCH) %>%
    tidyr::pivot_longer(
      .,
      cols = c("ATOXDSCL", "ATOXDSCH"),
      names_to = "VARNAME",
      values_to = "ATOXDSCLH"
    )
) %>%
  filter(!is.na(ATOXDSCLH))

#Convert dataframe into label_map that can be used with the a_freq_j afun function
xlabel_map <- toxterms %>%
  mutate(
    var = "ATOXGRLHx",
    value = "Y",
    label = as.character(ATOXDSCLH)
  ) %>%
  select(ATOXDSCLH, value, label)

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
    all_of(flagvars),
    LBSEQ,
    AVAL,
    AVALC
  ) %>%
  inner_join(adsl) %>%
  mutate(
    ATOXGRL = as.character(ATOXGRL),
    ATOXGRH = as.character(ATOXGRH)
  ) %>%
  relocate(
    .,
    USUBJID,
    ANL04FL,
    ANL05FL,
    ONTRTFL,
    TRTEMFL,
    AVISIT,
    ATOXGRL,
    ATOXGRH,
    ATOXDSCL,
    ATOXDSCH,
    PARAMCD,
    AVISIT,
    AVAL
  )

#Important: previous actions lost the label of variables
filtered_adlb <- var_relabel_list(adlb, var_labels(adlb_complete, fill = T))

filtered_adlb <- filtered_adlb %>% mutate(PARCAT1 = tools::toTitleCase(tolower(PARCAT1)))


#Low grades : ATOXDSCL ATOXGRL ANL04FL
filtered_adlb_low <- filtered_adlb %>%
  filter(!is.na(ATOXDSCL) & !is.na(ATOXGRL)) %>%
  mutate(
    ATOXDSCLH = ATOXDSCL,
    ATOXGRLH = ATOXGRL,
    ATOXDIR = "LOW",
    ANL045FL = ANL04FL
  ) %>%
  select(
    USUBJID,
    starts_with("PAR"),
    starts_with("ATOX"),
    all_of(flagvars),
    ANL04FL,
    ANL05FL,
    ANL045FL,
    LBSEQ,
    AVAL,
    AVALC
  ) %>%
  select(-c(ATOXGRL, ATOXGRH, ATOXDSCL, ATOXDSCH))


#High grades: ATOXDSCH ATOXGRH ANL05FL
filtered_adlb_high <- filtered_adlb %>%
  filter(!is.na(ATOXDSCH) & !is.na(ATOXGRH)) %>%
  mutate(
    ATOXDSCLH = ATOXDSCH,
    ATOXGRLH = ATOXGRH,
    ATOXDIR = "HIGH",
    ANL045FL = ANL05FL
  ) %>%
  select(
    USUBJID,
    starts_with("PAR"),
    starts_with("ATOX"),
    all_of(flagvars),
    ANL04FL,
    ANL05FL,
    ANL045FL,
    LBSEQ,
    AVAL,
    AVALC
  ) %>%
  select(-c(ATOXGRL, ATOXGRH, ATOXDSCL, ATOXDSCH))


#Combine Low and high into adlb_tox
filtered_adlb_tox <-
  bind_rows(
    filtered_adlb_low,
    filtered_adlb_high
  ) %>%
  select(-c(ATOXGR, ATOXGRN)) %>%
  mutate(ATOXGRLHN = as.numeric(ATOXGRLH)) %>%
  inner_join(adsl)

filtered_adlb_tox <- unique(
  filtered_adlb_tox
)

filtered_adlb_tox_1 <- filtered_adlb_tox %>%
  ## On treatment
  filter(ONTRTFL == "Y")

# Note on On-treatment

#Note: by filter ANL04FL/ANL05FL, this table is restricted to On-treatment values, per definition of ANL04FL/ANL05FL
#Same for LVOTFL therefor, no need to add ONTRTFL in filter
#If derivation of ANL04FL/ANL05FL/LVOTFL is not restricted to ONTRTFL records, adding ONTRTFL here will not give the correct answer either
#as mixing worst with other period is not giving the proper selection !!!

if (toupper(last_any) == "ANY") {
  filtered_adlb_tox_1 <- filtered_adlb_tox_1 %>%
    ## Optional : Any : ensure to have one record per subject for direction
    filter(ANL04FL == "Y" | ANL05FL == "Y")
}

if (toupper(last_any) == "LAST") {
  filtered_adlb_tox_1 <- filtered_adlb_tox_1 %>%
    ## Optional : last on treatment record only
    filter(LVOTFL == "Y")
}


#DO NOT USE TRTEMFL = Y in filter, as this will remove subjects from both numerator and denominator
#instead : set ATOXGRLH to a non-reportable value (ie Grade 0) and keep in dataset
if (trtemfl) {
  filtered_adlb_tox_1 <- filtered_adlb_tox_1 %>%
    mutate(
      ATOXGRLHN = case_when(
        is.na(TRTEMFL) | TRTEMFL != "Y" ~ 0,
        TRUE ~ ATOXGRLHN
      )
    )
}

#Sorting: alphabetically by base term, LOW before HIGH within same term
atoxdsclh_levels <- filtered_adlb_tox_1 %>%
  distinct(ATOXDSCLH, ATOXDIR) %>%
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
  arrange(sort_test, sort_order) %>%
  pull(ATOXDSCLH) %>%
  as.character()

filtered_adlb_tox_1 <- filtered_adlb_tox_1 %>%
  mutate(
    ATOXGRLHx = case_when(
      ATOXGRLHN >= as.numeric(grade_threshold) ~ "Y",
      TRUE ~ "N"
    )
  ) %>%
  mutate(ATOXGRLHx = factor(ATOXGRLHx, levels = c("Y", "N"))) %>%
  mutate(ATOXDSCLH = factor(ATOXDSCLH, levels = atoxdsclh_levels))


# check uniqueness
check_non_unique_subject <- filtered_adlb_tox_1 %>%
  group_by(USUBJID, ATOXDSCLH) %>%
  mutate(n_subject = n()) %>%
  filter(n_subject > 1)

if (nrow(check_non_unique_subject)) {
  message(
    "Please review your data selection process, subject has multiple records"
  )
}


################################################################################
# Define layout and build table:
################################################################################

extra_args_freq <- list(
  denom = "n_df",
  .stats = c("count_unique_denom_fraction")
)

lyt <- basic_table(show_colcounts = TRUE, colcount_format = "N=xx") %>%
  append_topleft("Laboratory Category") %>%
  ### first columns
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
    label_pos = "topleft",
    child_labels = "visible",
    split_label = " ",
    section_div = " "
  ) %>%
  split_rows_by(
    "ATOXDSCLH",
    label_pos = "topleft",
    split_label = paste0(
      "Laboratory Test >= Grade ",
      grade_threshold,
      ", n/N (%)"
    ),
    child_labels = "hidden",
    split_fun = drop_split_levels
  ) %>%
  analyze(
    "ATOXGRLHx",
    afun = a_freq_j,
    extra_args = append(
      extra_args_freq,
      list(
        val = c("Y"),
        label_map = xlabel_map
      )
    ),
    show_labels = "hidden",
    indent_mod = 0L
  )

result <- build_table(lyt, filtered_adlb_tox_1, alt_counts_df = adsl, round_type = "sas")

#######################################################Y#########################
# Post-Processing:
################################################################################

result <- result %>% prune_table(prune_func = keep_rows(keep_non_null_rows))

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)


################################################################################
# Convert to tbl file and output table
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
