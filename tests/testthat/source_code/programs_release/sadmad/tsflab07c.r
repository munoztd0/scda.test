###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsflab07c.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create tsflab07c:
##                            Shift in Laboratory Values
##                            From Baseline to Worst Grade
##                            [During Time Period] Based on DAIDS Criteria – [SAD/MAD] [Part 1]
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adlb or adlbc or adlc
## Output:                    tsflab07c.rtf
## Remarks:                   Only produce this table if your trial has used lbtoxgrade file for adlb
##
## Modification History:
##  Rev #:
##  Modified By:
##  Reporting Effort:
##  Date:
##  Description:
################################################################################

#### BEFORE YOU START USING THIS PROGRAM ENSURE THE FOLLOWING: For your trial you should EITHER use lab toxicity grading (lbtoxgrade file)
#### As ANL04FL is a single variable to indicate worst, you CANNOT work with both classification methods
#### Only produce this table if your trial has used lbtoxgrade file for adlb
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


tblid <- "TSFLAB07c"
fileid <- write_path(opath, tblid)
# Load SADMAD-specific titles from a dedicated CSV file.
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

popfl <- "SAFFL"

# Actual treatment variable (default=TRT01A).
trtvar <- "TRT01A"

ctrl_grp <- "Pooled Placebo"

ad_domain <- "adlb"

#Flag variable for worst low toxicity grade (default=ANL04FL).
anl_low_fl <- "ANL04FL"
#Flag variable for worst high toxicity grade (default=ANL05FL).
anl_high_fl <- "ANL05FL"


studypart_var <- "PARTC" # need to be changed to STUDYPRT variable
studyprt_val <- "PART 1"

actarm_str <- "SAD"


flagvars <- c("ONTRTFL", "TRTEMFL", "LVOTFL", "ABLFL")

#Show percentage alongside count (default=FALSE).
show_pct <- FALSE

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
# Read LBTOXGRADE definitions
################################################################################

lbtoxgrade_file <- read_path(dpspath, "lbtoxgrade.xlsx")
lbtoxgrade_sheets <- readxl::excel_sheets(path = lbtoxgrade_file)

#Check if DAIDS toxicity grading sheet is available
if (!"DAIDS21c" %in% lbtoxgrade_sheets) {
  stop("DAIDS toxicity grading sheet missing")
}


lbtoxgrade_defs <- readxl::read_excel(
  lbtoxgrade_file,
  sheet = "DAIDS21c"
)

lbtoxgrade_defs <- unique(
  lbtoxgrade_defs %>%
    select(TOXTERM, TOXGRD, INDICATR) %>%
    mutate(
      ATOXDSCLH = TOXTERM,
      ATOXGRLH = paste("Grade", TOXGRD)
    ) %>%
    rename(ATOXDIR = INDICATR) %>%
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
  filter(.data[[popfl]] == "Y") %>%
  mutate(
    !!trtvar := ifelse(
      as.character(.data[[trtvar]]) == "Placebo",
      ctrl_grp,
      as.character(.data[[trtvar]])
    )
  ) %>%
  {
    trt_levels <- c(
      sort(setdiff(unique(.[[trtvar]]), ctrl_grp)),
      ctrl_grp
    )

    mutate(
      .,
      !!rlang::sym(trtvar) := factor(
        .data[[trtvar]],
        levels = trt_levels
      )
    )
  } %>%
  select(
    USUBJID,
    all_of(c(popfl, trtvar))
  )


adlb <- adlb_complete %>%
  filter(!is.na(ATOXGR)) %>%
  filter(!is.na(BTOXGR)) %>%
  select(
    USUBJID,
    AVISITN,
    AVISIT,
    PARAMCD,
    PARAM,
    PARAMN,
    PARCAT1,
    PARCAT3,
    PARCAT4,
    PARCAT5,
    PARCAT6,
    starts_with("ATOX"),
    starts_with("BTOX"),
    starts_with("ANL"),
    all_of(flagvars),
    LBSEQ,
    AVAL,
    AVALC
  ) %>%
  inner_join(adsl) %>%
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
    ABLFL,
    BTOXGRL,
    BTOXGRH,
    PARAMCD,
    PARCAT3,
    AVISIT,
    AVAL
  )

adlb <- var_relabel_list(adlb, var_labels(adlb_complete, fill = T))

filtered_adlb <- adlb

filtered_adlb <- filtered_adlb %>% mutate(PARCAT1 = tools::toTitleCase(tolower(PARCAT1)))

#Low grades: ATOXDSCL ATOXGRL ANL04FL
filtered_adlb_low <- filtered_adlb %>%
  filter(!is.na(ATOXDSCL) & !is.na(ATOXGRL)) %>%
  mutate(
    ATOXDSCLH = ATOXDSCL,
    ATOXGRLH = ATOXGRL,
    BTOXGRLH = BTOXGRL,
    ATOXDIR = "LOW",
    ANLHLFL = .data[[anl_low_fl]]
  ) %>%
  select(-c(ATOXGRL, ATOXGRH, ATOXDSCL, ATOXDSCH, BTOXGRL, BTOXGRH))

#High grades: ATOXDSCH ATOXGRH ANL05FL
filtered_adlb_high <- filtered_adlb %>%
  filter(!is.na(ATOXDSCH) & !is.na(ATOXGRH)) %>%
  mutate(
    ATOXDSCLH = ATOXDSCH,
    ATOXGRLH = ATOXGRH,
    BTOXGRLH = BTOXGRH,
    ATOXDIR = "HIGH",
    ANLHLFL = .data[[anl_high_fl]]
  ) %>%
  select(-c(ATOXGRL, ATOXGRH, ATOXDSCL, ATOXDSCH, BTOXGRL, BTOXGRH))


#Combine Low and high into adlb_tox
filtered_adlb_tox <-
  bind_rows(
    filtered_adlb_low,
    filtered_adlb_high
  ) %>%
  select(-c(ATOXGR, ATOXGRN)) %>%
  inner_join(adsl)


#Filter for timepoints (either per visit or worst over a time period (overall,period,phase,...))
#here, use worst overall: ANL04FL/ANL05FL --- ANLHLFL
#if worst per period/phase:use ANL07/8 and ANL9/10 --- and ensure to create a combined version ANL078FL, ANL0910FL in the above code
filtered_adlb_tox <-
  filtered_adlb_tox %>%
  filter(ANLHLFL == "Y")


#Add the word Grade to ATOXGRLH,BTOXGRLH, For BTOXGRLH add an extra level (for display purpose)
filtered_adlb_tox <-
  filtered_adlb_tox %>%
  mutate(
    ATOXGRLH = factor(paste("Grade", ATOXGRLH), levels = paste("Grade", 0:4)),
    BTOXGRLH = factor(
      paste("Grade", BTOXGRLH),
      levels = c("N", paste("Grade", 0:4))
    )
  )

# Filter to non-missing PARCAT4 only
filtered_adlb_tox <- filtered_adlb_tox %>%
  filter(PARCAT4 == "Graded tests")


#Sort ATOXDSCLH alphabetically
#Sorting: alphabetically by base term, LOW before HIGH within same term
atoxdsclh_levels <- filtered_adlb_tox %>%
  distinct(ATOXDSCLH, ATOXDIR) %>%
  mutate(
    base_term = sub(", (low|high)$", "", ATOXDSCLH, ignore.case = TRUE),
    ATOXDIR_ord = factor(ATOXDIR, levels = c("LOW", "HIGH"))
  ) %>%
  arrange(base_term, ATOXDIR_ord) %>%
  pull(ATOXDSCLH) %>%
  as.character()


filtered_adlb_tox <- filtered_adlb_tox %>%
  mutate(ATOXDSCLH = factor(as.character(ATOXDSCLH), levels = atoxdsclh_levels))


#Trick for alt_counts_df to work with col splitting
#add BNRIND to adsl, all assign to extra level N (column will be used for N counts)
adsl <- adsl %>%
  mutate(BTOXGRLH = "N") %>%
  mutate(BTOXGRLH = factor(BTOXGRLH, levels = c("N", paste("Grade", 0:4))))

#Add variable for column split header
filtered_adlb_tox$BTOXGRLH_header <- as.factor("Baseline Toxicity Grade")
adsl$BTOXGRLH_header <- as.factor("Baseline Toxicity Grade")

filtered_adlb_tox$BTOXGRLH_header2 <- as.factor(" ") ## first column N should not appear under Baseline column span
adsl$BTOXGRLH_header2 <- as.factor(" ") ## first column N should not appear under Baseline column span

################################################################################
# Create parameter map using LBTOXGRADE and PARCAT1
################################################################################
### add relevant extra vars to lbtoxgrade_defs, only restrict to those actually in trial
lbtoxgrade_defs <- lbtoxgrade_defs %>%
  inner_join(
    .,
    unique(
      filtered_adlb_tox %>%
        select(PARAMCD, PARAM, ATOXDIR, ATOXDSCLH, PARCAT1)
    ),
    relationship = "many-to-many"
  )


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
  #!!!! no factors are allowed in this split_fun map definition
  mutate(
    PARAMCD = as.character(PARAMCD),
    PARAM = as.character(PARAM),
    ATOXDIR = as.character(ATOXDIR),
    ATOXDSCLH = as.character(ATOXDSCLH),
    ATOXGRLH = as.character(ATOXGRLH),
    PARCAT1 = as.character(PARCAT1)
  )


#trick for alt_counts_df to work with col splitting
#add BNRIND to adsl, all assign to extra level N (column will be used for N counts)
adsl <- adsl %>%
  mutate(BTOXGRLH = "N") %>%
  mutate(BTOXGRLH = factor(BTOXGRLH, levels = c("N", paste("Grade", 0:4))))

#add variable for column split header
filtered_adlb_tox$BTOXGRLH_header <- as.factor("Baseline Toxicity Grade")
adsl$BTOXGRLH_header <- as.factor("Baseline Toxicity Grade")

filtered_adlb_tox$BTOXGRLH_header2 <- as.factor(" ") # first column N should not appear under Baseline column span
adsl$BTOXGRLH_header2 <- as.factor(" ") # first column N should not appear under Baseline column span


################################################################################
# Define layout and build table:
################################################################################

lyt <- basic_table(show_colcounts = FALSE) %>%
  ## to ensure N column is not under the Baseline column span header
  split_cols_by("BTOXGRLH_header2") %>%
  split_cols_by("BTOXGRLH", split_fun = keep_split_levels("N")) %>%
  split_cols_by("BTOXGRLH_header", nested = FALSE) %>%
  split_cols_by(
    "BTOXGRLH",
    split_fun = make_split_fun(
      pre = list(rm_levels(excl = "N")),
      post = list(add_overall_facet("TOTAL", "Total"))
    )
  ) %>%
  # Replace the split_rows() + summarize() workflow with a single analyze() call.
  # a_freq_j works because special arguments can be used: denomf = adslx, .stats = count_unique
  # Treatment group counts should come from ADSL rather than the input dataset; therefore, use countsource = altdf.
  analyze(
    vars = trtvar,
    afun = a_freq_j,
    extra_args = list(
      restr_columns = "N",
      .stats = "count_unique",
      countsource = "altdf",
      extrablankline = TRUE
    ),
    indent_mod = -1L
  ) %>%
  #main part of table, restart row-split so nested = FALSE
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
    split_fun = trim_levels_to_map(param_map),
    section_div = " "
  ) %>%
  #--------------------
  split_rows_by(
    trtvar,
    label_pos = "hidden",
    split_label = "Treatment Group",
    section_div = " "
  ) %>%
  summarize_row_groups(
    trtvar,
    cfun = a_freq_j,
    extra_args = list(
      .stats = "n_rowdf",
      restr_columns = c("N")
    ),
  ) %>%
  #add extra level TOTAL using new_levels, rather than earlier technique
  # Create the TOTAL level with new_levels() rather than the legacy method.
  # This allows denominator values to be derived from n_rowdf, making it straightforward to display fractions/percentages.
  # Switch .stats to count_unique_denom_fraction or count_unique_fraction.
  analyze(
    "ATOXGRLH",
    var_labels = if (show_pct) "Worst toxicity grade, n (%)" else "Worst toxicity grade, n",
    show_labels = "visible",
    indent_mod = 0L,
    afun = a_freq_j,
    extra_args = list(
      .stats = if (show_pct) "count_unique_fraction" else "count_unique",
      denom = "n_rowdf",
      new_levels = list(
        c("Total"),
        list(c("Grade 0", "Grade 1", "Grade 2", "Grade 3", "Grade 4"))
      ),
      new_levels_after = TRUE,
      .indent_mods = 0L,
      restr_columns = c(
        c("GRADE 0", "GRADE 1", "GRADE 2", "GRADE 3", "GRADE 4", "TOTAL")
      )
    )
  )

result <- build_table(lyt, filtered_adlb_tox, alt_counts_df = adsl, round_type = "sas")

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

################################################################################
# Convert to tbl file and output table
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, label_width_ins = 2.4, orientation = "landscape")
