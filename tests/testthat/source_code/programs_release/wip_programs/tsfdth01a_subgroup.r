###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsfdth01a_subgroup.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tsfdth01: Table of deaths
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL
## Output:                    TSFDTH01.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:                     1
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
## Description:
################################################################################

### example on how a shell can be converted into a subgroup table, when it cannot be done as a single table
### analyze call followed by split_rows_by are causing problems with single table approach
## adding split_rows_by around the table shell does not result in the proper structure

## current approach is a work around, if revising the core shell is too complicated (eg tsfdth01)

## construct wrapper function to construct table for each level of subgroup variable
## this wrapper function is based on the core shell with split_rows_by(subgroup) variable
## loop over the different subgroup levels - mapply
## combine the different subgroup levels tabletrees, using rtables::rbind

# Question : rtables::rbind, do we have something like a page_by?

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
# - Define output ID and file location
# - Define treatment variable used (default=TRT01A)
# - Define population flag used (default=SAFFL)
# - Define the number of days required for death within xx days (ie how DTHTRTFL and DTHAFTFL were defined)
# - Define the cause of death variable to be used. DTHCAUS is defaulted as the primary cause of death.
# - Choose whether or not you want to present a combined active treatment column (default=TRUE)
# - Choose whether or not you want to present the risk difference columns (default=TRUE)
# - Choose which risk difference method you would like (default=Wald)
# - Define what the control treatment group is for your study (e.g Placebo)
# - Define how to create combined treatment columns (if required)
################################################################################

tblid <- "TSFDTH01"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()


trtvar <- "TRT01A"
popfl <- "SAFFL"

days <- 30
dthcausevar <- "DTHCAUS"

combined_colspan_trt <- TRUE
risk_diff <- TRUE
rr_method <- "wald"
ctrl_grp <- "Placebo"

if (combined_colspan_trt == TRUE) {
  # Set up levels and label for the required combined columns
  add_combo <- add_combo_facet(
    "Combined",
    label = "Combined",
    levels = c("Xanomeline High Dose", "Xanomeline Low Dose")
  )

  # choose if any facets need to be removed - e.g remove the combined column for placebo
  rm_combo_from_placebo <- cond_rm_facets(
    facets = "Combined",
    ancestor_pos = NA,
    value = " ",
    split = "colspan_trt"
  )

  mysplit <- make_split_fun(post = list(add_combo, rm_combo_from_placebo))
}

################################################################################
# Process Data:
################################################################################

sas_vs_rds_check("adsl", a_in)
adsl <- readRDS(read_path(a_in, "adsl.rds")) %>%
  filter(!!rlang::sym(popfl) == "Y") %>%
  # remove DTH60TFL variable from below if Deaths within 60 days section is not required
  select(
    STUDYID,
    USUBJID,
    all_of(trtvar),
    all_of(popfl),
    DTHFL,
    DTHTRTFL,
    DTH60TFL,
    DTHAFTFL,
    all_of(dthcausevar),
    AGEGR1
  )

adsl$colspan_trt <- factor(
  ifelse(adsl[[trtvar]] == "Placebo", " ", "Active Study Agent"),
  levels = c("Active Study Agent", " ")
)

if (risk_diff == TRUE) {
  adsl$rrisk_header <- "Risk Difference (%) (95% CI)"
  adsl$rrisk_label <- paste(adsl[[trtvar]], paste("vs", ctrl_grp))
}

adsl <- adsl %>%
  mutate(
    DTHVAR = !!as.name(dthcausevar),
    DTHVAR = as.factor(DTHVAR)
  )

# re-order so OTHER becomes the last category in the table (if there is an OTHER)
is_other <- adsl %>%
  filter(DTHVAR == "OTHER")

if (length(is_other$DTHVAR) != 0) {
  adsl <- adsl %>%
    mutate(DTHVAR = forcats::fct_relevel(DTHVAR, "OTHER", after = Inf))
}

# Convert Reasons to sentence case
adsl$DTHVAR <- stringr::str_to_sentence(adsl$DTHVAR)
adsl$DTHVAR <- as.factor(adsl$DTHVAR)

colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = ctrl_grp,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

################################################################################
# Define layout and build table:
################################################################################
ref_path <- c("colspan_trt", " ", trtvar, ctrl_grp)


### construct wrapper function for a single subgroup level
### this is based on core template + at one place split_rows_by + summarize_row_groups of subgroup
subgrptable <- function(
  df = adsl,
  subgrpvar,
  subgrplvl,
  label_fstr = "Age Group: %s years",
  ref_path,
  colspan_trt_map,
  combined_colspan_trt,
  rr_method,
  mysplit,
  trtvar,
  risk_diff
) {
  adslx <- subset(df, df[[subgrpvar]] == subgrplvl)
  adslx[[subgrpvar]] <- factor(
    as.character(adslx[[subgrpvar]]),
    levels = subgrplvl
  )

  extra_args_1 <- list(
    riskdiff = FALSE,
    .stats = c("n_altdf"),
    denom = "n_altdf",
    extrablankline = TRUE,
    label_fstr = label_fstr
  )

  extra_args_rr <- list(
    method = rr_method,
    #  drop_levels = TRUE,
    .stats = c("count_unique_denom_fraction"),
    ref_path = ref_path,
    denom = "n_altdf"
  )

  # check if we actually have any deaths so this can be used for the layout
  anydth <- adslx %>%
    filter(DTHFL == "Y")

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

  if (risk_diff == TRUE) {
    lyt <- lyt %>%
      split_cols_by("rrisk_header", nested = FALSE) %>%
      split_cols_by(
        trtvar,
        labels_var = "rrisk_label",
        split_fun = remove_split_levels("Placebo")
      )
  }

  if (length(anydth$DTHFL) != 0) {
    lyt <- lyt %>%
      ### here the subgrpvar processing
      split_rows_by(subgrpvar) %>%
      summarize_row_groups(
        subgrpvar,
        cfun = a_freq_j,
        extra_args = extra_args_1
      ) %>%
      split_rows_by(
        "DTHFL",
        split_fun = keep_split_levels("Y"),
        split_label = "Deaths",
        label_pos = "topleft",
        section_div = " "
      ) %>%
      summarize_row_groups(
        "DTHFL",
        cfun = a_freq_j,
        extra_args = append(extra_args_rr, list(label = "Total deaths"))
      ) %>%
      analyze(
        "DTHVAR",
        var_labels = " ",
        a_freq_j,
        extra_args = append(extra_args_rr, list(drop_levels = TRUE)),
        indent_mod = 0,
        show_labels = "hidden"
      ) %>%
      ### note that this split is back to root, from here onwards add indent_mod to split_rows_by
      split_rows_by(
        "DTHTRTFL",
        split_fun = keep_split_levels("Y"),
        # split_fun = drop_split_levels,
        section_div = " ",
        indent_mod = 1L,
        nested = TRUE
      ) %>%
      summarize_row_groups(
        "DTHTRTFL",
        cfun = a_freq_j,
        extra_args = append(
          extra_args_rr,
          list(
            val = "Y",
            label = paste0("Deaths within ", days, " days of last dose")
          )
        )
      ) %>%
      analyze(
        "DTHVAR",
        var_labels = " ",
        afun = a_freq_j,
        extra_args = append(extra_args_rr, list(drop_levels = TRUE)),
        indent_mod = 0,
        show_labels = "hidden"
      ) %>%
      split_rows_by(
        "DTHAFTFL",
        split_fun = keep_split_levels("Y"),
        section_div = " ",
        indent_mod = 1L
      ) %>%
      summarize_row_groups(
        "DTHAFTFL",
        cfun = a_freq_j,
        extra_args = append(
          extra_args_rr,
          list(label = paste0("Deaths beyond ", days, " days of last dose"))
        )
      ) %>%
      analyze(
        "DTHVAR",
        var_labels = " ",
        afun = a_freq_j,
        extra_args = extra_args_rr,
        indent_mod = 0,
        show_labels = "hidden"
      ) %>%
      # Remove below section if the Deaths within 60 days section is not required
      split_rows_by(
        "DTH60TFL",
        split_fun = keep_split_levels("Y"),
        section_div = " ",
        indent_mod = 1L
      ) %>%
      summarize_row_groups(
        "DTH60TFL",
        cfun = a_freq_j,
        extra_args = append(
          extra_args_rr,
          list(label = "Deaths within 60 days of first dose")
        )
      ) %>%
      analyze(
        "DTHVAR",
        var_labels = " ",
        afun = a_freq_j,
        extra_args = extra_args_rr,
        indent_mod = 0,
        show_labels = "hidden"
      )
  } else {
    lyt <- lyt %>%
      analyze(
        "DTHFL",
        a_freq_j,
        show_labels = "hidden",
        extra_args = append(extra_args_rr, list(label = "Total deaths"))
      )
  }

  lyt <- lyt %>%
    append_topleft("  Cause of Death, n (%)")

  ### ensure the original dataframe is used for alt_counts_df, all subtables will have the correct colcounts
  result <- build_table(lyt, adslx, df)

  # If there is no deaths remove top row and display "No data to display" text
  if (length(anydth$DTHFL) == 0) {
    result <- safe_prune_table(
      result,
      prune_func = remove_rows(removerowtext = "Total deaths")
    )
  }

  #########################################################################################
  # Post-Processing step to sort by descending count on chosen active treatment columns.
  # Default is the last treatment (inc. Combined if applicable) under the active treatment
  # spanning header (defaulted to colspan_trt variable). See function documentation for
  # jj_complex_scorefun should your require a different sorting behavior.
  #########################################################################################

  if (length(anydth$DTHFL) != 0) {
    result <- sort_at_path(
      result,
      c(subgrpvar, "*", "DTHFL", "*", "DTHVAR"),
      scorefun = jj_complex_scorefun()
    )
    ### due to layout, from here back to root
    result <- sort_at_path(
      result,
      c("root", "DTHTRTFL", "*", "DTHVAR"),
      scorefun = jj_complex_scorefun()
    )
    result <- sort_at_path(
      result,
      c("root", "DTHAFTFL", "*", "DTHVAR"),
      scorefun = jj_complex_scorefun()
    )
    result <- sort_at_path(
      result,
      c("root", "DTH60TFL", "*", "DTHVAR"),
      scorefun = jj_complex_scorefun()
    )
  }

  return(result)
}


### set up for subgroup
subgrpvar <- "AGEGR1"
subgrplvl <- levels(adsl[[subgrpvar]])
label_fstr <- "Age Group: %s years"

## get all arguments
fn_args <- list(
  df = adsl,
  subgrpvar = subgrpvar,
  label_fstr = label_fstr,
  ref_path = ref_path,
  colspan_trt_map = colspan_trt_map,
  combined_colspan_trt = combined_colspan_trt,
  rr_method = rr_method,
  mysplit = mysplit,
  trtvar = trtvar,
  risk_diff = risk_diff
)

## construct each subgroup table
all_results <- mapply(
  subgrptable,
  subgrplvl,
  MoreArgs = fn_args,
  SIMPLIFY = FALSE
)

# if problem with mapply
# all_results <- list()
# for (i in 1:length(subgrplvl)){
#   all_results[[i]] <- do.call(subgrptable, args = append(fn_args, list(subgrplvl = subgrplvl[[i]])))
# }

## combine results
result <- all_results[[1]]

for (i in 2:length(all_results)) {
  result <- rbind(result, all_results[[i]])
}


## Remove the N=xx column headers for the risk difference columns
result <- remove_col_count(result)

head(result, 30)


################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

################################################################################
# Convert to tbl file and output table
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")


wald_quick_check <- function(
  n1_yes,
  n1_tot,
  n2_yes,
  n2_tot,
  conf_level = 0.95
) {
  stat_propdiff_ci(
    x = list(n1_yes),
    y = list(n2_yes),
    N_x = n1_tot,
    N_y = n2_tot,
    conf_level = conf_level
  )
}

wald_quick_check(n1_yes = 0, n1_tot = 132, n2_yes = 2, n2_tot = 117)
wald_quick_check(n1_yes = 4, n1_tot = 132, n2_yes = 6, n2_tot = 117)
