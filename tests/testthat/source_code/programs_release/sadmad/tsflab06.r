###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsflab06.r
## R version:                 4.5.2
## junco version:             0.1.3
## Short Description:         Shift in [Laboratory Category] Laboratory Values From Baseline to [Time Point] –
##                            [SAD/MAD][Part 1]
## Disclaimer:                This script is a direct copy of the corresponding Core Standard output identifier. For
##                            SAD/MAD specific changes, refer to tsfvit02b.r, lsidm05.r, and gsfvit02.r for examples of
##                            STUDYPRT filtering, COHORT handling, treatment column structure modifications, pooled
##                            placebo derivations, combined treatment columns, dose-level updates, and other
##                            output-specific structural differences as applicable.
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adlb or adlbc or adlc
## Output:                    tsflab06chm.rtf
##                            tsflab06hem.rtf
## Remarks:
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

tblid <- "tsflab06"
fileid <- write_path(opath, tblid)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

popfl <- "SAFFL"

# Actual treatment variable (default=TRT01A).
trtvar <- "TRT01A"

ctrl_grp <- "Placebo"

ad_domain <- "ADLB"

## select a Single time point only — set to NULL to include all visits
selvisit <- NULL #c("Cycle 02")

# PARCAT1 categories to produce: CHEMISTRY -> CHM, HEMATOLOGY -> HM
parcat1_categories <- list(
  chm = "CHEMISTRY",
  hem = "HEMATOLOGY"
)

tblid_chm <- paste0(tblid, "chm")
tblid_hem <- paste0(tblid, "hem")

################################################################################
# Process Data:
################################################################################

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  df_na() |>
  filter(.data[[popfl]] == "Y") |>
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c("Xanomeline Low Dose", "Xanomeline High Dose", "Placebo")
    )
  ) |>
  select(STUDYID, USUBJID, all_of(c(popfl, trtvar)))

adlb_complete <- haven::read_sas(read_path(a_in, paste0(tolower(ad_domain), ".sas7bdat"))) |>
  df_na()

flagvars <- c("ONTRTFL", "TRTEMFL", "LVOTFL")
adlb00 <- adlb_complete |>
  # Filter to CHEMISTRY and HEMATOLOGY only
  filter(toupper(PARCAT1) %in% toupper(unlist(parcat1_categories))) |>
  filter(ANL02FL == "Y", if (!is.null(selvisit)) AVISIT %in% selvisit else TRUE) |>
  filter(!is.na(ANRIND)) |>
  filter(!is.na(BNRIND)) |>
  mutate(AVISIT = stringr::str_to_sentence(AVISIT)) |>
  mutate(
    AVISIT = factor(
      .data[['AVISIT']],
      levels = unique(.data[['AVISIT']])[order(unique(.data[['AVISITN']]))]
    )
  ) |>
  select(
    USUBJID,
    AVISITN,
    AVISIT,
    PARAMCD,
    PARAM,
    PARAMN,
    PARCAT1,
    PARCAT3,
    PARCAT3N,
    ANRIND,
    BNRIND,
    ANRLO,
    ANRHI,
    starts_with("ANL"),
    all_of(flagvars),
    LBSEQ,
    AVAL,
    AVALC
  )

# Sort by PARCAT1, PARCAT3N, then PARAM alphabetically
adlb00 <- adlb00 |>
  arrange(PARCAT1, PARCAT3N, PARCAT3, PARAM)

# Set PARAM factor levels in PARCAT1 -> PARCAT3N -> PARAM order
params_ordered <- unique(as.character(adlb00$PARAM))

adlb00 <- adlb00 |>
  mutate(PARAM = factor(as.character(PARAM), levels = params_ordered))

adlb00 <- var_relabel_list(adlb00, var_labels(adlb_complete, fill = T))

filtered_adlb <- adlb00 |>
  filter(ANL02FL == "Y", if (!is.null(selvisit)) AVISIT %in% selvisit else TRUE) |>
  filter(!is.na(ANRIND)) |>
  filter(!is.na(BNRIND)) |>
  inner_join(adsl)

## trick for alt_counts_df to work with col splitting
# add BNRIND to adsl, all assign to extra level N (column will be used for N counts)
adslx <- adsl |>
  mutate(BNRIND = "N") |>
  mutate(
    BNRIND = factor(
      BNRIND,
      levels = c("N", "LOW", "NORMAL", "HIGH"),
      labels = c("N", "Low", "Normal", "High")
    )
  )

### make factor var and add extra level to BNRIND to be used as N column
filtered_adlb <- filtered_adlb |>
  mutate(
    BNRIND = factor(
      as.character(BNRIND),
      levels = c("N", "LOW", "NORMAL", "HIGH"),
      labels = c("N", "Low", "Normal", "High")
    )
  ) |>
  mutate(
    ANRIND = factor(
      as.character(ANRIND),
      levels = c("LOW", "NORMAL", "HIGH"),
      labels = c("Low", "Normal", "High")
    )
  )

## add variable for column split header
filtered_adlb$BNRIND_header <- "Baseline"
adslx$BNRIND_header <- "Baseline"

filtered_adlb$BNRIND_header2 <- " " ## first column N should not appear under Baseline column span
adslx$BNRIND_header2 <- " " ## first column N should not appear under Baseline column span

################################################################################
# define a  mapping for low/normal/high, tests that do not have both directions are to be identified
################################################################################
low_high_map <- unique(
  filtered_adlb |>
    mutate(
      xANRLO = !is.na(ANRLO),
      xANRHI = !is.na(ANRHI)
    ) |>
    select(PARCAT1, PARCAT3N, PARCAT3, PARAMCD, PARAM, xANRLO, xANRHI) |>
    group_by(PARAMCD) |>
    mutate(xANRLO = any(xANRLO), xANRHI = any(xANRHI)) |>
    ungroup()
)

# check if there are tests that only have 1 direction
lh_1 <- nrow(
  low_high_map |>
    filter(!(xANRLO & xANRHI))
) >
  0

if (lh_1) {
  low_high <- low_high_map |>
    mutate(ANRIND = "LOW") |>
    mutate(ANRIND = factor(ANRIND, levels = c("LOW", "NORMAL", "HIGH"))) |>
    mutate(PARAMCD = droplevels(PARAMCD)) |>
    tidyr::expand(., PARAMCD, ANRIND)

  low_high_map <- low_high_map |>
    full_join(low_high) |>
    mutate(
      todel = case_when(
        ANRIND == "LOW" & !xANRLO ~ TRUE,
        ANRIND == "HIGH" & !xANRHI ~ TRUE,
        TRUE ~ FALSE
      )
    ) |>
    filter(!todel) |>
    select(-c(todel, xANRLO, xANRHI))
}

### here no such tests, so there is no need to work with a mapping table for now

################################################################################
# Define layout and build table:
################################################################################

################################################################################
# Core function to produce shell for specific parcat1 selection
################################################################################

build_result_parcat1 <- function(
  df = filtered_adlb,
  PARCAT1sel = NULL,
  .adsl = adslx,
  map = low_high_map,
  tblid,
  .trtvar = trtvar
) {
  lyt_filter <- function(PARCAT1sel = NULL, map) {
    if (!is.null(PARCAT1sel)) {
      map <- map |>
        filter(toupper(PARCAT1) %in% toupper(PARCAT1sel))

      df_filtered <- df |>
        filter(toupper(PARCAT1) %in% toupper(PARCAT1sel))
    } else {
      df_filtered <- df
    }

    lyt <- basic_table(show_colcounts = FALSE) |>
      ## to ensure N column is not under the Baseline column span header
      split_cols_by("BNRIND_header2") |>
      split_cols_by("BNRIND", split_fun = keep_split_levels("N")) |>
      split_cols_by("BNRIND_header", nested = FALSE) |>
      split_cols_by(
        "BNRIND",
        split_fun = make_split_fun(
          pre = list(rm_levels(excl = "N")),
          post = list(add_overall_facet("TOTAL", "Total"))
        )
      ) |>
      #### replace split_rows and summarize by single analyze call
      ### a_freq_j only works due to
      ### special arguments can do the trick : denomf = adslx & .stats = count_unique
      ### we want counts of treatment group coming from adsl, not from input dataset, therefor, countsource = altdf
      analyze(
        vars = .trtvar,
        afun = a_freq_j,
        extra_args = list(
          restr_columns = "N",
          .stats = "count_unique",
          countsource = "altdf",
          extrablankline = TRUE
        ),
        indent_mod = -1L
      ) |>
      ## main part of table, restart row-split so nested = FALSE
      split_rows_by(
        "PARAM",
        nested = FALSE,
        label_pos = "topleft",
        child_labels = "visible",
        split_label = "Laboratory Test",
        split_fun = drop_split_levels,
        section_div = " "
      ) |>
      split_rows_by(
        "AVISIT",
        label_pos = "topleft",
        split_label = "Study Visit",
        section_div = " "
      ) |>
      split_rows_by(
        .trtvar,
        label_pos = "hidden",
        split_label = "Treatment Group",
        section_div = " "
      ) |>
      ### a_freq_j
      ### the special statistic "n_rowdf" option does the trick here of getting the proper value for the N column
      summarize_row_groups(
        .trtvar,
        cfun = a_freq_j,
        extra_args = list(
          .stats = "n_rowdf",
          restr_columns = c("N")
        )
      ) |>
      ## add extra level TOTAL using new_levels, rather than earlier technique
      ## advantage for denominator derivation -- n_rowdf can be used, if we'd like to present fraction as well
      ## switch .stats to count_unique_denom_fraction or count_unique_fraction
      analyze(
        "ANRIND",
        afun = a_freq_j,
        extra_args = list(
          .stats = "count_unique",
          denom = "n_rowdf",
          new_levels = list(c("Total"), list(c("Low", "Normal", "High"))),
          new_levels_after = TRUE,
          .indent_mods = 1L,
          restr_columns = c(
            c("LOW", "NORMAL", "HIGH", "TOTAL")
          )
        )
      )
  }

  lyt <- lyt_filter(PARCAT1sel = PARCAT1sel, map = map)

  if (!is.null(PARCAT1sel)) {
    df <- df |>
      filter(toupper(PARCAT1) %in% toupper(PARCAT1sel))
  }

  if (nrow(df) > 0) {
    result <- build_table(lyt, df, alt_counts_df = .adsl, round_type = "sas")
  } else {
    result <- NULL
    message(paste0(
      "Parcat1 [",
      PARCAT1sel,
      "] is not present on input dataset"
    ))
    return(result)
  }

  ################################################################################
  # Set title
  ################################################################################

  return(result)
}

################################################################################
# Apply core function: one RTF per PARCAT1 category
#  - Chemistry (CHM)
#  - Hematology (HM)
# Data sorted by PARCAT1 -> PARCAT3N -> PARAM (alphabetical)
################################################################################
### note: the same core tblid (tsflab06) is used; abbreviation (CHM/HM) is appended via get_tblid
result_chm <- build_result_parcat1(
  PARCAT1sel = parcat1_categories[["chm"]],
  tblid = tblid
)

result_hem <- build_result_parcat1(
  PARCAT1sel = parcat1_categories[["hem"]],
  tblid = tblid
)

################################################################################
# Add titles and footnotes:
################################################################################

result_chm <- set_titles(result_chm, tab_titles)
result_hem <- set_titles(result_hem, tab_titles)

# If we have table id specific tile and footer uncomment this
# result_chm <- set_titles(result_chm, get_titles_internal(tblid_chm))
# result_hem <- set_titles(result_hem, get_titles_internal(tblid_hem))

################################################################################
# Convert to tbl file and output table
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result_chm, file = write_path(opath, tblid_chm), orientation = "landscape")
tt_to_tlgrtf(string_map = string_map, tt = result_hem, file = write_path(opath, tblid_hem), orientation = "landscape")
