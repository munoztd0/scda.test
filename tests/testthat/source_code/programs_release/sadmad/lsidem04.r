###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lsidem04.r
## R version:                 4.5.2
## junco version:             0.1.6
## Short Description:         Listing of Drug Screening/Pregnancy Test/Serology/
##                            Alcohol Test – [SAD/MAD] [Part 1]
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adlb
## Output:                    lsidem04.rtf
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

###############################################################################
# Prep environment
###############################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)
library(dplyr)
library(tidyr)
library(rtables)
library(rlistings)
library(junco)

###############################################################################
# Define script level parameters
###############################################################################

tblid <- "LSIDEM04"
fileid <- write_path(opath, tblid)
popfl <- "SAFFL"
trtvar <- "TRT01A"
key_cols <- paste0("COL", 0:2)
concat_sep <- " / "
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
studypart_var <- "PARTC" # need to be changed to STUDYPRT variable
studyprt_val <- "PART 1"
cohort_var <- "COHORTC"
actarm_str <- "SAD"

# change as per study requirements
req_paramcd <- c("OPIATE", "HCG", "ETHANOLU")

disp_cols <- paste0("COL", 0:(3 + length(req_paramcd)))


###############################################################################
# Process data
###############################################################################

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) %>%
  df_na() %>%
  dplyr::filter(
    !!rlang::sym(studypart_var) == studyprt_val,
    !!rlang::sym(cohort_var) != cohort_val,
    grepl(actarm_str, ACTARM, ignore.case = TRUE),
    !!rlang::sym(popfl) == "Y"
  ) %>%
  dplyr::mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = unique(sort(.data[[trtvar]]))
    ),
    !!rlang::sym(cohort_var) := factor(
      .data[[cohort_var]],
      levels = unique(.data[[cohort_var]])
    )
  ) %>%
  dplyr::select(
    STUDYID,
    USUBJID,
    TRT01A,
    COHORT
  )

adlb <- haven::read_sas(read_path(a_in, "adlb.sas7bdat")) %>%
  df_na() %>%
  dplyr::select(
    STUDYID,
    USUBJID,
    AVISIT,
    ADTM,
    PARAM,
    PARAMCD,
    PARAMCON,
    PARCAT1,
    AVAL,
    AVALC
  )

new_cols <- setNames(
  lapply(req_paramcd, function(x) {
    rlang::expr(explicit_na(!!rlang::sym(x), ""))
  }),
  paste0("COL", seq_along(req_paramcd) + 3)
)

lsting <- adlb %>%
  dplyr::filter(
    !(is.na(AVAL) & is.na(AVALC)),
    PARAMCD %in% req_paramcd
  ) %>%
  dplyr::mutate(RESULT = dplyr::coalesce(AVALC, as.character(AVAL))) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(
    id_cols = c(USUBJID, AVISIT, ADTM, STUDYID),
    names_from = PARAMCD,
    values_from = RESULT
  ) %>%
  dplyr::inner_join(adsl, by = c("STUDYID", "USUBJID")) %>%
  dplyr::mutate(
    ATDT = as.Date(ADTM, format = "%Y-%m-%d"),
    ADT_TM = substr(as.character(ADTM), 12, 16),
    ADT_Y = dplyr::if_else(
      stringr::str_length(ADTM) >= 4,
      substr(ADTM, 1, 4),
      NA
    ),
    ADT_M = month.abb[
      as.numeric(dplyr::if_else(
        stringr::str_length(sub("T.*", "", ADTM)) >= 7 &
          substr(sub("T.*", "", ADTM), 6, 7) != "--",
        substr(sub("T.*", "", ADTM), 6, 7),
        NA
      ))
    ],
    ADT_DY = dplyr::if_else(
      stringr::str_length(ADTM) >= 10,
      substr(ADTM, 9, 10),
      NA
    )
  ) %>%
  unite(
    "ADTML",
    ADT_DY,
    ADT_M,
    ADT_Y,
    sep = "",
    na.rm = TRUE,
    remove = FALSE
  ) %>%
  dplyr::mutate(
    COL0 = explicit_na(.data[[trtvar]], ""), # or cohort_var
    COL1 = explicit_na(USUBJID, ""),
    COL2 = explicit_na(AVISIT, ""),
    COL3 = case_when(
      (ADTML != "" | !is.na(ADTML)) & ADT_TM != "" ~
        paste0(toupper(ADTML), concat_sep, ADT_TM),
      (ADTML != "" | !is.na(ADTML)) & ADT_TM == "" ~
        paste0(toupper(ADTML), concat_sep, "--:--"),
      TRUE ~ ""
    ),
    !!!new_cols
  ) %>%
  dplyr::arrange(COL0, COL1, ADTM)

param_labels <- sapply(
  req_paramcd,
  function(x) as.character(unique(adlb$PARAM[adlb$PARAMCD == x]))[1]
)

relabel_list <- c(
  list(
    COL0 = "Treatment Group",
    COL1 = "Subject ID",
    COL2 = "Visit",
    COL3 = "Date/Time of Collection"
  ),
  setNames(
    as.list(param_labels),
    paste0("COL", seq_along(req_paramcd) + 3)
  )
)

lsting <- do.call(var_relabel, c(list(lsting), relabel_list))

###############################################################################
# Build listing
###############################################################################

result <- rlistings::as_listing(
  df = lsting,
  key_cols = key_cols,
  disp_cols = disp_cols,
  sort_cols = NULL,
  round_type = "sas"
)

###############################################################################
# Add titles and footnotes
###############################################################################

result <- set_titles(result, tab_titles)

###############################################################################
# Output listing
###############################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
