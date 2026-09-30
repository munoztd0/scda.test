###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lsiex02.r
## R version:                 4.5.2
## junco version:             0.1.3
## Short Description:         Listing of Study Treatment Batch Lot Number – [SAD/MAD] [Part 1]
## Disclaimer:                This script is a direct copy of the corresponding Core Standard output identifier. For
##                            SAD/MAD specific changes, refer to tsfvit02b.r, lsidm05.r, and gsfvit02.r for examples of
##                            STUDYPRT filtering, COHORT handling, treatment column structure modifications, pooled
##                            placebo derivations, combined treatment columns, dose-level updates, and other
##                            output-specific structural differences as applicable.
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adex
## Output:                    lsiex02.rtf
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

tblid <- "lsiex02"
fileid <- write_path(opath, tblid)
popfl <- "SAFFL"
trtvar <- "TRT01A"
key_cols <- c("COL0", "COL1")
sort_cols <- c("COL0", "COL1")
disp_cols <- paste0("COL", 0:6)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()


###############################################################################
# Process data
###############################################################################

adex <- haven::read_sas(envsetup::read_path(a_in, "adex.sas7bdat")) |>
  df_na() |>
  filter(!!rlang::sym(popfl) == "Y" & ADOSE > 0) |>
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c("Xanomeline Low Dose", "Xanomeline High Dose", "Placebo")
    )
  )

lsting <- adex |>
  mutate(
    DAEXPDT = as.Date(DAEXPDTC, format = "%Y-%m-%d"),
    DAEXPYR = ifelse(
      stringr::str_length(sub("T.*", "", DAEXPDTC)) >= 4 &
        substr(sub("T.*", "", DAEXPDTC), 1, 4) != "----",
      substr(sub("T.*", "", DAEXPDTC), 1, 4),
      NA
    ),
    DAEXPMO = toupper(month.abb[
      as.numeric(ifelse(
        stringr::str_length(sub("T.*", "", DAEXPDTC)) >= 7 &
          substr(sub("T.*", "", DAEXPDTC), 6, 7) != "--",
        substr(sub("T.*", "", DAEXPDTC), 6, 7),
        NA
      ))
    ]),
    DAEXPDAY = ifelse(
      stringr::str_length(sub("T.*", "", DAEXPDTC)) >= 10 &
        substr(sub("T.*", "", DAEXPDTC), 9, 10) != "--",
      substr(sub("T.*", "", DAEXPDTC), 9, 10),
      NA
    ),
  ) |>
  unite(
    "DAEXPDTL",
    DAEXPDAY,
    DAEXPMO,
    DAEXPYR,
    sep = "",
    na.rm = TRUE,
    remove = FALSE
  ) |>
  mutate(
    COL0 = explicit_na(.data[[trtvar]], ""),
    COL1 = explicit_na(USUBJID, ""),
    # Optional Column: COL2/AVISIT
    COL2 = explicit_na(stringr::str_to_sentence(AVISIT), ""),
    # Optional Column: COL3/ASTDT
    COL3 = ifelse(!is.na(ASTDT), paste0(toupper(format(ASTDT, "%d%b%Y"))), ""),
    # Optional Column: COL4/ASTDTM
    COL4 = ifelse(!is.na(ASTDTM), substr(ASTDTM, 12, 16), ""),
    COL5 = explicit_na(EXLOT, ""),
    # Optional Column: COL6/DAEXPDTC
    COL6 = explicit_na(DAEXPDTL, "")
  ) |>
  arrange(
    COL0,
    COL1,
    !is.na(ASTDT),
    ASTDT,
    !is.na(ASTDTM),
    ASTDTM,
    AVISITN,
    COL2
  )

lsting <- var_relabel(
  lsting,
  COL0 = "Treatment Group",
  COL1 = "Subject ID",
  # Optional Column: COL2/AVISIT
  COL2 = "Visit",
  # Optional Column: COL3/ASTDT
  # Select appropriate column header label
  # COL3 = "Date Dispensed",
  COL3 = "Date Administered",
  # Optional Column: COL4/ASTDTM
  # Select appropriate column header label
  # COL4 = "Time Dispensed",
  COL4 = "Time Administered",
  COL5 = "Batch Lot Number",
  # Optional Column: COL6/DAEXPDTC
  COL6 = "Batch Lot Expiration Date"
)

###############################################################################
# Build listing
###############################################################################

result <- rlistings::as_listing(
  df = lsting,
  key_cols = key_cols,
  sort_cols = sort_cols,
  disp_cols = disp_cols,
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
