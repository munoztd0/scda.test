###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

###############################################################################
## Original Reporting Effort: Standards
## Program Name:              lsids05.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Listing of Overall Study Start and End Dates – [SAD/MAD] [Part 1]
## Disclaimer:                This script is a direct copy of the corresponding Core Standard output identifier. For
##                            SAD/MAD specific changes, refer to tsfvit02b.r, lsidm05.r, and gsfvit02.r for examples of
##                            STUDYPRT filtering, COHORT handling, treatment column structure modifications, pooled
##                            placebo derivations, combined treatment columns, dose-level updates, and other
##                            output-specific structural differences as applicable.
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl
## Output:                    lsids05.rtf
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
library(rtables)
library(rlistings)
library(junco)

###############################################################################
# Define script level parameters
###############################################################################

tblid <- "lsids05"
fileid <- write_path(opath, tblid)
popfl <- "SCRNFL"
key_cols <- c("COL1")
sort_cols <- c("COL1")
disp_cols <- paste0("COL", 1:2)
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

###############################################################################
# Process data
###############################################################################

adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) |>
  df_na() |>
  filter(!!rlang::sym(popfl) == "Y") |>
  summarise(
    rficdt_min = toupper(format(min(RFICDT, na.rm = TRUE), "%d%b%Y")),
    lstsvdt_max = toupper(format(max(LASTCTDT, na.rm = TRUE), "%d%b%Y"))
  )

lsting <- adsl |>
  mutate(
    COL1 = rficdt_min,
    COL2 = lstsvdt_max,
  )

lsting <- var_relabel(
  lsting,
  COL1 = "First Contact in the Study~[super a]",
  COL2 = "Last Contact/Visit in the Study~[super b]"
)

###############################################################################
# Build listing
###############################################################################

result <- rlistings::as_listing(
  df = lsting,
  key_cols = key_cols,
  disp_cols = disp_cols,
  sort_cols = sort_cols,
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
