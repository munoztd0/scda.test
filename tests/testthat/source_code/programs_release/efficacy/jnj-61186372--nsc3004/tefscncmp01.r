###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefscncmp01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefscncmp01:  Compliance of Scans;
##                            Full Analysis Set (Study mmy, bc, lc)
## Author:                    Technology Solutions
## Date:                      2026-09-302025
## Input:                     ADSL, ADPRO, SV
## Output:                    TEFSCNCMP01.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:                     1
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
## Description:
################################################################################

# Environment ----

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)


library(dplyr)
library(rtables)
library(stringr)
library(junco)
library(haven)

# Parameters ----

# Development version or production version to be used?
# Determines path for use of SV data set below.
dev_version <- TRUE

# Define output ID and file location.
tblid <- "TEFSCNCMP01"
fileid <- write_path(opath, tblid)

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
tab_titles$title <- tab_titles$title[1]

# Define treatment variable used (default=TRT01P).
trtvar <- "TRT01P"

# Define control group label used in the treatment variable.
ctrlab <- "CP"

# Define population flags used.
popfl <- "FASFL"

# Define the PARAMCD value where we want to check compliance.
comppar <- "PGI0101"

# Additional domain flags.
domfl <- quote(ABLFL == "Y" | APOBLFL == "Y")

# Formatting options.
formats <- list(
  count_percent = jjcsformat_count_fraction
)

# Data ----

## ADSL ----

adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y")) |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl))

adsl$colspan_trt <- factor(
  ifelse(adsl[[trtvar]] == ctrlab, " ", "Active Study Agent"),
  levels = c("Active Study Agent", " ")
)

colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = ctrlab,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

## ADPRO ----

adpro <- haven::read_sas(read_path(a_in, "adpro.sas7bdat")) |>
  filter(if_all(all_of(popfl), ~ .x == "Y"), !!domfl) |>
  filter(PARAMCD == comppar) |>
  mutate(AVISIT = ifelse(ABLFL == "Y", "Baseline", AVISIT)) |>
  select(
    STUDYID,
    USUBJID,
    AVISIT,
    AVISITN,
    all_of(popfl),
    PARAMCD,
    ABLFL,
    APOBLFL
  )

## SV ----

# Note: This code might have to be adapted significantly for the study.

sv_path <- if (dev_version) {
  file.path(
    "/adr/PREPROD/pharma/globalcode_development/jnj-61186372--nsc3004/dbr_csr_final/_source",
    "sv.sas7bdat"
  )
} else {
  read_path(d_in, "sv.sas7bdat")
}
sv <- haven::read_sas(sv_path)

# Step 1: Identify Expected Visits.
exp <- adsl |>
  select(USUBJID, all_of(popfl)) |>
  inner_join(sv, by = "USUBJID") |>
  mutate(
    PARAM = "PRO",
    PARAMCD = "PRO",
    AVISITN = VISITNUM,
    AVISIT = stringi::stri_trans_totitle(VISIT),
    AVAL = 1,
    AVALC = "1",
    SVTC = stringr::str_sub(SVSTDTC, 1, 10),
    ADTC = case_when(
      nchar(SVTC) == 10 & !str_starts(SVTC, "SITE|HOME") ~ SVTC,
      nchar(SVTC) == 7 & !str_starts(SVTC, "UN") ~
        paste0(
          stringr::str_remove_all(SVTC, "-"),
          "-01"
        ),
      nchar(SVTC) == 4 & !str_starts(SVTC, "III") ~
        paste0(
          stringr::str_remove_all(SVTC, "-"),
          "-01-01"
        ),
      TRUE ~ NA_character_
    ),
    ADT = lubridate::ymd(ADTC)
  ) |>
  select(USUBJID, AVISITN, AVISIT, PARAM, PARAMCD, AVAL, AVALC, ADT)

# Step 2: Filter Expected Visits.
exp1 <- exp |>
  filter(AVISITN %in% c(2001, 4001, 30000)) |>
  mutate(
    AVISIT = if_else(AVISITN == 2001, "Baseline", AVISIT),
    ABLFL = if_else(AVISITN == 2001, "Y", "N"),
    APOBLFL = if_else(AVISITN > 2001, "Y", "N")
  )

# Step 3: Count for Cycle x Visit.
chk <- exp1 |>
  filter(
    stringr::str_detect(AVISIT, "Cycle"),
    !str_detect(AVISIT, "UNSCHED")
  ) |>
  select(AVISITN, AVISIT)

nvis <- max(chk$AVISITN)

# Step 4: Generate Visit data.
exp2 <- exp1 |>
  mutate(
    AVISIT = if_else(ABLFL == "Y" & APOBLFL != "Y", "Baseline", AVISIT),
    AVISITN = if_else(ABLFL == "Y" & APOBLFL != "Y", 2001, AVISITN)
  ) |>
  filter(
    AVISIT == "Baseline" |
      (AVISITN >= 4001 & AVISITN <= nvis & (AVISITN - 4001) %% 1000 == 0)
  )

# Sort visit data.
evist <- exp2 |>
  arrange(USUBJID, AVISITN, AVISIT)

# Get End of Treatment Visits.
eot <- sv |>
  filter(tolower(VISIT) == "end of treatment") |>
  mutate(
    AVISITN = 30000,
    AVISIT = "End Of Treatment"
  ) |>
  distinct(USUBJID, AVISITN, AVISIT)

# Filter down to population subjects.
eot <- adsl |>
  select(USUBJID, all_of(popfl)) |>
  inner_join(eot, by = "USUBJID")

# Join with expected visits.
evist <- evist |>
  full_join(eot, by = c("USUBJID", "AVISITN", "AVISIT"))

# Step 5: Identify Received Visits.
pro <- adpro |>
  mutate(
    AVISIT = if_else(ABLFL == "Y", "Baseline", AVISIT),
    AVISITN = if_else(ABLFL == "Y", 2001, AVISITN)
  ) |>
  filter(
    AVISIT == "Baseline" |
      (AVISITN >= 4001 & AVISITN <= nvis & (AVISITN - 4001) %% 1000 == 0) |
      AVISITN == 30000
  ) |>
  select(USUBJID, AVISIT, AVISITN)

# Step 6: Flag Expected and Received at Each Visit.
proexp <- evist |>
  left_join(
    pro |> mutate(from_pro = TRUE),
    by = c("USUBJID", "AVISITN", "AVISIT")
  ) |>
  mutate(
    f1fl = "Y", # Expected visit.
    f2fl = if_else(is.na(from_pro), "N", if_else(from_pro, "Y", "N")) # Received visit.
  ) |>
  left_join(adsl |> select(USUBJID, all_of(trtvar)), by = "USUBJID") |>
  mutate(
    ord = case_when(
      AVISITN == 2001 ~ 1,
      AVISITN == 4001 ~ 2,
      AVISITN == 30000 ~ 3
    )
  )

# Step 7: Flag Subjects Expected and Received in Population
apro <- pro |>
  select(USUBJID) |>
  distinct()

aproexp <- adsl |>
  select(USUBJID, all_of(trtvar)) |>
  left_join(
    apro |> mutate(from_apro = TRUE),
    by = "USUBJID"
  ) |>
  mutate(
    a1fl = "Y", # Expected subject
    a2fl = if_else(is.na(from_apro), "N", if_else(from_apro, "Y", "N")) # Received subject
  )

## Analysis ----

# Double check numbers:
table(aproexp$TRT01P)
table(subset(proexp, AVISIT == "Baseline")$TRT01P)
# These are the same, so we can just proceed with proexp below.

ana <- proexp |>
  mutate(
    expected = (f1fl == "Y"),
    received = (f2fl == "Y"),
    missing = (f1fl == "Y" & f2fl == "N"),
    AVISIT = factor(AVISIT),
    CONST = "CONST"
  ) |>
  select(USUBJID, AVISIT, CONST, expected, received, missing, all_of(trtvar))

# Layout ----

lyt <- rtables::basic_table() |>
  split_cols_by(
    trtvar,
    show_colcounts = TRUE,
    colcount_format = "N=xx"
  ) |>
  split_cols_by("CONST", split_fun = cmp_split_fun) |>
  split_rows_by("AVISIT") |>
  summarize_row_groups(
    "USUBJID",
    cfun = cmp_cfun,
    extra_args = list(
      formats = formats,
      variables = list(
        expected = "expected",
        received = "received",
        missing = "missing"
      )
    )
  ) |>
  append_topleft("Visit")

# Output ----

result <- build_table(lyt, ana, alt_counts_df = adsl)

# Add titles and footnotes.

result <- set_titles(result, tab_titles)

# Convert to tbl file and output table.
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
