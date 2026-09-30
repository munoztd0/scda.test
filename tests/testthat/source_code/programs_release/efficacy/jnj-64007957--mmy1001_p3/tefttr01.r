###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tefttr01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefttr01: Time to Response; Full Analysis Set
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADTTEEF.
## Output:                    TEFTTR01.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
##  Rev #:                    1
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
library(haven)

################################################################################
# Define script level parameters:
################################################################################

################################################################################
# - Define output ID and file location
# - Define treatment variable used (default=TRT01P)
# - Define population flag used (default=FASFL)
# - Define control treatment arm
# - Define PARAMCD assignments for time to response and time to best response
################################################################################

tblid <- "TEFTTR01"
fileid <- write_path(opath, tblid)

# Current workaround needed to get correct title:
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
# Warning: Title file should contain exactly one title record per Table ID
# Therefore need to make sure we only have one title record:
tab_titles$title <- tab_titles$title[1]
tab_titles$main_footer <- tab_titles$main_footer[c(1, 2, 3, 4, 5)]

trtvar <- "TRT01A"
trtvarn <- "TRT01AN"
popfl <- "AEFFICFL"
ctrl_grp <- "Phase 1"

time_to_response_paramcd_adeff <- "IRCBRESP"

time_to_response_paramcd_adtte <- "IRCTTRM"
time_to_best_response_paramcd_adtte <- "IRCTTBRM"
time_to_vgpr_or_better_paramcd_adtte <- "IRCTTVRM"
time_to_cr_or_better_paramcd_adtte <- "IRCTTCRM"

################################################################################
# Read in ADSL, ADTTE SAS datasets
################################################################################
adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  select(STUDYID, USUBJID, all_of(trtvar), all_of(popfl))

# Read in ADEFF to get number of subjects who were responders
adeff <- haven::read_sas(read_path(a_in, "adeff.sas7bdat")) |>
  filter(
    PARAMCD == time_to_response_paramcd_adeff &
      (grepl("sCR)", AVALC) |
        grepl("CR)", AVALC) |
        grepl("VGPR)", AVALC) |
        grepl("PR)", AVALC))
  ) |> # Modify as needed
  mutate(RESPONDER = "Y") |>
  select(USUBJID, RESPONDER, AVALC)

# Read in ADTTE
adtte <- haven::read_sas(read_path(a_in, "adtteef.sas7bdat")) |>
  mutate(PARAMCD = as.factor(PARAMCD))

# join data together
ttr <- adtte |>
  right_join(adeff, by = c("USUBJID")) |>
  mutate(
    RESPONSE = RESPONDER == "Y",
    AVAL1 = case_when(
      (PARAMCD == time_to_response_paramcd_adtte &
        EVNTDESC == "first response of PR or better" &
        RESPONDER == "Y") ~
        AVAL
    ),
    AVAL2 = case_when(
      PARAMCD == time_to_best_response_paramcd_adtte &
        CNSR == 0 &
        RESPONDER == "Y" ~
        AVAL
    ),
    AVAL3 = case_when(
      PARAMCD == time_to_vgpr_or_better_paramcd_adtte &
        CNSR == 0 &
        RESPONDER == "Y" ~
        AVAL
    ),
    AVAL4 = case_when(
      PARAMCD == time_to_cr_or_better_paramcd_adtte &
        CNSR == 0 &
        RESPONDER == "Y" ~
        AVAL
    )
  )

# Extra processing, most likely study specific
ttr <- ttr |>
  ## study specific filtering
  filter(STDYPRTN != 3 | (STDYPRTN == 3 & COHTAFL == "Y")) |>
  ###########################################################################################
  ##### study specific section - can be removed if study is just using TRT01A (or any TRT variable that exists in ADSL) as it is in ADSL
  ###########################################################################################
  mutate(
    !!trtvar := case_when(
      STDYPRTN %in% c(1, 2) & RP2DFL == "Y" ~ "Phase 1", # Phase 1 RP2D
      STDYPRTN %in% c(3) & !!rlang::sym(trtvarn) == 301 ~ "Phase 2 Cohort A", # Phase 2 Cohort A
      TRUE ~ NA
    ),
    !!trtvarn := case_when(
      STDYPRTN %in% c(1, 2) & RP2DFL == "Y" ~ 11, # Phase 1 RP2D
      STDYPRTN %in% c(3) & !!rlang::sym(trtvarn) == 301 ~ 12, # Phase 2 Cohort A
      TRUE ~ NA
    )
  ) |>
  mutate(!!trtvar := as.factor(.data[[trtvar]]))

# replace treatment var from ADSL now they have been re-defined and keep all subjects
ttrx <- adtte |>
  filter(STDYPRTN != 3 | (STDYPRTN == 3 & COHTAFL == "Y")) |>
  mutate(
    !!trtvar := case_when(
      STDYPRTN %in% c(1, 2) & RP2DFL == "Y" ~ "Phase 1", # Phase 1 RP2D
      STDYPRTN %in% c(3) & !!rlang::sym(trtvarn) == 301 ~ "Phase 2 Cohort A", # Phase 2 Cohort A
      TRUE ~ NA
    ),
    !!trtvarn := case_when(
      STDYPRTN %in% c(1, 2) & RP2DFL == "Y" ~ 11, # Phase 1 RP2D
      STDYPRTN %in% c(3) & !!rlang::sym(trtvarn) == 301 ~ 12, # Phase 2 Cohort A
      TRUE ~ NA
    )
  ) |>
  mutate(!!trtvar := as.factor(.data[[trtvar]]))

newtrts <- ttrx |>
  group_by(USUBJID) |>
  slice(1) |>
  ungroup() |>
  filter(!is.na(!!rlang::sym(trtvar))) |>
  select(USUBJID, all_of(trtvar), all_of(trtvarn))

adsl <- adsl |>
  select(-all_of(trtvar))

adsl <- left_join(newtrts, adsl, by = c("USUBJID"))
###########################################################################################

# keep only required variables
ttr <- ttr |>
  select(USUBJID, PARAMCD, RESPONSE, AVAL1, AVAL2, AVAL3, AVAL4)

adsl$colspan_trt <- factor(
  ifelse(adsl[[trtvar]] == ctrl_grp, " ", "RP2D"),
  levels = c(" ", "RP2D")
)

adslttr <- ttr |> inner_join(adsl, by = c("USUBJID"))

# colspan_trt_map <- create_colspan_map(adsl,
#                                       non_active_grp = ctrl_grp,
#                                       non_active_grp_span_lbl = " ",
#                                       active_grp_span_lbl = "Active Study Agent",
#                                       colspan_var = "colspan_trt",
#                                       trt_var = trtvar)

# using this instead of standard active vs control spanner
colspan_trt_map <- data.frame(
  colspan_trt = c(" ", "RP2D"),
  TRT01A = c("Phase 1", "Phase 2 Cohort A")
)

################################################################################
# Define layout and build table:
################################################################################

lyt <- rtables::basic_table(
  top_level_section_div = " ",
  show_colcounts = TRUE,
  colcount_format = "N=xx"
) |>
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) |>
  split_cols_by(trtvar) |>
  add_overall_col("Total") |>
  analyze(
    "RESPONSE",
    show_labels = "hidden",
    afun = function(df) {
      in_rows(
        length(unique(df$USUBJID[df$RESPONSE])),
        .formats = "xx.",
        .labels = "Subjects who achieved response"
      )
    }
  ) |>
  analyze(
    "AVAL1",
    nested = FALSE,
    table_names = "AVAL1x",
    var_labels = "Time to First Response (months)~[super a]",
    show_labels = "visible",
    indent_mod = 0L,
    afun = function(x) {
      list(
        "N" = rcell(length(x), format = jjcsformat_xx("xx"))
      )
    }
  ) |>
  analyze(
    "AVAL1",
    nested = TRUE,
    var_labels = "Time to First Response (months)~[super a]",
    show_labels = "hidden",
    indent_mod = 2L,
    afun = function(x) {
      list(
        "Mean (SD)" = rcell(
          c(mean(x), sd(x)),
          format = jjcsformat_xx("xx.xx (xx.xxx)")
        ),
        "Median" = rcell(median(x), format = jjcsformat_xx("xx.xx")),
        "Min, max" = rcell(
          c(min(x), max(x)),
          format = jjcsformat_xx("xx.x, xx.x")
        )
      )
    }
  ) |>
  analyze(
    "AVAL2",
    nested = FALSE,
    table_names = "AVAL2x",
    var_labels = "Time to Best Response (months)~[super a]",
    show_labels = "visible",
    indent_mod = 0L,
    afun = function(x) {
      list(
        "N" = rcell(length(x), format = jjcsformat_xx("xx"))
      )
    }
  ) |>
  analyze(
    "AVAL2",
    nested = TRUE,
    var_labels = "Time to Best Response (months)~[super a]",
    show_labels = "hidden",
    indent_mod = 2L,
    afun = function(x) {
      list(
        "Mean (SD)" = rcell(
          c(mean(x), sd(x)),
          format = jjcsformat_xx("xx.x (xx.xxx)")
        ),
        "Median" = rcell(median(x), format = jjcsformat_xx("xx.xx")),
        "Min, max" = rcell(
          c(min(x), max(x)),
          format = jjcsformat_xx("xx.x, xx.x")
        )
      )
    }
  )
## add optional parameters if required

# Optional parameter 1
lyt <- lyt |>
  analyze(
    "AVAL3",
    nested = FALSE,
    table_names = "AVAL3x",
    var_labels = "Time to VGPR or better (months)~[super a]",
    show_labels = "visible",
    indent_mod = 0L,
    afun = function(x) {
      list(
        "N" = rcell(length(x), format = jjcsformat_xx("xx"))
      )
    }
  ) |>
  analyze(
    "AVAL3",
    nested = TRUE,
    var_labels = "Time to VGPR or better (months)~[super a]",
    show_labels = "hidden",
    indent_mod = 2L,
    afun = function(x) {
      list(
        "Mean (SD)" = rcell(
          c(mean(x), sd(x)),
          format = jjcsformat_xx("xx.xx (xx.xxx)")
        ),
        "Median" = rcell(median(x), format = jjcsformat_xx("xx.xx")),
        "Min, max" = rcell(
          c(min(x), max(x)),
          format = jjcsformat_xx("xx.x, xx.x")
        )
      )
    }
  )

# Optional parameter 2
lyt <- lyt |>
  analyze(
    "AVAL4",
    nested = FALSE,
    table_names = "AVAL4x",
    var_labels = "Time to CR or better (months)~[super a]",
    show_labels = "visible",
    indent_mod = 0L,
    afun = function(x) {
      list(
        "N" = rcell(length(x), format = jjcsformat_xx("xx"))
      )
    }
  ) |>
  analyze(
    "AVAL4",
    nested = TRUE,
    var_labels = "Time to CR or better (months)~[super a]",
    show_labels = "hidden",
    indent_mod = 2L,
    afun = function(x) {
      list(
        "Mean (SD)" = rcell(
          c(mean(x), sd(x)),
          format = jjcsformat_xx("xx.xx (xx.xxx)")
        ),
        "Median" = rcell(median(x), format = jjcsformat_xx("xx.xx")),
        "Min, max" = rcell(
          c(min(x), max(x)),
          format = jjcsformat_xx("xx.x, xx.x")
        )
      )
    }
  )


lyt <- lyt |>
  append_topleft("Parameter")

result <- build_table(lyt, adslttr, alt_counts_df = adsl)

################################################################################
# Add titles and footnotes:
################################################################################

# Add title and main footnotes.
result <- set_titles(result, tab_titles)


################################################################################
# Convert to tbl file and output table
################################################################################
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
