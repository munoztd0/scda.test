###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort:     Standards
## Program Name:              tefttr02.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Program to create tefttr02: Time to Response; Full Analysis Set
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     ADSL, ADTTEEF.
## Output:                    TEFTTR02.rtf
## Remarks:                   Template R script version using rtables framework
##
## Modification History:
## Rev #:                     1
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
## Description:
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

tblid <- "TEFTTR02"
fileid <- write_path(opath, tblid)

tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

trtvar <- "TRT01P"
popfl <- "FASFL"
ctrl_grp <- "CP"

time_to_response_paramcd_adeff <- "BORIRC"
time_to_response_paramcd_adtte <- "TTR"
time_to_best_response_paramcd_adtte <- "TTBR"
time_to_vgpr_or_better_paramcd_adtte <- "TTVGPR"
time_to_cr_or_better_paramcd_adtte <- "TTCR"

# If above parameters are not available in ADTTE efficacy dataset then you may
# need to derive this somehow by using another parameter and then re-deriving based on this
alternative_time_to_response_paramcd_adtte <- "DORIRC"
alternative_time_to_best_response_paramcd_adtte <- "DORINV"

################################################################################
# Read in ADSL, ADTTE SAS datasets
################################################################################
adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(!!rlang::sym(popfl) == "Y") |>
  mutate(!!trtvar := as.factor(.data[[trtvar]])) |>
  select(STUDYID, USUBJID, RANDDT, all_of(trtvar), all_of(popfl))

adsl$colspan_trt <- factor(
  ifelse(adsl[[trtvar]] == ctrl_grp, " ", "Active Study Agent"),
  levels = c("Active Study Agent", " ")
)

# Read in ADEFF to get number of subjects who were responders
adeff <- haven::read_sas(read_path(a_in, "adeff.sas7bdat")) |>
  filter(
    PARAMCD == time_to_response_paramcd_adeff &
      AVALC %in% c("Partial Response (PR)", "Complete Response (CR)")
  ) |>
  mutate(RESPONDER = "Y") |>
  select(USUBJID, RESPONDER)

adtte <- haven::read_sas(read_path(a_in, "adtteef.sas7bdat")) |>
  mutate(PARAMCD = as.factor(PARAMCD)) |>
  select(USUBJID, PARAMCD, CNSR, AVAL, STARTDT)

# join data together
adtte <- adtte |> right_join(adeff, by = c("USUBJID"))
adtte <- adtte |> right_join(adsl, by = c("USUBJID"))

# check to see if parameters exist to produce mandatory sections Time to First Response and Time to Best Response
# Time to First Response
required1 <- adtte |>
  filter(PARAMCD == time_to_response_paramcd_adtte)

# if it does use it, if not derive
if (length(required1$PARAMCD) == 0) {
  required1 <- adtte |>
    filter(PARAMCD == alternative_time_to_response_paramcd_adtte) |>
    mutate(AVAL = as.numeric(((STARTDT - RANDDT) + 1) / 30.4375))
}

# Time to Best Response
required2 <- adtte |>
  filter(PARAMCD == time_to_best_response_paramcd_adtte)

# if it does use it, if not derive
if (length(required2$PARAMCD) == 0) {
  required2 <- adtte |>
    filter(PARAMCD == alternative_time_to_best_response_paramcd_adtte) |>
    mutate(AVAL = as.numeric(((STARTDT - RANDDT) + 1) / 30.4375))
}

adtte_req <- rbind(required1, required2)

# check to see if parameters exist to produce optional sections Time to [VGPR] or better and/or Time to [CR or better]
# Time to [VGPR] or better
optional1 <- adtte |>
  filter(PARAMCD == time_to_vgpr_or_better_paramcd_adtte)

if (length(optional1$PARAMCD) != 0) {
  adtte_req <- rbind(adtte_req, optional1)
}

# Time to [CR or better]
optional2 <- adtte |>
  filter(PARAMCD == time_to_cr_or_better_paramcd_adtte)

if (length(optional2$PARAMCD) != 0) {
  adtte_req <- rbind(adtte_req, optional2)
}

all <- adtte_req |>
  mutate(
    RESPONSE = RESPONDER == "Y",
    AVAL1 = case_when(
      PARAMCD == time_to_response_paramcd_adtte |
        PARAMCD == alternative_time_to_response_paramcd_adtte ~
        AVAL
    ),
    AVAL2 = case_when(
      PARAMCD == time_to_best_response_paramcd_adtte |
        PARAMCD == alternative_time_to_best_response_paramcd_adtte ~
        AVAL
    )
  )

## add optional parameters if required
if (length(optional1$PARAMCD) != 0) {
  all <- all |>
    mutate(
      AVAL3 = case_when(PARAMCD == time_to_vgpr_or_better_paramcd_adtte ~ AVAL)
    )
}
if (length(optional2$PARAMCD) != 0) {
  all <- all |>
    mutate(
      AVAL4 = case_when(PARAMCD == time_to_cr_or_better_paramcd_adtte ~ AVAL)
    )
}

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
if (length(optional1$PARAMCD) != 0) {
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
}
if (length(optional2$PARAMCD) != 0) {
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
}

lyt <- lyt |>
  append_topleft("Parameter")

result <- build_table(lyt, all, alt_counts_df = adsl)

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

################################################################################
# Convert to tbl file and output table
################################################################################
tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "portrait")
