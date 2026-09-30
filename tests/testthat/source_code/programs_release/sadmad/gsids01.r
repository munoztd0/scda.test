###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsids01.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Disposition of Subjects – [SAD/MAD] [Part 1]
## Disclaimer:                This script is a direct copy of the corresponding Core Standard output identifier. For
##                            SAD/MAD specific changes, refer to tsfvit02b.r, lsidm05.r, and gsfvit02.r for examples of
##                            STUDYPRT filtering, COHORT handling, treatment column structure modifications, pooled
##                            placebo derivations, combined treatment columns, dose-level updates, and other
##                            output-specific structural differences as applicable.
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, addisp
## Output:                    gsids01.rtf
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
# Prep environment -------------------------------------------------------------
################################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(dplyr)
library(tidytlg)
library(junco)
library(ggplot2)
library(stringr)
library(tibble)
library(magick)

################################################################################
# Define output ID  ------------------------------------------------------------
################################################################################

tblid <- "gsids01"

################################################################################
# Get titles and footnotes  ----------------------------------------------------
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Define script level parameters  ----------------------------------------------
################################################################################

trtvar <- "TRT01P"
popfl <- "SCRNFL"
combination_trt <- FALSE


if (combination_trt) {
  stdy_agts <- haven::read_sas(envsetup::read_path(a_in, "addisp.sas7bdat")) |>
    filter(!!rlang::sym(popfl) == "Y" & !is.na(DSSCAT) & DSSCAT != "") |>
    select(DSSCAT) |>
    unique() |>
    unlist() |>
    as.character() |>
    str_to_sentence()
  num_words <- c(
    "zero",
    "one",
    "both",
    "three",
    "four",
    "five",
    "six",
    "seven",
    "eight",
    "nine",
    "ten"
  )

  count_label <- num_words[length(stdy_agts) + 1]
} else {
  stdy_agts <- NULL
}


################################################################################
# Process data -----------------------------------------------------------------
################################################################################

adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) |>
  df_na() |>
  filter(.data[[popfl]] == "Y") |>
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c(
        "Xanomeline Low Dose",
        "Xanomeline High Dose",
        "Placebo"
      )
    )
  ) |>
  select(
    USUBJID,
    all_of(trtvar),
    all_of(popfl),
    EOTSTT,
    DCTREAS,
    EOSSTT,
    DCSREAS,
    RANDFL,
    SCRNFL,
    SCRFFL,
    ENRLFL,
    DCSCREEN,
    RESCRNFL
  )


if (combination_trt) {
  addisp_ <- haven::read_sas(envsetup::read_path(a_in, "addisp.sas7bdat")) |>
    df_na() |>
    filter(!!rlang::sym(popfl) == "Y") |>
    filter(grepl("^(EOTS|DCTS)", PARAMCD)) |>
    select(USUBJID, PARAMCD, AVALC)

  addispa <- addisp_ |>
    tidyr::pivot_wider(
      names_from = c(PARAMCD),
      values_from = AVALC
    ) |>
    df_na() |>
    mutate(
      ONGOING_SUB = case_when(
        rowSums(dplyr::across(dplyr::starts_with("EOTS"), ~ toupper(.x) == "ONGOING"), na.rm = TRUE) > 0 ~ "Y",
        TRUE ~ NA_character_
      ),
      COMPL_SUB = dplyr::case_when(
        rowSums(dplyr::across(dplyr::starts_with("EOTS"), ~ toupper(.x) == "COMPLETED"), na.rm = TRUE) ==
          rowSums(!is.na(dplyr::across(dplyr::starts_with("EOTS")))) &
          rowSums(!is.na(dplyr::across(dplyr::starts_with("EOTS")))) > 0 ~
          "Y",
        TRUE ~ NA_character_
      ),
      DISC_BOTH_TRT = dplyr::case_when(
        rowSums(dplyr::across(dplyr::starts_with("EOTS"), ~ toupper(.x) == "DISCONTINUED"), na.rm = TRUE) ==
          rowSums(!is.na(dplyr::across(dplyr::starts_with("EOTS")))) &
          rowSums(!is.na(dplyr::across(dplyr::starts_with("EOTS")))) > 0 ~
          "Y",
        TRUE ~ NA_character_
      ),
      DISC_ONE_TRT = dplyr::case_when(
        rowSums(dplyr::across(dplyr::starts_with("EOTS"), ~ toupper(.x) == "DISCONTINUED"), na.rm = TRUE) > 0 ~ "Y",
        TRUE ~ NA_character_
      )
    ) |>
    mutate(
      dplyr::across(
        dplyr::starts_with("EOTS"),
        ~ dplyr::case_when(toupper(.x) == "COMPLETED" ~ "Y", TRUE ~ NA_character_),
        .names = "COMPLETED_{gsub('EOTS|STT', '', .col)}"
      ),
      dplyr::across(
        dplyr::starts_with("EOTS"),
        ~ dplyr::case_when(toupper(.x) == "DISCONTINUED" ~ "Y", TRUE ~ NA_character_),
        .names = "DISCONTINUED_{gsub('EOTS|STT', '', .col)}"
      )
    )

  addisp <- inner_join(addispa, adsl, by = c("USUBJID"))
} else {
  addisp <- adsl
}

################################################################################
# CONSORT CONFIGURATION---------------------------------------------------------
################################################################################

CONSORT_CONFIG <- list(
  # data----
  data = addisp,

  # screened----
  screened_flag = "SCRNFL",

  # screen failure----
  show_screen_failure = TRUE,
  screen_failure_flag = "SCRFFL",
  screen_failure_reason = "DCSCREEN",
  show_screen_failure_reasons = TRUE,

  # re-screen----
  show_rescreen = TRUE,
  rescreen_flag = "RESCRNFL",
  show_rescreen_screen_failure = TRUE,
  show_rescreen_rand_enr = TRUE,

  # randomized / enrolled----
  #  rand_enr_var = "RANDFL"/ "ENRLFL"
  #  rand_enr_lbl = "Randomized"/ "Enrolled"
  # -----------------------------------------
  rand_enr_var = "RANDFL",
  rand_enr_lbl = "Randomized",

  # treatment----
  treatment_variable = trtvar,
  combination_trt = combination_trt,

  # non-combination treatment discontinuation reason----

  show_treatment_discontinued_reason = TRUE,
  treatment_discontinued_reason = "DCTREAS",

  # study agent display names----
  study_agent_names = stdy_agts,

  # study agent components----
  show_agent_ongoing = TRUE,
  show_agent_completed = TRUE,
  show_agent_discontinued = TRUE,
  show_agent_discontinued_reason = TRUE,
  show_agent_reasons = TRUE,

  # agent treatment components----
  show_ongoing_sub = TRUE,
  ongoing_sub_variable = "ONGOING_SUB",

  show_compl_sub = TRUE,
  compl_sub_variable = "COMPL_SUB",

  show_disc_both_trt = TRUE,
  disc_both_trt_variable = "DISC_BOTH_TRT",

  show_disc_one_trt = TRUE,
  disc_one_trt_variable = "DISC_ONE_TRT",

  # study column----
  show_study_ongoing = TRUE,
  show_study_completed = TRUE,
  show_study_discontinued = TRUE,
  show_study_reasons = TRUE,

  # percentages----
  show_percentages = TRUE,
  percentage_digits = 1,

  # upper N position----
  # "after": Randomized N =100
  # "before": N = 100 Randomized
  # ----------------------------
  upper_n_position = "after",

  # text wrapping----
  wrap_width = 30
)

################################################################################
# LAYOUT CONFIGURATION----------------------------------------------------------
################################################################################

LAYOUT_CONFIG <- list(
  # screened----
  screened = list(
    x = 0,
    y = 25,
    width = 3.0,
    height = 1.0,
    fill = "white",
    colour = "black",
    linewidth = 0.5,
    text_size = 7,
    lineheight = 0.9
  ),

  # screening split----
  screening_split = list(
    y = 24,
    x_spacing = 8,
    linewidth = 0.5
  ),

  # screening boxes----
  screening = list(
    y = 23,
    width = 3.0,
    height = 1.0,
    fill = "white",
    colour = "black",
    linewidth = 0.5,
    text_size = 7,
    lineheight = 0.9
  ),

  # re-screen children----
  rescreen_children = list(
    y = 21,
    x_spacing = 4,
    width = 2.8,
    height = 0.9,
    fill = "white",
    colour = "black",
    linewidth = 0.5,
    text_size = 7,
    lineheight = 0.9,
    connector_y = 22
  ),

  # screen failure reason box----
  screen_failure_reason_box = list(
    y = 21,
    width = 5,
    min_height = 0,
    line_spacing = 0.3,
    padding = 0.75,
    fill = "white",
    colour = "black",
    linewidth = 0.5,
    text_size = 7,
    lineheight = 0.75
  ),

  # treatment split----
  treatment_split = list(
    connector_y = 19.5,
    x_spacing = 9,
    linewidth = 0.5
  ),

  # treatment box----
  treatment = list(
    y = 18.5,
    width = 5.0,
    height = 1.0,
    fill = "white",
    colour = "black",
    linewidth = 0.5,
    text_size = 7,
    lineheight = 0.9
  ),

  # lower treatment blocks----
  lower_blocks = list(
    x_spacing = 10,
    connector_y = 17.5,
    vertical_distance = 1.5,
    spine_gap = 0.25
  ),

  # lower columns----
  columns = list(
    x_spacing = 4.5,
    vertical_distance = 1.2
  ),

  # agent treatment box----
  subject = list(
    width = 2.0,
    height = 0.85,
    y_spacing = 1.25,
    fill = "white",
    colour = "black",
    linewidth = 0.5,
    text_size = 4,
    lineheight = 0.9
  ),

  # study agent box----
  agent = list(
    width = 2.2,
    height = 0.85,
    y_spacing = 1.25,
    fill = "white",
    colour = "black",
    linewidth = 0.5,
    text_size = 4,
    lineheight = 0.9
  ),

  # treatment status box----
  treatment_status = list(
    width = 3.5,
    height = 0.85,
    y_spacing = 1.25,
    fill = "white",
    colour = "black",
    linewidth = 0.5,
    text_size = 6,
    lineheight = 0.9
  ),

  # study box----
  study = list(
    width = 3.5,
    height = 0.85,
    y_spacing = 1.25,
    fill = "white",
    colour = "black",
    linewidth = 0.5,
    text_size = 6,
    lineheight = 0.9
  ),

  # reason box----
  reason = list(
    width = 4.2,
    min_height = 0,
    line_spacing = 0.3,
    padding = 1,
    gap = 0.5,
    fill = "white",
    colour = "black",
    linewidth = 0.5,
    text_size = 6,
    lineheight = 0.75
  ),

  # lines----
  lines = list(
    linewidth = 0.5,
    lineend = "square"
  ),

  # plot----
  plot = list(
    expand = 2,
    margin = c(20, 20, 20, 20)
  )
)

################################################################################
# VALIDATE OPTIONAL VARIABLE----------------------------------------------------
################################################################################

validate_optional_variable <- function(
  data,
  variable,
  show,
  variable_label
) {
  if (!isTRUE(show)) {
    return(invisible(NULL))
  }

  if (
    is.null(variable) ||
      length(variable) != 1 ||
      is.na(variable) ||
      !variable %in% names(data)
  ) {
    stop(
      paste0(
        "The variable '",
        paste(variable, collapse = ", "),
        "' specified for ",
        variable_label,
        " was not found in the dataset."
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

################################################################################
# VALIDATE REQUIRED VARIABLE----------------------------------------------------
################################################################################

validate_required_variable <- function(
  data,
  variable,
  variable_label
) {
  if (
    is.null(variable) ||
      length(variable) != 1 ||
      is.na(variable) ||
      !variable %in% names(data)
  ) {
    stop(
      paste0(
        "The variable '",
        paste(variable, collapse = ", "),
        "' specified for ",
        variable_label,
        " was not found in the dataset."
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

################################################################################
# WRAP TEXT---------------------------------------------------------------------
################################################################################

wrap_consort_text <- function(
  text,
  width = CONSORT_CONFIG$wrap_width
) {
  stringr::str_wrap(
    as.character(text),
    width = width
  )
}

################################################################################
# UPPER LABEL-------------------------------------------------------------------
################################################################################

make_upper_label <- function(
  label,
  n,
  n_position = CONSORT_CONFIG$upper_n_position
) {
  label <- wrap_consort_text(label)

  n_text <- paste0("N = ", format(n, big.mark = ",", trim = TRUE))

  if (n_position == "before") {
    paste(n_text, label, sep = "\n")
  } else {
    paste(label, n_text, sep = "\n")
  }
}

################################################################################
# LOWER LABEL-------------------------------------------------------------------
################################################################################

make_lower_label <- function(
  label,
  n,
  denominator,
  show_percentage = CONSORT_CONFIG$show_percentages,
  digits = CONSORT_CONFIG$percentage_digits
) {
  label <- wrap_consort_text(label)

  n_text <- format(n, big.mark = ",", trim = TRUE)

  if (
    !isTRUE(show_percentage) ||
      is.null(denominator) ||
      is.na(denominator) ||
      denominator == 0
  ) {
    return(paste(label, n_text, sep = "\n"))
  }

  pct <- round(100 * n / denominator, digits)

  pct_text <- format(pct, nsmall = digits, trim = TRUE)

  paste(label, paste0(n_text, " (", pct_text, "%)"), sep = "\n")
}

################################################################################
# GET REASON COUNTS-------------------------------------------------------------
# Rules:
# 1. Sort descending by count.
# 2. "Other" is always last.
################################################################################

get_reason_counts <- function(
  data,
  flag_variable,
  reason_variable,
  flag_value = "Y"
) {
  if (is.null(flag_variable) || is.null(reason_variable)) {
    return(
      tibble(
        reason = character(),
        n = integer()
      )
    )
  }

  result <- data |>
    filter(
      !is.na(.data[[flag_variable]]),
      .data[[flag_variable]] == flag_value
    ) |>
    mutate(
      reason = as.character(
        .data[[reason_variable]]
      )
    ) |>
    filter(
      !is.na(reason),
      reason != ""
    ) |>
    count(
      reason,
      name = "n"
    ) |>
    arrange(
      desc(n),
      reason
    )

  if (any(result$reason == "Other")) {
    result <- bind_rows(
      result |>
        filter(reason != "Other"),

      result |>
        filter(reason == "Other")
    )
  }

  result
}

################################################################################
# WRAP REASON WITH INTACT COUNT / PERCENTAGE GROUP------------------------------
#
# For long reasons the text will break at wrap_width length, but content within
# (x, x.x%) will not break
#
# Examples:
# Death 4 (6.2%)
# Withdrawal by Subject
# 4 (6.2%)
################################################################################

wrap_reason_with_count <- function(
  reason,
  count_text,
  width = CONSORT_CONFIG$wrap_width
) {
  reason_text <- stringr::str_squish(as.character(reason))

  count_text <- stringr::str_squish(as.character(count_text))

  count_length <- nchar(count_text, type = "width")

  reason_words <- unlist(strsplit(reason_text, "\\s+"))

  reason_words <- reason_words[nzchar(reason_words)]

  lines <- character()
  current_line <- ""

  for (word in reason_words) {
    candidate <- if (nzchar(current_line)) {
      paste(current_line, word)
    } else {
      word
    }

    if (nchar(candidate, type = "width") < width) {
      current_line <- candidate
    } else {
      if (nzchar(current_line)) {
        lines <- c(lines, current_line)
      }
      current_line <- word
    }
  }

  if (nzchar(current_line)) {
    lines <- c(lines, current_line)
  }

  if (length(lines) == 0) {
    return(count_text)
  }

  final_candidate <- paste(lines[length(lines)], count_text)

  if (nchar(final_candidate, type = "width") < width) {
    lines[length(lines)] <- final_candidate
  } else {
    lines <- c(lines, count_text)
  }

  paste(lines, collapse = "\n")
}

################################################################################
# MAKE ONE COMBINED REASON BOX LABEL--------------------------------------------
################################################################################

make_reason_box_label <- function(
  reason_counts,
  denominator = NULL,
  show_percentage = CONSORT_CONFIG$show_percentages,
  digits = CONSORT_CONFIG$percentage_digits
) {
  if (nrow(reason_counts) == 0) {
    return("")
  }

  labels <- character(nrow(reason_counts))

  for (i in seq_len(nrow(reason_counts))) {
    reason <- reason_counts$reason[i]

    n <- reason_counts$n[i]

    if (
      isTRUE(show_percentage) &&
        !is.null(denominator) &&
        !is.na(denominator) &&
        denominator > 0
    ) {
      pct <- round(100 * n / denominator, digits)
      pct_text <- format(pct, nsmall = digits, trim = TRUE)
      count_text <- paste0(format(n, big.mark = ",", trim = TRUE), " (", pct_text, "%)")
    } else {
      count_text <- paste0("(", format(n, big.mark = ",", trim = TRUE), ")")
    }

    labels[i] <- wrap_reason_with_count(reason = reason, count_text = count_text)
  }

  paste(labels, collapse = "\n\n")
}

################################################################################
# REASON BOX HEIGHT-------------------------------------------------------------
################################################################################

get_reason_box_height <- function(
  n_reasons,
  min_height,
  line_spacing,
  padding
) {
  if (n_reasons <= 0) {
    return(0)
  }

  max(min_height, 2 * padding + n_reasons * line_spacing)
}

################################################################################
# GET AGENT NUMBERS-------------------------------------------------------------
################################################################################

get_agent_numbers <- function(
  data
) {
  status_variables <- names(data)[
    stringr::str_detect(
      names(data),
      "^EOTS[0-9]+STT$"
    )
  ]

  if (length(status_variables) == 0) {
    return(integer())
  }

  numbers <- stringr::str_extract(
    status_variables,
    "(?<=EOTS)[0-9]+(?=STT)"
  )

  sort(as.integer(numbers))
}

################################################################################
# GET AGENT DISPLAY NAMES-------------------------------------------------------
################################################################################

get_agent_names <- function(
  config
) {
  if (is.null(config$study_agent_names)) {
    return(character())
  }

  as.character(config$study_agent_names)
}

################################################################################
# MAKE NODE---------------------------------------------------------------------
################################################################################

make_node <- function(
  x,
  y,
  label,
  width,
  height,
  fill = "white",
  colour = "black",
  linewidth = 0.5,
  text_size = 4,
  lineheight = 0.9
) {
  tibble(
    x = x,
    y = y,
    label = label,
    width = width,
    height = height,
    fill = fill,
    colour = colour,
    linewidth = linewidth,
    text_size = text_size,
    lineheight = lineheight
  )
}

################################################################################
# CENTERED POSITIONS
################################################################################

get_centered_positions <- function(
  n,
  center = 0,
  spacing = 3
) {
  if (n <= 0) {
    return(numeric())
  }

  if (n == 1) {
    return(center)
  }

  seq(
    from = center - spacing * (n - 1) / 2,
    to = center + spacing * (n - 1) / 2,
    length.out = n
  )
}

################################################################################
# TREATMENT POSITIONS-----------------------------------------------------------
################################################################################

get_treatment_positions <- function(
  n,
  center,
  spacing
) {
  if (n <= 0) {
    return(numeric())
  }

  if (n == 1) {
    return(center)
  }

  center - spacing * rev(seq_len(n) - 1)
}

################################################################################
# ORTHOGONAL EDGE
################################################################################

make_orthogonal_edge <- function(
  x1,
  y1,
  x2,
  y2,
  connector_y = NULL
) {
  if (is.null(connector_y)) {
    connector_y <- (y1 + y2) / 2
  }

  # Same X----

  if (isTRUE(all.equal(x1, x2))) {
    return(
      tibble(
        x = x1,
        y = y1,
        xend = x2,
        yend = y2
      )
    )
  }

  # Different X----

  tibble(
    x = c(x1, x1, x2),
    y = c(y1, connector_y, connector_y),
    xend = c(x1, x2, x2),
    yend = c(connector_y, connector_y, y2)
  )
}

################################################################################
# ADD ORTHOGONAL EDGE-----------------------------------------------------------
################################################################################

add_orthogonal_edge <- function(
  edges,
  x1,
  y1,
  x2,
  y2,
  connector_y = NULL
) {
  bind_rows(
    edges,

    make_orthogonal_edge(
      x1 = x1,
      y1 = y1,
      x2 = x2,
      y2 = y2,
      connector_y = connector_y
    )
  )
}

################################################################################
# CONSORT DIAGRAM FUNCTION------------------------------------------------------
################################################################################

consort_diagram <- function(
  config = CONSORT_CONFIG,
  layout = LAYOUT_CONFIG
) {
  data <- config$data

  ################
  # VALIDATION----
  ################

  # Screened----
  validate_required_variable(
    data = data,
    variable = config$screened_flag,
    variable_label = "Screened flag"
  )

  # Screen Failure----
  validate_optional_variable(
    data = data,
    variable = config$screen_failure_flag,
    show = config$show_screen_failure,
    variable_label = "Screen Failure flag"
  )

  validate_optional_variable(
    data = data,
    variable = config$screen_failure_reason,
    show = config$show_screen_failure &&
      config$show_screen_failure_reasons,
    variable_label = "Screen Failure reason"
  )

  # Re-screen----
  validate_optional_variable(
    data = data,
    variable = config$rescreen_flag,
    show = config$show_rescreen,
    variable_label = "Re-screen flag"
  )

  # Randomized / Enrolled----
  validate_required_variable(
    data = data,
    variable = config$rand_enr_var,
    variable_label = "Randomized / Enrolled variable"
  )

  # Treatment----
  validate_required_variable(
    data = data,
    variable = config$treatment_variable,
    variable_label = "Treatment variable"
  )

  # Treatment discontinuation reasons----
  # (used when combination_trt = FALSE)
  validate_optional_variable(
    data = data,
    variable = config$treatment_discontinued_reason,
    show = !config$combination_trt &&
      config$show_treatment_discontinued_reason,
    variable_label = "Treatment discontinuation reason"
  )

  # Agent Treatment----
  validate_optional_variable(
    data = data,
    variable = config$ongoing_sub_variable,
    show = config$combination_trt &&
      config$show_ongoing_sub,
    variable_label = "Ongoing Treatment variable"
  )

  validate_optional_variable(
    data = data,
    variable = config$compl_sub_variable,
    show = config$combination_trt &&
      config$show_compl_sub,
    variable_label = "Completed Treatment variable"
  )

  validate_optional_variable(
    data = data,
    variable = config$disc_both_trt_variable,
    show = config$combination_trt &&
      config$show_disc_both_trt,
    variable_label = "Discontinued Both Treatments variable"
  )

  validate_optional_variable(
    data = data,
    variable = config$disc_one_trt_variable,
    show = config$combination_trt &&
      config$show_disc_one_trt,
    variable_label = "Discontinued One Treatment variable"
  )

  # Study----
  if (
    config$show_study_ongoing ||
      config$show_study_completed ||
      config$show_study_discontinued ||
      config$show_study_reasons
  ) {
    validate_required_variable(
      data = data,
      variable = "EOSSTT",
      variable_label = "Study status variable EOSSTT"
    )
  }

  # Study reasons----
  validate_optional_variable(
    data = data,
    variable = "DCSREAS",
    show = config$show_study_reasons,
    variable_label = "Study discontinuation reason"
  )

  # agent numbers and names----
  agent_numbers <-
    get_agent_numbers(data)

  agent_names <-
    get_agent_names(config)

  # Check display names----
  if (
    length(agent_numbers) > 0 &&
      length(agent_names) < max(agent_numbers)
  ) {
    stop(
      paste0(
        "study_agent_names contains only ",
        length(agent_names),
        " name(s), but the dataset contains ",
        max(agent_numbers),
        " study agent(s)."
      ),
      call. = FALSE
    )
  }

  # Validate agent variables when needed----
  if (length(agent_numbers) > 0) {
    for (agent_number in agent_numbers) {
      status_variable <- paste0("EOTS", agent_number, "STT")

      validate_required_variable(
        data = data,
        variable = status_variable,
        variable_label = paste0(
          "Study Agent ",
          agent_number,
          " status"
        )
      )

      if (
        config$show_agent_discontinued_reason &&
          config$show_agent_reasons
      ) {
        reason_variable <- paste0("DCTS", agent_number, "RS")

        validate_required_variable(
          data = data,
          variable = reason_variable,
          variable_label = paste0(
            "Study Agent ",
            agent_number,
            " discontinuation reason"
          )
        )
      }
    }
  }

  # screened N----
  screened_n <- sum(data[[config$screened_flag]] == "Y", na.rm = TRUE)

  # randomized / enrolled data----
  rand_enr_data <-
    data[
      !is.na(
        data[[config$rand_enr_var]]
      ) &
        data[[config$rand_enr_var]] == "Y",
      ,
      drop = FALSE
    ]

  rand_enr_n <- nrow(rand_enr_data)

  # treatment levels----
  treatment_levels <- levels(rand_enr_data[[config$treatment_variable]])

  n_treatments <- length(treatment_levels)

  # storage----
  nodes <- tibble()
  edges <- tibble()

  # screened node----
  nodes <- bind_rows(
    nodes,

    make_node(
      x = layout$screened$x,
      y = layout$screened$y,

      label = make_upper_label(
        label = "Screened",
        n = screened_n,
        n_position = config$upper_n_position
      ),

      width = layout$screened$width,
      height = layout$screened$height,
      fill = layout$screened$fill,
      colour = layout$screened$colour,
      linewidth = layout$screened$linewidth,
      text_size = layout$screened$text_size,
      lineheight = layout$screened$lineheight
    )
  )

  # screening items----
  screening_items <- character()

  if (config$show_screen_failure) {
    screening_items <- c(screening_items, "screen_failure")
  }

  if (config$show_rescreen) {
    screening_items <- c(screening_items, "rescreen")
  }

  screening_items <- c(screening_items, "rand_enr")

  screening_x <-
    get_centered_positions(
      n = length(screening_items),
      center = layout$screened$x,
      spacing = layout$screening_split$x_spacing
    )

  # screen failure node----
  if (config$show_screen_failure) {
    screen_failure_x <- screening_x[which(screening_items == "screen_failure")]

    screen_failure_n <- sum(data[[config$screen_failure_flag]] == "Y", na.rm = TRUE)

    nodes <- bind_rows(
      nodes,

      make_node(
        x = screen_failure_x,
        y = layout$screening$y,

        label = make_upper_label(
          label = "Screen Failure",

          n = screen_failure_n,

          n_position = config$upper_n_position
        ),

        width = layout$screening$width,
        height = layout$screening$height,
        fill = layout$screening$fill,
        colour = layout$screening$colour,
        linewidth = layout$screening$linewidth,
        text_size = layout$screening$text_size,
        lineheight = layout$screening$lineheight
      )
    )
  }

  # re-screen node----
  if (config$show_rescreen) {
    rescreen_x <- screening_x[which(screening_items == "rescreen")]

    rescreen_n <- sum(data[[config$rescreen_flag]] == "Y", na.rm = TRUE)

    nodes <- bind_rows(
      nodes,

      make_node(
        x = rescreen_x,
        y = layout$screening$y,

        label = make_upper_label(
          label = "Re-screened",

          n = rescreen_n,

          n_position = config$upper_n_position
        ),

        width = layout$screening$width,
        height = layout$screening$height,
        fill = layout$screening$fill,
        colour = layout$screening$colour,
        linewidth = layout$screening$linewidth,
        text_size = layout$screening$text_size,
        lineheight = layout$screening$lineheight
      )
    )
  }

  # randomized / enrolled node----
  rand_enr_x <- screening_x[which(screening_items == "rand_enr")]

  nodes <- bind_rows(
    nodes,

    make_node(
      x = rand_enr_x,
      y = layout$screening$y,

      label = make_upper_label(
        label = config$rand_enr_lbl,

        n = rand_enr_n,

        n_position = config$upper_n_position
      ),

      width = layout$screening$width,
      height = layout$screening$height,
      fill = layout$screening$fill,
      colour = layout$screening$colour,
      linewidth = layout$screening$linewidth,
      text_size = layout$screening$text_size,
      lineheight = layout$screening$lineheight
    )
  )

  # screened -> screening----
  for (i in seq_along(screening_items)) {
    edges <-
      add_orthogonal_edge(
        edges = edges,
        x1 = layout$screened$x,
        y1 = layout$screened$y - layout$screened$height / 2,
        x2 = screening_x[i],
        y2 = layout$screening$y + layout$screening$height / 2,
        connector_y = layout$screening_split$y
      )
  }

  # screen failure reasons----
  if (
    config$show_screen_failure &&
      config$show_screen_failure_reasons
  ) {
    reason_counts <-
      get_reason_counts(
        data = data,
        flag_variable = config$screen_failure_flag,
        reason_variable = config$screen_failure_reason,
        flag_value = "Y"
      )

    if (nrow(reason_counts) > 0) {
      screen_failure_x <- screening_x[which(screening_items == "screen_failure")]

      reason_label <-
        make_reason_box_label(
          reason_counts = reason_counts,
          denominator = NULL,
          show_percentage = FALSE
        )

      reason_height <-
        get_reason_box_height(
          n_reasons = nrow(reason_counts),

          min_height = layout$screen_failure_reason_box$min_height,

          line_spacing = layout$screen_failure_reason_box$line_spacing,

          padding = layout$screen_failure_reason_box$padding
        )

      reason_y <- layout$screen_failure_reason_box$y

      nodes <- bind_rows(
        nodes,

        make_node(
          x = screen_failure_x,
          y = reason_y,
          label = reason_label,

          width = layout$screen_failure_reason_box$width,

          height = reason_height,

          fill = layout$screen_failure_reason_box$fill,

          colour = layout$screen_failure_reason_box$colour,

          linewidth = layout$screen_failure_reason_box$linewidth,

          text_size = layout$screen_failure_reason_box$text_size,

          lineheight = layout$screen_failure_reason_box$lineheight
        )
      )

      # Screen failure -> one reason box----

      edges <-
        add_orthogonal_edge(
          edges = edges,
          x1 = screen_failure_x,
          y1 = layout$screening$y - layout$screening$height / 2,
          x2 = screen_failure_x,
          y2 = reason_y + reason_height / 2
        )
    }
  }

  # re-screen children----

  if (config$show_rescreen) {
    rescreen_x <- screening_x[which(screening_items == "rescreen")]
    child_items <- character()
    child_labels <- character()
    child_counts <- numeric()

    # re-screen -> screen failure----

    if (config$show_rescreen_screen_failure) {
      child_items <- c(child_items, "screen_failure")
      child_labels <- c(child_labels, "Screen Failure")

      child_n <- sum(
        data[[config$rescreen_flag]] == "Y" &
          data[[config$screen_failure_flag]] == "Y",
        na.rm = TRUE
      )

      child_counts <- c(child_counts, child_n)
    }

    # re-screen -> randomized / enrolled----

    if (config$show_rescreen_rand_enr) {
      child_items <- c(child_items, "rand_enr")

      child_labels <- c(child_labels, config$rand_enr_lbl)

      child_n <- sum(
        data[[config$rescreen_flag]] == "Y" &
          data[[config$rand_enr_var]] == "Y",
        na.rm = TRUE
      )

      child_counts <- c(child_counts, child_n)
    }

    # Create children----

    if (length(child_items) > 0) {
      child_x <-
        get_centered_positions(
          n = length(child_items),
          center = rescreen_x,
          spacing = layout$rescreen_children$x_spacing
        )

      for (i in seq_along(child_items)) {
        nodes <- bind_rows(
          nodes,

          make_node(
            x = child_x[i],
            y = layout$rescreen_children$y,

            label = make_upper_label(
              label = child_labels[i],
              n = child_counts[i],
              n_position = config$upper_n_position
            ),

            width = layout$rescreen_children$width,

            height = layout$rescreen_children$height,

            fill = layout$rescreen_children$fill,

            colour = layout$rescreen_children$colour,

            linewidth = layout$rescreen_children$linewidth,

            text_size = layout$rescreen_children$text_size,

            lineheight = layout$rescreen_children$lineheight
          )
        )

        edges <-
          add_orthogonal_edge(
            edges = edges,
            x1 = rescreen_x,
            y1 = layout$screening$y - layout$screening$height / 2,
            x2 = child_x[i],
            y2 = layout$rescreen_children$y +
              layout$rescreen_children$height / 2,

            connector_y = layout$rescreen_children$connector_y
          )
      }

      # Horizontal connector to Enrolled/Randomized vertical connector----
      #
      # When the Re-screened Randomized / Enrolled child is
      # displayed, a horizontal connector from the right
      # edge of that child box to the vertical line descending
      # from the main Randomized / Enrolled box is added.

      if ("rand_enr" %in% child_items) {
        rescreen_rand_enr_index <- which(child_items == "rand_enr")[1]

        rescreen_rand_enr_x <- child_x[rescreen_rand_enr_index]

        rescreen_rand_enr_y <- layout$rescreen_children$y

        edges <-
          bind_rows(
            edges,

            tibble(
              x = rescreen_rand_enr_x + layout$rescreen_children$width / 2,
              y = rescreen_rand_enr_y,
              xend = rand_enr_x,
              yend = rescreen_rand_enr_y
            )
          )
      }
    }
  }

  # treatment positions----
  if (
    config$combination_trt &&
      n_treatments > 1
  ) {
    n_lower_columns_for_spacing <- 1 + length(agent_numbers) + 1

    lower_block_width_for_spacing <-
      layout$subject$width / 2 + (n_lower_columns_for_spacing - 1) * layout$columns$x_spacing + layout$study$width / 2

    treatment_spacing <-
      max(
        layout$treatment_split$x_spacing,
        lower_block_width_for_spacing + 0.5
      )
  } else {
    treatment_spacing <- layout$treatment_split$x_spacing
  }

  treatment_x <-
    get_treatment_positions(
      n = n_treatments,
      center = rand_enr_x,
      spacing = treatment_spacing
    )

  # treatment boxes----
  if (n_treatments > 0) {
    for (i in seq_len(n_treatments)) {
      trt <- treatment_levels[i]

      trt_data <-
        rand_enr_data[
          as.character(
            rand_enr_data[[config$treatment_variable]]
          ) ==
            trt,
          ,
          drop = FALSE
        ]

      trt_n <- nrow(trt_data)

      nodes <- bind_rows(
        nodes,

        make_node(
          x = treatment_x[i],
          y = layout$treatment$y,

          label = make_upper_label(
            label = trt,

            n = trt_n,

            n_position = config$upper_n_position
          ),

          width = layout$treatment$width,
          height = layout$treatment$height,
          fill = layout$treatment$fill,
          colour = layout$treatment$colour,
          linewidth = layout$treatment$linewidth,
          text_size = layout$treatment$text_size,
          lineheight = layout$treatment$lineheight
        )
      )
    }

    # Randomized / Enrolled -> Treatment----

    for (i in seq_len(n_treatments)) {
      edges <-
        add_orthogonal_edge(
          edges = edges,

          x1 = rand_enr_x,
          y1 = layout$screening$y - layout$screening$height / 2,
          x2 = treatment_x[i],
          y2 = layout$treatment$y + layout$treatment$height / 2,

          connector_y = layout$treatment_split$connector_y
        )
    }
  }

  # combination treatment mode----
  if (
    config$combination_trt &&
      n_treatments > 0
  ) {
    n_lower_columns <- 1 + length(agent_numbers) + 1

    # Column positions within a treatment bloc

    column_x_relative <-
      get_centered_positions(
        n = n_lower_columns,
        center = 0,
        spacing = layout$columns$x_spacing
      )

    # The lower block must be wide enough for all of its columns
    lower_block_width <-
      layout$subject$width / 2 + (n_lower_columns - 1) * layout$columns$x_spacing + layout$study$width / 2

    # Keep the user-configured spacing, but never allow adjacent
    # lower treatment blocks to overlap
    block_spacing <-
      max(
        layout$lower_blocks$x_spacing,
        lower_block_width + 0.5
      )

    # Treatment boxes are centered over their lower blocks.
    treatment_x <-
      get_treatment_positions(
        n = n_treatments,
        center = rand_enr_x,
        spacing = block_spacing
      )

    # First lower-row position----

    first_box_height <-
      max(
        layout$subject$height,
        layout$agent$height,
        layout$study$height
      )

    lower_start_y <-
      layout$treatment$y -
      layout$treatment$height / 2 -
      layout$lower_blocks$vertical_distance -
      first_box_height / 2

    # Global status-row positions keep the same status aligned
    # horizontally across Treatment, Agents and Study
    lower_status_rows <-
      c(
        "ongoing",
        "completed",
        "discontinued"
      )

    lower_y_spacing <-
      max(
        layout$subject$y_spacing,
        layout$agent$y_spacing,
        layout$study$y_spacing
      )

    lower_status_y <-
      setNames(
        lower_start_y -
          (seq_along(lower_status_rows) - 1) * lower_y_spacing,
        lower_status_rows
      )

    # Helper to connect status boxes vertically inside one column.
    connect_statuses <- function(
      x,
      active_statuses,
      box_height
    ) {
      active_statuses <-
        active_statuses[
          active_statuses %in% names(lower_status_y)
        ]

      if (length(active_statuses) < 2) {
        return(tibble())
      }

      bind_rows(
        lapply(
          seq_len(length(active_statuses) - 1),
          function(k) {
            tibble(
              x = x,
              y = lower_status_y[active_statuses[k]] - box_height / 2,
              xend = x,
              yend = lower_status_y[active_statuses[k + 1]] + box_height / 2
            )
          }
        )
      )
    }

    # loop through treatments----

    for (trt_index in seq_len(n_treatments)) {
      trt <- treatment_levels[trt_index]

      trt_data <-
        rand_enr_data[
          as.character(
            rand_enr_data[[config$treatment_variable]]
          ) ==
            trt,
          ,
          drop = FALSE
        ]

      trt_n <- nrow(trt_data)
      block_center <- treatment_x[trt_index]
      column_x <- column_x_relative + block_center
      subject_x <- column_x[1]
      study_x <- column_x[length(column_x)]

      # agent treatment column----

      subject_nodes <- list()
      subject_active_statuses <- character()

      if (isTRUE(config$show_ongoing_sub)) {
        n <- sum(trt_data[[config$ongoing_sub_variable]] == "Y", na.rm = TRUE)

        subject_nodes$ongoing <-
          list(
            y = lower_status_y["ongoing"],
            label = make_lower_label(
              label = "Ongoing treatment",
              n = n,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )
          )

        subject_active_statuses <- c(subject_active_statuses, "ongoing")
      }

      if (isTRUE(config$show_compl_sub)) {
        n <- sum(trt_data[[config$compl_sub_variable]] == "Y", na.rm = TRUE)

        subject_nodes$completed <-
          list(
            y = lower_status_y["completed"],
            label = make_lower_label(
              label = "Completed treatment",
              n = n,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )
          )

        subject_active_statuses <- c(subject_active_statuses, "completed")
      }

      # The two discontinuation subject/treatment measures remain
      # independent optional components and therefore occupy their
      # own rows after the ongoing/completed rows.
      subject_extra_rows <- character()

      if (isTRUE(config$show_disc_both_trt)) {
        subject_extra_rows <-
          c(subject_extra_rows, "disc_both")
      }

      if (isTRUE(config$show_disc_one_trt)) {
        subject_extra_rows <-
          c(subject_extra_rows, "disc_one")
      }

      subject_next_row <-
        sum(c(isTRUE(config$show_ongoing_sub), isTRUE(config$show_compl_sub))) + 1

      if (length(subject_extra_rows) > 0) {
        for (j in seq_along(subject_extra_rows)) {
          row_type <- subject_extra_rows[j]

          subject_y <-
            lower_start_y -
            (subject_next_row + j - 2) * layout$subject$y_spacing

          if (row_type == "disc_both") {
            n <- sum(trt_data[[config$disc_both_trt_variable]] == "Y", na.rm = TRUE)

            label <-
              make_lower_label(
                label = paste("Discontinued", count_label, "treatments"),
                n = n,
                denominator = trt_n,
                show_percentage = config$show_percentages,
                digits = config$percentage_digits
              )
          } else {
            n <- sum(trt_data[[config$disc_one_trt_variable]] == "Y", na.rm = TRUE)

            label <-
              make_lower_label(
                label = "Discontinued one treatment",
                n = n,
                denominator = trt_n,
                show_percentage = config$show_percentages,
                digits = config$percentage_digits
              )
          }

          subject_nodes[[row_type]] <- list(y = subject_y, label = label)
        }

        # Treat the final subject/treatment row as the discontinued
        # branch for connector purposes
        subject_active_statuses <- c(subject_active_statuses, "discontinued")
      }

      for (node_key in names(subject_nodes)) {
        node_info <- subject_nodes[[node_key]]

        nodes <-
          bind_rows(
            nodes,
            make_node(
              x = subject_x,
              y = node_info$y,
              label = node_info$label,
              width = layout$subject$width,
              height = layout$subject$height,
              fill = layout$subject$fill,
              colour = layout$subject$colour,
              linewidth = layout$subject$linewidth,
              text_size = layout$subject$text_size,
              lineheight = layout$subject$lineheight
            )
          )
      }

      # study agent columns----

      agent_statuses_by_column <- list()

      if (length(agent_numbers) > 0) {
        for (agent_index in seq_along(agent_numbers)) {
          agent_number <- agent_numbers[agent_index]
          agent_x <- column_x[agent_index + 1]
          agent_name <- agent_names[agent_number]
          status_variable <- paste0("EOTS", agent_number, "STT")
          agent_nodes <- list()
          agent_active_statuses <- character()

          if (isTRUE(config$show_agent_ongoing)) {
            n <- sum(trt_data[[status_variable]] == "ONGOING", na.rm = TRUE)

            agent_nodes$ongoing <-
              list(
                y = lower_status_y["ongoing"],
                label = make_lower_label(
                  label = paste("Ongoing", agent_name),
                  n = n,
                  denominator = trt_n,
                  show_percentage = config$show_percentages,
                  digits = config$percentage_digits
                )
              )

            agent_active_statuses <- c(agent_active_statuses, "ongoing")
          }

          if (isTRUE(config$show_agent_completed)) {
            n <- sum(trt_data[[status_variable]] == "COMPLETED", na.rm = TRUE)

            agent_nodes$completed <-
              list(
                y = lower_status_y["completed"],
                label = make_lower_label(
                  label = paste("Completed", agent_name),
                  n = n,
                  denominator = trt_n,
                  show_percentage = config$show_percentages,
                  digits = config$percentage_digits
                )
              )

            agent_active_statuses <- c(agent_active_statuses, "completed")
          }

          discontinued_agent_y <- NULL

          if (isTRUE(config$show_agent_discontinued)) {
            n <- sum(trt_data[[status_variable]] == "DISCONTINUED", na.rm = TRUE)

            discontinued_agent_y <- lower_status_y["discontinued"]

            agent_nodes$discontinued <-
              list(
                y = discontinued_agent_y,
                label = make_lower_label(
                  label = paste("Discontinued", agent_name),
                  n = n,
                  denominator = trt_n,
                  show_percentage = config$show_percentages,
                  digits = config$percentage_digits
                )
              )

            agent_active_statuses <- c(agent_active_statuses, "discontinued")
          }

          for (node_key in names(agent_nodes)) {
            node_info <- agent_nodes[[node_key]]

            nodes <-
              bind_rows(
                nodes,
                make_node(
                  x = agent_x,
                  y = node_info$y,
                  label = node_info$label,
                  width = layout$agent$width,
                  height = layout$agent$height,
                  fill = layout$agent$fill,
                  colour = layout$agent$colour,
                  linewidth = layout$agent$linewidth,
                  text_size = layout$agent$text_size,
                  lineheight = layout$agent$lineheight
                )
              )
          }

          # Agent discontinuation reasons: one combined box----

          if (
            isTRUE(config$show_agent_discontinued) &&
              isTRUE(config$show_agent_discontinued_reason) &&
              isTRUE(config$show_agent_reasons)
          ) {
            reason_variable <- paste0("DCTS", agent_number, "RS")

            reason_counts <-
              get_reason_counts(
                data = trt_data,
                flag_variable = status_variable,
                reason_variable = reason_variable,
                flag_value = "DISCONTINUED"
              )

            if (nrow(reason_counts) > 0) {
              reason_label <-
                make_reason_box_label(
                  reason_counts = reason_counts,
                  denominator = trt_n,
                  show_percentage = config$show_percentages,
                  digits = config$percentage_digits
                )

              reason_height <-
                get_reason_box_height(
                  n_reasons = nrow(reason_counts),
                  min_height = layout$reason$min_height,
                  line_spacing = layout$reason$line_spacing,
                  padding = layout$reason$padding
                )

              reason_y <-
                discontinued_agent_y -
                layout$agent$height / 2 -
                layout$reason$gap -
                reason_height / 2

              nodes <-
                bind_rows(
                  nodes,
                  make_node(
                    x = agent_x,
                    y = reason_y,
                    label = reason_label,
                    width = layout$reason$width,
                    height = reason_height,
                    fill = layout$reason$fill,
                    colour = layout$reason$colour,
                    linewidth = layout$reason$linewidth,
                    text_size = layout$reason$text_size,
                    lineheight = layout$reason$lineheight
                  )
                )

              edges <-
                add_orthogonal_edge(
                  edges = edges,
                  x1 = agent_x,
                  y1 = discontinued_agent_y - layout$agent$height / 2,
                  x2 = agent_x,
                  y2 = reason_y + reason_height / 2
                )
            }
          }

          agent_statuses_by_column[[as.character(agent_number)]] <-
            agent_active_statuses
        }
      }

      # study column----

      study_nodes <- list()
      study_active_statuses <- character()

      if (isTRUE(config$show_study_ongoing)) {
        study_n <- sum(trt_data$EOSSTT == "ONGOING", na.rm = TRUE)

        study_nodes$ongoing <-
          list(
            y = lower_status_y["ongoing"],
            label = make_lower_label(
              label = "Ongoing Study",
              n = study_n,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )
          )

        study_active_statuses <- c(study_active_statuses, "ongoing")
      }

      if (isTRUE(config$show_study_completed)) {
        study_n <- sum(trt_data$EOSSTT == "COMPLETED", na.rm = TRUE)

        study_nodes$completed <-
          list(
            y = lower_status_y["completed"],
            label = make_lower_label(
              label = "Completed Study",
              n = study_n,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )
          )

        study_active_statuses <- c(study_active_statuses, "completed")
      }

      discontinued_study_y <- NULL

      if (isTRUE(config$show_study_discontinued)) {
        study_n <- sum(trt_data$EOSSTT == "DISCONTINUED", na.rm = TRUE)

        discontinued_study_y <- lower_status_y["discontinued"]

        study_nodes$discontinued <-
          list(
            y = discontinued_study_y,
            label = make_lower_label(
              label = "Discontinued Study",
              n = study_n,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )
          )

        study_active_statuses <- c(study_active_statuses, "discontinued")
      }

      for (node_key in names(study_nodes)) {
        node_info <- study_nodes[[node_key]]

        nodes <-
          bind_rows(
            nodes,
            make_node(
              x = study_x,
              y = node_info$y,
              label = node_info$label,
              width = layout$study$width,
              height = layout$study$height,
              fill = layout$study$fill,
              colour = layout$study$colour,
              linewidth = layout$study$linewidth,
              text_size = layout$study$text_size,
              lineheight = layout$study$lineheight
            )
          )
      }

      # Study discontinuation reasons: one combined box----

      if (
        isTRUE(config$show_study_discontinued) &&
          isTRUE(config$show_study_reasons)
      ) {
        reason_counts <-
          get_reason_counts(
            data = trt_data,
            flag_variable = "EOSSTT",
            reason_variable = "DCSREAS",
            flag_value = "DISCONTINUED"
          )

        if (nrow(reason_counts) > 0) {
          reason_label <-
            make_reason_box_label(
              reason_counts = reason_counts,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )

          reason_height <-
            get_reason_box_height(
              n_reasons = nrow(reason_counts),
              min_height = layout$reason$min_height,
              line_spacing = layout$reason$line_spacing,
              padding = layout$reason$padding
            )

          reason_y <-
            discontinued_study_y -
            layout$study$height / 2 -
            layout$reason$gap -
            reason_height / 2

          nodes <-
            bind_rows(
              nodes,
              make_node(
                x = study_x,
                y = reason_y,
                label = reason_label,
                width = layout$reason$width,
                height = reason_height,
                fill = layout$reason$fill,
                colour = layout$reason$colour,
                linewidth = layout$reason$linewidth,
                text_size = layout$reason$text_size,
                lineheight = layout$reason$lineheight
              )
            )

          edges <-
            add_orthogonal_edge(
              edges = edges,
              x1 = study_x,
              y1 = discontinued_study_y - layout$study$height / 2,
              x2 = study_x,
              y2 = reason_y + reason_height / 2
            )
        }
      }

      # treatment -> lower block----

      lower_connector_y <- layout$lower_blocks$connector_y

      # Lower-column spine positions----

      subject_spine_x <-
        subject_x -
        layout$subject$width / 2 -
        layout$lower_blocks$spine_gap

      agent_spine_x <-
        if (length(agent_numbers) > 0) {
          column_x[2:(length(agent_numbers) + 1)] -
            layout$agent$width / 2 -
            layout$lower_blocks$spine_gap
        } else {
          numeric()
        }

      study_spine_x <- study_x - layout$study$width / 2 - layout$lower_blocks$spine_gap

      lower_spine_x <-
        c(
          if (length(subject_nodes) > 0) subject_spine_x,
          agent_spine_x,
          study_spine_x
        )

      lower_spine_x <- sort(unique(lower_spine_x))

      # Treatment -> horizontal lower trunk----

      edges <-
        add_orthogonal_edge(
          edges = edges,
          x1 = block_center,
          y1 = layout$treatment$y - layout$treatment$height / 2,
          x2 = block_center,
          y2 = lower_connector_y
        )

      # The horizontal trunk runs between the LOWER COLUMN SPINES,
      # not between the centers of the boxes.
      if (length(lower_spine_x) > 1) {
        edges <-
          bind_rows(
            edges,
            tibble(
              x = min(lower_spine_x),
              y = lower_connector_y,
              xend = max(lower_spine_x),
              yend = lower_connector_y
            )
          )
      }

      # Helper: create one lower-column spine and its branches----

      add_lower_column_spine <- function(
        edges,
        spine_x,
        box_x,
        box_width,
        status_nodes,
        status_order
      ) {
        if (length(status_nodes) == 0) {
          return(edges)
        }

        status_order <-
          status_order[
            status_order %in% names(status_nodes)
          ]

        if (length(status_order) == 0) {
          return(edges)
        }

        status_y <-
          vapply(
            status_nodes[status_order],
            function(node) node$y,
            numeric(1)
          )

        # The main spine starts at the lower horizontal trunk and
        # ENDS at the center of the LAST status branch. It never
        # continues below the Discontinued branch.
        spine_end_y <-
          status_y[length(status_y)]

        edges <-
          bind_rows(
            edges,
            tibble(
              x = spine_x,
              y = lower_connector_y,
              xend = spine_x,
              yend = spine_end_y
            )
          )

        # Horizontal branches from the spine to the LEFT edge of
        # each status box.
        for (status_name in status_order) {
          y_value <-
            status_nodes[[status_name]]$y

          edges <-
            bind_rows(
              edges,
              tibble(
                x = spine_x,
                y = y_value,
                xend = box_x - box_width / 2,
                yend = y_value
              )
            )
        }

        edges
      }

      # Agent treatment spine
      subject_status_order <-
        c(
          "ongoing",
          "completed",
          "disc_both",
          "disc_one"
        )

      edges <-
        add_lower_column_spine(
          edges = edges,
          spine_x = subject_spine_x,
          box_x = subject_x,
          box_width = layout$subject$width,
          status_nodes = subject_nodes,
          status_order = subject_status_order
        )

      # study agent spines
      if (length(agent_numbers) > 0) {
        for (agent_index in seq_along(agent_numbers)) {
          agent_number <- agent_numbers[agent_index]

          agent_x <- column_x[agent_index + 1]

          # Rebuild the status nodes for connector geometry from
          # the same aligned rows used to create the boxes.
          agent_status_nodes <- list()

          status_variable <- paste0("EOTS", agent_number, "STT")

          if (isTRUE(config$show_agent_ongoing)) {
            agent_status_nodes$ongoing <- list(y = lower_status_y["ongoing"])
          }

          if (isTRUE(config$show_agent_completed)) {
            agent_status_nodes$completed <- list(y = lower_status_y["completed"])
          }

          if (isTRUE(config$show_agent_discontinued)) {
            agent_status_nodes$discontinued <- list(y = lower_status_y["discontinued"])
          }

          edges <-
            add_lower_column_spine(
              edges = edges,
              spine_x = agent_spine_x[agent_index],
              box_x = agent_x,
              box_width = layout$agent$width,
              status_nodes = agent_status_nodes,
              status_order = c(
                "ongoing",
                "completed",
                "discontinued"
              )
            )
        }
      }

      # study spine----

      study_status_nodes <- list()

      if (isTRUE(config$show_study_ongoing)) {
        study_status_nodes$ongoing <- list(y = lower_status_y["ongoing"])
      }

      if (isTRUE(config$show_study_completed)) {
        study_status_nodes$completed <- list(y = lower_status_y["completed"])
      }

      if (isTRUE(config$show_study_discontinued)) {
        study_status_nodes$discontinued <- list(y = lower_status_y["discontinued"])
      }

      edges <-
        add_lower_column_spine(
          edges = edges,
          spine_x = study_spine_x,
          box_x = study_x,
          box_width = layout$study$width,
          status_nodes = study_status_nodes,
          status_order = c(
            "ongoing",
            "completed",
            "discontinued"
          )
        )
    }
  }

  # non-combination treatment mode----

  if (!config$combination_trt) {
    # Helper for non-combination lower-column spines

    add_lower_column_spine_noncombination <- function(
      edges,
      spine_x,
      box_x,
      box_width,
      status_nodes,
      status_order,
      connector_y
    ) {
      if (length(status_nodes) == 0) {
        return(edges)
      }

      status_order <- status_order[status_order %in% names(status_nodes)]

      if (length(status_order) == 0) {
        return(edges)
      }

      status_y <-
        vapply(
          status_nodes[status_order],
          function(node) node$y,
          numeric(1)
        )

      spine_end_y <- status_y[length(status_y)]

      edges <-
        bind_rows(
          edges,
          tibble(
            x = spine_x,
            y = connector_y,
            xend = spine_x,
            yend = spine_end_y
          )
        )

      for (status_name in status_order) {
        y_value <- status_nodes[[status_name]]$y

        edges <-
          bind_rows(
            edges,
            tibble(
              x = spine_x,
              y = y_value,
              xend = box_x - box_width / 2,
              yend = y_value
            )
          )
      }

      edges
    }

    # Lower-row geometry----

    first_box_height <-
      max(
        layout$treatment_status$height,
        layout$study$height
      )

    lower_start_y <-
      layout$treatment$y -
      layout$treatment$height / 2 -
      layout$lower_blocks$vertical_distance -
      first_box_height / 2

    lower_status_rows <-
      c(
        "ongoing",
        "completed",
        "discontinued"
      )

    lower_y_spacing <-
      max(
        layout$treatment_status$y_spacing,
        layout$study$y_spacing
      )

    lower_status_y <-
      setNames(
        lower_start_y -
          (seq_along(lower_status_rows) - 1) * lower_y_spacing,
        lower_status_rows
      )

    # The lower horizontal trunk sits above the first status row.
    # The lower-column spine helper uses this value for its
    # vertical connection down to the status branches.
    lower_connector_y <- layout$lower_blocks$connector_y

    # One lower Treatment and Study block per treatment arm

    for (trt_index in seq_len(n_treatments)) {
      trt <- treatment_levels[trt_index]

      trt_data <-
        rand_enr_data[
          as.character(
            rand_enr_data[[config$treatment_variable]]
          ) ==
            trt,
          ,
          drop = FALSE
        ]

      trt_n <- nrow(trt_data)
      block_center <- treatment_x[trt_index]

      # The Treatment and Study lower columns are centered as a pair
      # beneath the corresponding Treatment arm
      treatment_status_x <-
        block_center -
        layout$columns$x_spacing / 2

      study_x <-
        block_center +
        layout$columns$x_spacing / 2

      # treatment status----

      treatment_status_nodes <- list()
      discontinued_treatment_y <- NULL

      treatment_status_rows <-
        c(
          "ongoing",
          "completed",
          "discontinued"
        )

      for (status_type in treatment_status_rows) {
        if (status_type == "ongoing") {
          treatment_n <- sum(trt_data$EOSSTT == "ONGOING", na.rm = TRUE)

          label <-
            make_lower_label(
              label = "Ongoing Treatment",
              n = treatment_n,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )
        } else if (status_type == "completed") {
          treatment_n <- sum(trt_data$EOSSTT == "COMPLETED", na.rm = TRUE)

          label <-
            make_lower_label(
              label = "Completed Treatment",
              n = treatment_n,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )
        } else {
          treatment_n <- sum(trt_data$EOSSTT == "DISCONTINUED", na.rm = TRUE)

          label <-
            make_lower_label(
              label = "Discontinued Treatment",
              n = treatment_n,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )

          discontinued_treatment_y <- lower_status_y["discontinued"]
        }

        treatment_status_nodes[[status_type]] <-
          list(
            y = lower_status_y[status_type],
            label = label
          )
      }

      for (node_key in names(treatment_status_nodes)) {
        node_info <- treatment_status_nodes[[node_key]]

        nodes <-
          bind_rows(
            nodes,
            make_node(
              x = treatment_status_x,
              y = node_info$y,
              label = node_info$label,
              width = layout$treatment_status$width,
              height = layout$treatment_status$height,
              fill = layout$treatment_status$fill,
              colour = layout$treatment_status$colour,
              linewidth = layout$treatment_status$linewidth,
              text_size = layout$treatment_status$text_size,
              lineheight = layout$treatment_status$lineheight
            )
          )
      }

      # Treatment lower-column spine.
      treatment_spine_x <-
        treatment_status_x -
        layout$treatment_status$width / 2 -
        layout$lower_blocks$spine_gap

      edges <-
        add_lower_column_spine_noncombination(
          edges = edges,
          spine_x = treatment_spine_x,
          box_x = treatment_status_x,
          box_width = layout$treatment_status$width,
          status_nodes = treatment_status_nodes,
          status_order = treatment_status_rows,
          connector_y = lower_connector_y
        )

      # Treatment discontinuation reasons----

      if (
        !is.null(discontinued_treatment_y) &&
          isTRUE(config$show_treatment_discontinued_reason)
      ) {
        reason_counts <-
          get_reason_counts(
            data = trt_data,
            flag_variable = "EOTSTT",
            reason_variable = config$treatment_discontinued_reason,
            flag_value = "DISCONTINUED"
          )

        if (nrow(reason_counts) > 0) {
          reason_label <-
            make_reason_box_label(
              reason_counts = reason_counts,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )

          reason_height <-
            get_reason_box_height(
              n_reasons = nrow(reason_counts),
              min_height = layout$reason$min_height,
              line_spacing = layout$reason$line_spacing,
              padding = layout$reason$padding
            )

          reason_y <-
            discontinued_treatment_y -
            layout$treatment_status$height / 2 -
            layout$reason$gap -
            reason_height / 2

          nodes <-
            bind_rows(
              nodes,
              make_node(
                x = treatment_status_x,
                y = reason_y,
                label = reason_label,
                width = layout$reason$width,
                height = reason_height,
                fill = layout$reason$fill,
                colour = layout$reason$colour,
                linewidth = layout$reason$linewidth,
                text_size = layout$reason$text_size,
                lineheight = layout$reason$lineheight
              )
            )

          edges <-
            add_orthogonal_edge(
              edges = edges,
              x1 = treatment_status_x,
              y1 = discontinued_treatment_y -
                layout$treatment_status$height / 2,
              x2 = treatment_status_x,
              y2 = reason_y + reason_height / 2
            )
        }
      }

      # study status----
      study_status_nodes <- list()

      if (isTRUE(config$show_study_ongoing)) {
        study_n <- sum(trt_data$EOSSTT == "ONGOING", na.rm = TRUE)

        study_status_nodes$ongoing <-
          list(
            y = lower_status_y["ongoing"],
            label = make_lower_label(
              label = "Ongoing Study",
              n = study_n,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )
          )
      }

      if (isTRUE(config$show_study_completed)) {
        study_n <- sum(trt_data$EOSSTT == "COMPLETED", na.rm = TRUE)

        study_status_nodes$completed <-
          list(
            y = lower_status_y["completed"],
            label = make_lower_label(
              label = "Completed Study",
              n = study_n,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )
          )
      }

      discontinued_study_y <- NULL

      if (isTRUE(config$show_study_discontinued)) {
        study_n <- sum(trt_data$EOSSTT == "DISCONTINUED", na.rm = TRUE)

        discontinued_study_y <- lower_status_y["discontinued"]

        study_status_nodes$discontinued <-
          list(
            y = discontinued_study_y,
            label = make_lower_label(
              label = "Discontinued Study",
              n = study_n,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )
          )
      }

      for (node_key in names(study_status_nodes)) {
        node_info <- study_status_nodes[[node_key]]

        nodes <-
          bind_rows(
            nodes,
            make_node(
              x = study_x,
              y = node_info$y,
              label = node_info$label,
              width = layout$study$width,
              height = layout$study$height,
              fill = layout$study$fill,
              colour = layout$study$colour,
              linewidth = layout$study$linewidth,
              text_size = layout$study$text_size,
              lineheight = layout$study$lineheight
            )
          )
      }

      # Study lower-column spine.
      study_spine_x <-
        study_x -
        layout$study$width / 2 -
        layout$lower_blocks$spine_gap -
        0.05

      edges <-
        add_lower_column_spine_noncombination(
          edges = edges,
          spine_x = study_spine_x,
          box_x = study_x,
          box_width = layout$study$width,
          status_nodes = study_status_nodes,
          status_order = c(
            "ongoing",
            "completed",
            "discontinued"
          ),
          connector_y = lower_connector_y
        )

      # Study discontinuation reasons----

      if (
        !is.null(discontinued_study_y) &&
          isTRUE(config$show_study_reasons)
      ) {
        reason_counts <-
          get_reason_counts(
            data = trt_data,
            flag_variable = "EOSSTT",
            reason_variable = "DCSREAS",
            flag_value = "DISCONTINUED"
          )

        if (nrow(reason_counts) > 0) {
          reason_label <-
            make_reason_box_label(
              reason_counts = reason_counts,
              denominator = trt_n,
              show_percentage = config$show_percentages,
              digits = config$percentage_digits
            )

          reason_height <-
            get_reason_box_height(
              n_reasons = nrow(reason_counts),
              min_height = layout$reason$min_height,
              line_spacing = layout$reason$line_spacing,
              padding = layout$reason$padding
            )

          reason_y <-
            discontinued_study_y -
            layout$study$height / 2 -
            layout$reason$gap -
            reason_height / 2

          nodes <-
            bind_rows(
              nodes,
              make_node(
                x = study_x,
                y = reason_y,
                label = reason_label,
                width = layout$reason$width,
                height = reason_height,
                fill = layout$reason$fill,
                colour = layout$reason$colour,
                linewidth = layout$reason$linewidth,
                text_size = layout$reason$text_size,
                lineheight = layout$reason$lineheight
              )
            )

          edges <-
            add_orthogonal_edge(
              edges = edges,
              x1 = study_x,
              y1 = discontinued_study_y -
                layout$study$height / 2,
              x2 = study_x,
              y2 = reason_y + reason_height / 2
            )
        }
      }

      # treatment -> lower section----
      edges <-
        bind_rows(
          edges,
          tibble(
            x = treatment_x[trt_index],
            y = layout$treatment$y -
              layout$treatment$height / 2,
            xend = treatment_x[trt_index],
            yend = lower_connector_y - 0.02
          ),
          tibble(
            x = treatment_spine_x,
            y = lower_connector_y,
            xend = study_spine_x,
            yend = lower_connector_y
          )
        )
    }
  }

  # PLOT----
  p <- ggplot() +

    geom_segment(
      data = edges,

      aes(x = x, y = y, xend = xend, yend = yend),

      linewidth = layout$lines$linewidth,
      lineend = layout$lines$lineend
    ) +

    geom_rect(
      data = nodes,

      aes(
        xmin = x - width / 2,
        xmax = x + width / 2,
        ymin = y - height / 2,
        ymax = y + height / 2,

        fill = fill,
        colour = colour
      ),

      linewidth = 0.5
    ) +

    geom_text(
      data = nodes,

      aes(
        x = x,
        y = y,
        label = label,
        size = text_size,

        lineheight = lineheight
      ),

      show.legend = FALSE
    ) +

    scale_fill_identity() +
    scale_colour_identity() +
    scale_size_identity() +

    coord_equal(
      expand = TRUE,
      clip = "off"
    ) +

    theme_void() +

    theme(
      plot.margin = grid::unit(
        layout$plot$margin,
        "pt"
      )
    )

  list(
    plot = p,
    nodes = nodes,
    edges = edges,
    treatment_levels = treatment_levels,
    agent_numbers = agent_numbers,
    agent_names = agent_names
  )
}


################################################################################
# RUN CONSORT DIAGRAM-----------------------------------------------------------
################################################################################

consort_result <- consort_diagram()

################################################################################
# DISPLAY-----------------------------------------------------------------------
################################################################################

consort_result$plot

# Save as PNG----
pname <- paste0(tolower(tblid), ".png")

ggsave(
  filename = write_path(opath, pname),
  plot = consort_result$plot,
  width = 75,
  height = 40,
  units = "cm",
  dpi = 600,
  device = grDevices::png,
  type = "cairo",
  bg = "white"
)

plot_c <- magick::image_read(write_path(opath, pname))
plot_c <- magick::image_trim(plot_c)
plot_c <- magick::image_border(plot_c, "white", "15x15")
plot_c <- magick::image_resize(plot_c, geometry = "5250x")
magick::image_write(plot_c, write_path(opath, pname))

tidytlg::gentlg(
  tlf = "g",
  plotnames = write_path(opath, pname),
  orientation = "landscape",
  opath = write_path(opath),
  file = tblid,
  title = title_footer$title,
  footers = title_footer$main_footer
)
