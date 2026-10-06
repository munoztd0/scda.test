###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              tsfae01a.r
## R version:                 4.5.2
## junco Version:             0.1.3
## Short Description:         Program to create tsfae01a: Overall Summary of Subjects
##                            With Treatment-emergent Adverse Events
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     adsl, adae, adexsum
## Output:                    tsfae01a.rtf
## Remarks:                   This variant includes summary of AEs by max severity.
##
## Modification History:
##  Rev #:
##  Modified By:
##  Reporting Effort:
##  Date:
##  Description:
################################################################################

################################################################################
# Prep environment:
################################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)
library(dplyr)
library(rtables)
library(junco)
library(stringr)

################################################################################
# Define script level parameters:
################################################################################

tblid <- "TSFAE01a"
fileid <- write_path(opath, tblid)
popfl <- "SAFFL"
trtvar <- "TRT01A"
tab_titles <- get_titles_internal(tblid)
string_map <- make_jj_str_map()
combined_colspan_trt <- TRUE
risk_diff <- TRUE
ctrl_grp <- "Placebo"
durfup_paramcd <- "DURFUPD" # Options: DURFUPD, DURFUPW, DURFUPM, DURFUPY
# Section --- Optional row flags --------------------
show_ae_related_aes_row <- TRUE
show_ae_related_saes_row <- TRUE
show_ae_rel_ae_death_row <- TRUE
show_ae_infections_row <- TRUE
show_ae_serious_infections_row <- TRUE
show_ae_serious_rel_infections_row <- TRUE
show_ae_dosemod_row <- TRUE

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
# Process data:
################################################################################
adsl <- haven::read_sas(envsetup::read_path(a_in, "adsl.sas7bdat")) %>%
  df_na() %>%
  filter(!!rlang::sym(popfl) == "Y") %>%
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c(
        "Xanomeline Low Dose",
        "Xanomeline High Dose",
        "Placebo"
      )
    )
  ) %>%
  create_colspan_var(
    non_active_grp = ctrl_grp,
    non_active_grp_span_lbl = " ",
    active_grp_span_lbl = "Active Study Agent",
    colspan_var = "colspan_trt",
    trt_var = trtvar
  ) %>%
  select(
    USUBJID,
    !!rlang::sym(popfl),
    !!rlang::sym(trtvar),
    colspan_trt
  )

if (risk_diff == TRUE) {
  adsl$rrisk_header <- "Risk Difference (%) (95% CI)"
  adsl$rrisk_label <- paste(adsl[[trtvar]], paste("vs", ctrl_grp))
}

# Section 1 - Read, filter to TEAEs, select needed columns ----------------------
adae <- haven::read_sas(envsetup::read_path(a_in, "adae.sas7bdat")) |>
  df_na() |>
  filter(toupper(TRTEMFL) == "Y") |>
  select(
    USUBJID,
    AESER,
    AESDTH,
    AESLIFE,
    AESHOSP,
    AESDISAB,
    AESCONG,
    AESMIE,
    AEACN,
    AESEV,
    ASEVN,
    AEREL,
    TRDISCFL,
    AEOUT,
    AOCCIFL,
    TRTEMFL,
    AESOC
  )

# Section 2 - Derive maxsev per subject ------------------------------------------
# Pick ASEVN from the first occurrence flag (AOCCIFL=Y), fill missing, factorize
adae <- adae |>
  group_by(USUBJID) |>
  mutate(maxsev = ASEVN[which(AOCCIFL == "Y")][1]) |>
  ungroup() |>
  mutate(maxsev = ifelse(is.na(maxsev), "Missing", maxsev)) |>
  mutate(
    maxsev = factor(
      maxsev,
      levels = c("1", "2", "3", "Missing"),
      labels = c("Mild", "Moderate", "Severe", "Missing")
    )
  )

# Section 3 - Derive related SAE, related AE-death, and infection flags ----------
adae <- adae |>
  mutate(
    Rel_SAEs = if_else(
      toupper(AESER) == "Y" & toupper(AEREL) == "RELATED",
      "Y",
      NA_character_
    ),
    Rel_AE_Death = if_else(
      toupper(AEOUT) == "FATAL" & toupper(AEREL) == "RELATED",
      "Y",
      NA_character_
    ),
    infection_flg = if_else(
      toupper(AESOC) == "INFECTIONS AND INFESTATIONS",
      "Y",
      NA_character_
    ),
    serious_infection_flg = if_else(
      toupper(AESOC) == "INFECTIONS AND INFESTATIONS" &
        toupper(AESER) == "Y",
      "Y",
      NA_character_
    ),
    serious_rel_infection_flg = if_else(
      toupper(AESOC) == "INFECTIONS AND INFESTATIONS" &
        toupper(AESER) == "Y" &
        toupper(AEREL) == "RELATED",
      "Y",
      NA_character_
    )
  )

# Section 4 - Sentence-case AEACN factor levels ----------------------------------
adae <- adae |>
  mutate(
    AEACN = {
      levels(AEACN) <- str_to_sentence(levels(AEACN))
      AEACN
    }
  )

adexsum <- haven::read_sas(envsetup::read_path(a_in, "adexsum.sas7bdat")) |>
  df_na() |>
  filter(!!rlang::sym(popfl) == "Y") |>
  select(
    USUBJID,
    PARAMCD,
    PARAM,
    AVAL
  )

# Section --- Average duration of follow-up filter --------------------
adexsum_durfup <- adexsum |>
  filter(toupper(PARAMCD) == toupper(durfup_paramcd))

durfup_unit <- adexsum_durfup |>
  pull(PARAM) |>
  unique() |>
  str_extract("\\([^)]+\\)")

if (is.na(durfup_unit)) {
  durfup_unit <- "(-)"
  warning("Unit not found in PARAM for selected PARAMCD. Defaulting to '(days)'.")
}

durfup_label <- sprintf(
  "Average duration of follow-up %s",
  durfup_unit
)

adexsum_durfup <- adexsum_durfup |>
  select(USUBJID, AVAL) |>
  rename(durfup_aval = AVAL)

# Section --- Average exposure filter --------------------
adexsum_exposure <- adexsum |>
  filter(toupper(PARAMCD) == "TNUMDOS") |>
  select(USUBJID, AVAL) |>
  rename(exposure_aval = AVAL)

adae <- adsl |>
  left_join(adexsum_durfup, by = "USUBJID") |>
  left_join(adexsum_exposure, by = "USUBJID") |>
  inner_join(adae, by = c("USUBJID"))

# Section --- Deduplicate subject-level metrics for correct mean --------
# durfup_aval and exposure_aval are subject-level; retain only 1 unique
# value per subject to avoid record-weighted mean.
adae <- adae |>
  group_by(USUBJID) |>
  mutate(
    durfup_aval = if_else(
      !duplicated(durfup_aval),
      durfup_aval,
      NA_real_
    ),
    exposure_aval = if_else(
      !duplicated(exposure_aval),
      exposure_aval,
      NA_real_
    )
  ) |>
  ungroup()

################################################################################
# Define layout and build table:
################################################################################
# Section --- Row labels --------------------
rel_ae_label <- "Related AEs"
rel_sae_label <- "Related SAEs"
rel_aedth_label <- "Related AEs leading to death"
acn_label <- "AE leading to dose modification of any study treatment~[super a,b]"

colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = ctrl_grp,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)

# Section --- Check the levels of AEACN --------------------
aeacn_levels <- levels(adae$AEACN) %>%
  str_to_sentence() %>%
  unique()

excl_aeacn_levels <- c("Drug withdrawn", "Dose not changed", "Not applicable", "Unknown")
dosemod_lvls <- aeacn_levels[!(aeacn_levels %in% excl_aeacn_levels)]

# Rearrange levels for AEACN
newsort_AEACN <- unique(c(
  "Drug interrupted",
  "Dose reduced",
  "Dose rate reduced",
  "Dose increased",
  aeacn_levels
))

adae <- adae %>%
  mutate(AEACN = forcats::fct_relevel(AEACN, newsort_AEACN))


dosemod_spf <- make_combo_splitfun(
  nm = "modified",
  label = "AE leading to dose modification of study",
  levels = c(
    "Dose reduced",
    "Dose increased",
    "Drug interrupted",
    "Dose rate reduced"
  )
)
aesevall_spf <- make_combo_splitfun(
  nm = "AESEV_ALL",
  label = "Worst severity",
  levels = NULL
)

sae_class <- make_combo_splitfun(
  nm = "SAE_CLASS",
  label = "SAE classification~[super b]",
  levels = "Y"
)


rr_method <- "wald"
ref_path <- c("colspan_trt", " ", trtvar, ctrl_grp)
extra_args_rr <- list(
  method = rr_method,
  ref_path = ref_path,
  .stats = c("count_unique_fraction")
)

lyt <- basic_table(
  show_colcounts = TRUE,
  colcount_format = "N=xx",
  top_level_section_div = " "
) %>%
  append_topleft(c(" ", " ", "Event, n (%)")) %>%
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

lyt <- lyt %>%
  analyze_vars(
    "durfup_aval",
    .stats = "mean",
    .formats = c(mean = "xx.x"),
    .labels = c(mean = durfup_label),
    show_labels = "hidden",
    section_div = " "
  ) %>%
  analyze_vars(
    "exposure_aval",
    .stats = "mean",
    .formats = c(mean = "xx.x"),
    .labels = c(mean = "Average exposure (number of administrations)"),
    show_labels = "hidden"
  ) %>%
  split_rows_by(
    "TRTEMFL",
    split_fun = keep_split_levels("Y"),
    section_div = " "
  ) %>%
  summarize_row_groups(
    "TRTEMFL",
    cfun = a_freq_j,
    extra_args = list(
      label = "AEs",
      method = rr_method,
      ref_path = ref_path,
      .stats = c("count_unique_fraction")
    )
  ) %>%
  analyze(
    "AESER",
    afun = a_freq_j,
    show_labels = "hidden",
    extra_args = append(
      extra_args_rr,
      list(label = "SAEs", val = "Y", NULL)
    )
  )

if (show_ae_related_aes_row) {
  lyt <- lyt %>%
    analyze(
      "AEREL",
      afun = a_freq_j,
      show_labels = "hidden",
      extra_args = append(
        extra_args_rr,
        list(label = rel_ae_label, val = "RELATED", NULL)
      )
    )
}

if (show_ae_related_saes_row) {
  lyt <- lyt %>%
    analyze(
      "Rel_SAEs",
      afun = a_freq_j,
      show_labels = "hidden",
      extra_args = append(
        extra_args_rr,
        list(label = rel_sae_label, val = "Y", NULL)
      )
    )
}

lyt <- lyt %>%
  analyze(
    "TRDISCFL",
    afun = a_freq_j,
    show_labels = "hidden",
    extra_args = append(
      extra_args_rr,
      list(
        label = "AE leading to permanent discontinuation of study treatment",
        val = "Y",
        NULL
      )
    )
  )

if (show_ae_rel_ae_death_row) {
  lyt <- lyt %>%
    analyze(
      "Rel_AE_Death",
      afun = a_freq_j,
      show_labels = "hidden",
      extra_args = append(
        extra_args_rr,
        list(label = rel_aedth_label, val = "Y", NULL)
      )
    )
}

if (show_ae_infections_row) {
  lyt <- lyt %>%
    analyze(
      "infection_flg",
      afun = a_freq_j,
      show_labels = "hidden",
      extra_args = append(
        extra_args_rr,
        list(label = "Infections", val = "Y", NULL)
      )
    )
}

if (show_ae_serious_infections_row) {
  lyt <- lyt %>%
    analyze(
      "serious_infection_flg",
      afun = a_freq_j,
      show_labels = "hidden",
      extra_args = append(
        extra_args_rr,
        list(label = "Serious infections", val = "Y", NULL)
      )
    )
}

if (show_ae_serious_rel_infections_row) {
  lyt <- lyt %>%
    analyze(
      "serious_rel_infection_flg",
      afun = a_freq_j,
      show_labels = "hidden",
      extra_args = append(
        extra_args_rr,
        list(label = "Serious related infections", val = "Y", NULL)
      )
    )
}

lyt <- lyt %>%
  split_rows_by("maxsev", split_fun = aesevall_spf) %>%
  analyze("maxsev", afun = a_freq_j, extra_args = append(extra_args_rr, NULL))


if (show_ae_dosemod_row) {
  lyt <- lyt %>%
    split_rows_by(
      "AEACN",
      split_fun = dosemod_spf,
      section_div = " "
    ) %>%
    summarize_row_groups(
      "AEACN",
      cfun = a_freq_j,
      extra_args = list(
        label = acn_label,
        method = rr_method,
        ref_path = ref_path,
        .stats = c("count_unique_fraction")
      )
    ) %>%
    analyze(
      "AEACN",
      table_names = "AEACN",
      a_freq_j,
      show_labels = "hidden",
      extra_args = append(
        extra_args_rr,
        list(
          excl_levels = excl_aeacn_levels,
          drop_levels = TRUE
        )
      )
    )
}


lyt <- lyt %>%
  split_rows_by(
    "AESER",
    split_fun = sae_class,
    section_div = " ",
  ) %>%
  analyze(
    "AESDTH",
    afun = a_freq_j,
    show_labels = "hidden",
    extra_args = append(
      extra_args_rr,
      list(label = "Death", val = "Y", NULL)
    )
  ) %>%
  analyze(
    "AESLIFE",
    afun = a_freq_j,
    show_labels = "hidden",
    extra_args = append(
      extra_args_rr,
      list(label = "Life-threatening", val = "Y", NULL)
    )
  ) %>%
  analyze(
    "AESHOSP",
    afun = a_freq_j,
    show_labels = "hidden",
    extra_args = append(
      extra_args_rr,
      list(label = "Requires or prolongs hospitalization", val = "Y", NULL)
    )
  ) %>%
  analyze(
    "AESDISAB",
    afun = a_freq_j,
    show_labels = "hidden",
    extra_args = append(
      extra_args_rr,
      list(
        label = "Persistent or significant disability/incapacity",
        val = "Y",
        NULL
      )
    )
  ) %>%
  analyze(
    "AESCONG",
    afun = a_freq_j,
    show_labels = "hidden",
    extra_args = append(
      extra_args_rr,
      list(label = "Congenital anomaly or birth defect", val = "Y", NULL)
    )
  ) %>%
  analyze(
    "AESMIE",
    afun = a_freq_j,
    show_labels = "hidden",
    extra_args = append(
      extra_args_rr,
      list(label = "Other medically important event", val = "Y", NULL)
    )
  )

result <- build_table(lyt, adae, alt_counts_df = adsl, round_type = "sas")

################################################################################
# Post-Processing:
# - Remove N's from Risk cols
# - Prune any categories with all zeros.
################################################################################

result <- remove_col_count(result)

# Blank out risk difference columns for duration/exposure rows (per spec, these should be empty)
if (risk_diff == TRUE) {
  cpaths <- col_paths(result)
  rr_col_idx <- which(sapply(cpaths, function(p) "rrisk_header" %in% p))
  if (length(rr_col_idx) > 0) {
    for (ci in rr_col_idx) {
      result[1, ci] <- rcell(NULL)
      result[2, ci] <- rcell(NULL)
    }
  }
}

result <- suppressWarnings(safe_prune_table(
  result,
  prune_func = count_pruner(
    cat_exclude = c(
      durfup_label,
      "Average exposure (number of administrations)",
      "AEs",
      "SAEs",
      rel_ae_label,
      rel_sae_label,
      "AE leading to permanent discontinuation of study treatment",
      rel_aedth_label,
      "Infections",
      "Serious infections",
      "Serious related infections",
      "Mild",
      "Moderate",
      "Severe",
      "Congenital anomaly or birth defect",
      acn_label
    )
  )
))

################################################################################
# Add titles and footnotes:
################################################################################

result <- set_titles(result, tab_titles)

################################################################################
# Convert to tbl file and output table:
################################################################################

tt_to_tlgrtf(string_map = string_map, tt = result, file = fileid, orientation = "landscape")
