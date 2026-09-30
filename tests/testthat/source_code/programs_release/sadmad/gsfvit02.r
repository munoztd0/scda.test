###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gsfvit02.r
## R version:                 4.5.2
## junco Version:             0.1.6
## Short Description:         Individual Profile Over Time  – [SAD/MAD] [Part 1]
## Author:                    C&SP Methodology
## Date:                      2026-09-30
## Input:                     advs
## Output:                    gsfvit02.rtf
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
################################################################################

################################################################################
# Prep environment:
################################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(junco)
library(dplyr)
library(stringr)
library(ggplot2)

################################################################################
# Define output ID:
################################################################################

tblid <- "GSFVIT02"

popfl <- "SAFFL"

trtvar <- "TRT01A"

studypart_var <- "PARTC"
studyprt_val <- "PART 1"

actarm_str <- "SAD"

paramcd <- "SYSBPU"

# Split Timepoint into groups if not fit into single page
# default = NULL (no split)
# timepoint_split <- NULL

## user can comment this part if no split needed
# else User can define timepoints here
timepoint_split <- list(
  time_point_1 = c(
    "Period 1 Day 1, 1hr",
    "Period 1 Day 1, 2hr",
    "Period 1 Day 1, 3hr",
    "Period 1 Day 1, 4hr",
    "Period 1 Day 1, 5hr",
    "Period 1 Day 1, 6hr",
    "Period 1 Day 1, 8hr",
    "Period 1 Day 1, 12hr"
  ),
  time_point_2 = c(
    "Period 1 Day 2",
    "Period 1 Day 3 Unscheduled 1",
    "Period 1 Day 6",
    "FU/EOS/ET"
  )
)

################################################################################
# Get titles and footnotes:
################################################################################

title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################


advs <- haven::read_sas(envsetup::read_path(a_in, "advs.sas7bdat")) %>%
  df_na()

advs0 <- advs %>%
  dplyr::filter(
    !!rlang::sym(studypart_var) == studyprt_val,
    grepl(actarm_str, ACTARM, ignore.case = TRUE),
    !!rlang::sym(popfl) == "Y",
    !is.na(AVAL),
    PARAMCD %in% paramcd
  ) %>%
  dplyr::mutate(
    !!trtvar := if_else(
      .data[[trtvar]] == "Placebo",
      "Pooled Placebo",
      .data[[trtvar]]
    ),
    !!trtvar := factor(.data[[trtvar]])
  ) %>%
  dplyr::select(
    USUBJID,
    SUBJID,
    PARAM,
    all_of(trtvar),
    AVAL,
    AVISIT,
    AVISITN,
    ATPT,
    ATPTN
  )

# derive selvisit from data, ordered by AVISITN then ATPTN
selvisit <- advs %>%
  mutate(
    # if not required this can be removed, update it as per study requirement
    # removing "Period 1" and "Period 2" from AVISIT
    VISIT = str_trim(str_replace_all(AVISIT, "Period 1 |Period 2 ", "")),
    # replacing "FU/EOS/ET" string with "End of Study" string in VISIT
    VISIT = if_else(
      grepl("FU/EOS/ET", VISIT),
      sub("FU/EOS/ET", "End of Study", VISIT),
      VISIT
    )
  ) %>%
  arrange(AVISITN, ATPTN) %>%
  distinct(AVISIT, VISIT) %>%
  # VISITN is used on the x-axis to reflect the relative distance between visits
  mutate(VISITN = dplyr::row_number()) %>%
  select(AVISIT, VISIT, VISITN)

# join ordered visits with advs0
advs1 <- advs0 %>%
  inner_join(selvisit, by = c("AVISIT")) %>%
  group_by(USUBJID, !!rlang::sym(trtvar)) %>%
  arrange(USUBJID, VISITN) %>%
  ungroup()

################################################################################
# Generate plot:
################################################################################

plot_fun <- function(df, trt) {
  df %>%
    filter(!!rlang::sym(trtvar) %in% trt) %>%
    ggplot(
      aes(
        x = VISITN,
        y = AVAL,
        group = SUBJID,
        color = factor(SUBJID),
        linetype = factor(SUBJID),
        shape = factor(SUBJID)
      )
    ) +
    geom_line() +
    geom_point() +
    scale_shape_manual(values = 1:length(unique(df$SUBJID))) +
    scale_x_continuous(
      breaks = unique(df$VISITN),
      labels = unique(df$VISIT)
    ) +
    labs(
      title = unique(df$PARAM),
      subtitle = trt,
      x = "Timepoints",
      y = "Value",
      color = NULL,
      linetype = NULL,
      shape = NULL
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 9, color = "black", family = "Arial"),
      plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 9, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 9),
      axis.text.x = element_text(
        angle = -45,
        hjust = 0,
        vjust = 1
      ),
      plot.margin = margin(r = 20), # keep space on right for x axis labels
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 9, face = "bold")
    )
}

################################################################################
# Generate PNGs and RTF per category:
################################################################################

split_active <- exists("timepoint_split") &&
  !is.null(timepoint_split) &&
  length(timepoint_split) > 0

split_indices <- if (split_active) seq_along(timepoint_split) else 1

# List to store the generated plots
plot_list <- list()

for (trt in unique(advs1[[trtvar]])) {
  # Loop through each treatment
  for (j in split_indices) {
    # Loop through each timepoint

    if (split_active) {
      advs_mod <- advs1 %>%
        filter(AVISIT %in% timepoint_split[[j]])

      split_name <- names(timepoint_split)[j]
    } else {
      # Use all timepoints
      advs_mod <- advs1
      split_name <- "all_timepoints"
    }

    if (nrow(advs_mod) == 0) {
      warning(
        paste(
          "No data available for treatment:",
          trt,
          "and timepoint split:",
          split_name
        )
      )
      next
    }

    plot <- plot_fun(df = advs_mod, trt = trt)

    # Create treatment name for storing the plot in the list
    trt_name <- gsub("[^A-Za-z0-9]+", "_", trt)

    plot_name <- paste0(
      tolower(tblid),
      "_",
      trt_name,
      "_",
      split_name
    )

    # Store the plot
    plot_list[[plot_name]] <- plot
  }
}


plot_files <- purrr::imap_chr(
  plot_list,
  ~ {
    pname <- paste0(.y, ".png")

    png(
      filename = write_path(opath, pname),
      width = 22,
      height = 14,
      units = "cm",
      res = 300,
      type = "cairo"
    )

    print(write_path(opath, pname))
    print(.x)

    dev.off()

    write_path(opath, pname)
  }
)

tidytlg::gentlg(
  tlf = "g",
  plotnames = plot_files,
  plotwidth = 8,
  orientation = "landscape",
  opath = write_path(opath),
  file = tblid,
  title = title_footer$title,
  footers = title_footer$main_footerL
)
