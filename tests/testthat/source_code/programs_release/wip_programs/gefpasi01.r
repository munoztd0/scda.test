###############################################################
###                RELEASE VERSION - v2.0.0                 ###
###############################################################

################################################################################
## Original Reporting Effort: Standards
## Program Name:              gefpasi01.r
## R Version:                 4.5.2
## junco Version:             0.1.7
## Short Description:         Proportion of Subjects Achieving Response Through time Point
## Author:                    Technology Solutions
## Date:                      2026-09-302024
## Input:                     adsl.sas7bdat, adparspi.sas7bdat
## Output:                    GEFPASI01.png
## Remarks:
##
## Modification History:
## Rev #:
## Modified By:
## Reporting effort:
## Date:                      2026-09-30
## Description:
################################################################################

################################################################################
# Prep environment:
################################################################################

library(envsetup)
source(read_path(cl, 'utils_jjcs_internal.r'))
library(tern)


################################################################################
# Define output ID:
################################################################################

tblid <- "GEFPASI01"

################################################################################
# Get titles and footnotes:
################################################################################

# title_footer <- get_titles_internal(tblid)
string_map <- make_jj_str_map()

################################################################################
# Process data:
################################################################################

# reading data
adsl <- haven::read_sas(read_path(a_in, "adsl.sas7bdat")) |>
  filter(FASFL == "Y") |>
  group_by(TRT01P) |>
  summarise(n = n_distinct(USUBJID)) |>
  ungroup() |>
  mutate(
    trt = case_when(
      TRT01P == "PLACEBO" ~ "Placebo",
      TRT01P == "JNJ-77242113 25 MG QD" ~ "25 MG QD",
      TRT01P == "JNJ-77242113 50 MG QD" ~ "50 MG QD",
      TRT01P == "JNJ-77242113 25 MG BID" ~ "25 MG BID",
      TRT01P == "JNJ-77242113 100 MG QD" ~ "100 MG QD",
      TRT01P == "JNJ-77242113 100 MG BID" ~ "100 MG BID"
    ),
    trt = factor(
      paste0(trt, " (n=", n, ")"),
      # levels = c("Placebo (n=43)", "25 MG QD (n=43)", "50 MG QD (n=43)", "25 MG BID (n=41)", "100 MG QD (n=43)", "100 MG BID (n=42)"))
      levels = c(
        "Placebo (n=43)",
        "25 MG BID (n=41)",
        "25 MG QD (n=43)",
        "100 MG QD (n=43)",
        "50 MG QD (n=43)",
        "100 MG BID (n=42)"
      )
    )
  ) |>
  select(-n)

adparspi <- haven::read_sas(read_path(a_in, "adparspi.sas7bdat")) |>
  filter(
    FASFL == "Y",
    PARAMCD == "PASI75P",
    !is.na(AVALC),
    AVISIT %in% c("Week 1", "Week 2", "Week 4", "Week 8", "Week 12", "Week 16")
  ) |>
  mutate(
    visit = factor(
      AVISIT,
      levels = c("Week 1", "Week 2", "Week 4", "Week 8", "Week 12", "Week 16")
    )
  )

adparspi <- adparspi |>
  inner_join(adsl, by = "TRT01P")

# calculate percentage of subjects achieving response at each visit
denom <- adparspi |>
  group_by(trt, visit) |>
  summarise(denom = n_distinct(USUBJID)) |>
  ungroup()

count <- adparspi |>
  filter(AVALC == "Y") |>
  group_by(trt, visit) |>
  summarise(count = n()) |>
  ungroup() |>
  # fill in missing values for time points which don't have data
  tidyr::complete(trt, nesting(visit), fill = list(count = 0))

final <- denom |>
  left_join(count, by = c("trt", "visit")) |>
  mutate(
    visitn = case_when(
      visit == "Week 1" ~ 1,
      visit == "Week 2" ~ 2,
      visit == "Week 4" ~ 4,
      visit == "Week 8" ~ 8,
      visit == "Week 12" ~ 12,
      visit == "Week 16" ~ 16
    ),
    resp = count / denom,
    percent = tidytlg::roundSAS(resp * 100, 1), # Calculate percentage
    cil = tidytlg::roundSAS(
      100 * (resp - qnorm(0.975) * sqrt(resp * (1 - resp) / denom)),
      1
    ), # Lower CI
    ciu = tidytlg::roundSAS(
      100 * (resp + qnorm(0.975) * sqrt(resp * (1 - resp) / denom)),
      1
    ) # Upper CI
  ) |>
  mutate(
    cil = ifelse(cil < 0, 0, cil),
    ciu = ifelse(ciu > 100, 100, ciu)
  )

################################################################################
# Generate plot:
################################################################################

# define parameters for plotting:
# assign colorblind friendly palette: black, orange, green, dark blue, reddish-orange, pink
cbbPalette <- c(
  "#000000",
  "#E69F00",
  "#009E73",
  "#0072B2",
  "#D55E00",
  "#CC79A7"
)

pd <- position_dodge(0.4)

x_breaks <- c(1, 2, 4, 8, 12, 16)
y_breaks <- seq(0, 100, by = 10)
y_limits <- c(0, 100)
x_label <- "Week"
y_label <- "Percentage of Subjects (95% CI)"


# line plot ----------------------------------------------

pt <- final |>
  ggplot(aes(
    x = visitn,
    y = percent,
    group = trt,
    color = trt,
    linetype = trt,
    shape = trt
  )) +
  geom_errorbar(aes(ymin = cil, ymax = ciu), width = 0.9, position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, fill = "white", size = 2) + # white fill inside shape

  # assign colorblind friendly palette
  scale_color_manual(values = cbbPalette) +
  scale_shape_manual(values = c(15, 16, 17, 21, 22, 24)) + # Change shapes

  # define x-axis breaks and limits
  scale_x_continuous(breaks = x_breaks) +
  # define y-axis breaks and limits
  scale_y_continuous(breaks = y_breaks, limits = y_limits) +
  labs(
    x = x_label,
    y = y_label
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 9, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    # axis.title.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  )


################################################################################
# Create png and output file:
################################################################################

# create png file and output figure
pname <- paste0(tolower(tblid), ".png")

png(
  write_path(opath, pname),
  width = 22,
  height = 14,
  units = "cm",
  res = 300,
  type = "cairo"
)
print(write_path(opath, pname)) ### print png path and name in log
print(pt)
dev.off()

title <- "Subjects Achieving PASI 75 Response Through Week 16; Full Analysis Set (Study 77242113PSO2001)"
main_footer <- ""

tidytlg::gentlg(
  tlf = "g",
  plotnames = write_path(opath, pname),
  plotwidth = 8,
  orientation = "landscape",
  opath = write_path(opath),
  file = tblid,
  title = title, # title_footer$title
  footers = main_footer
) # title_footer$main_footer
