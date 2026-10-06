library(envsetup)
source(read_path(cl, "utils_jjcs_internal.r"))
source(read_path(cl, "qc_utils_jj.r"))
library(dplyr)
library(openxlsx)
library(compareDF)


# Program Background ----
# purpose of this program: creation of codelistinfo.xlsx file in dpspath folder for utilization in sas2rdsadxx.r programs

#### INTENDED USAGE of file codelistinfo.xlsx is that it will be shared in submission with FDA
#### and can be used by FDA to recreate rds datasets from ADaM sas source

# this program (create_codelistinfo_xlsx.r) should not be shared with FDA, as it needs APT/internal metadata as input

# the switch from APT/internal ADaM metadata files to codelistinfo.xlsx file is handled in your .Rprofile with
# method_sas2rds object
# note that the following choices are available :
# INTERNAL (using APT or internal JJ ADaM spreadsheets) or
# XLSX (using codelistinfo.xlsx)

# do not switch too early, as frequent reruns of this create_codelistinfo_xlsx program would be needed
# once your domains are in stable stage, you can consider switching the method_sas2rds in your .Rprofile

# due to write permission restrictions in SPACE folder dpspath this program can be executed in PDEV only, not in PREPROD/PROD)
# the codelistinfo.xlsx that this program generates will have to be promoted to PREPROD/PROD within the SPACE system
# the a_in_fixed_search setting (see below), is the way to control which envsetup environment (PDEV/PREPROD/PROD)
# for reading ADaM datasets from a_in for the creation of this file (for evaluation of observed values for char variables)
# it is recommended to use PREPROD and at time of ADaMs (sas version) available in PROD, switch to PROD

# Main Program Flow ----
# 1/ creation of codelistinfo.xlsx in dpspath
# 2/ confirmation that codelistinfo.xlsx can be used in sas2rds programs as alternative to calling APT
# Check that file "getcodelist_issues.xlsx" is NOT created in dpspath and that log contains the message
#
# Note: All required information from codelistinfo could be regenerated
# from excel file codelistinfo.xlsx for further processing.

# Actual Program start ----

a_in_fixed_search <- "PREPROD"
# STEP 1/ creation of codelistinfo.xlsx in dpspath ----

# retrieving codelist information from :
# 1/ input ADaM dataset
# 2/ SDTM metadata through metadata from source
# 3/ ADaM metadata through APT/or xlsx files
### this is mimicking adattrib and esub process

getcodelist_singledomain <- function(
  domain,
  envsetup_environ = Sys.getenv("ENVSETUP_ENVIRON")
) {
  domain <- tolower(domain)
  domain <- sub(".sas7bdat", "", domain, fixed = TRUE)

  # ignore non adam datasets
  if (!stringr::str_starts(domain, "ad")) {
    return(NULL)
  }

  df <- haven::read_sas(
    envsetup::read_path(a_in, paste0(domain, ".sas7bdat"), envsetup_environ = envsetup_environ)
  )

  code_decode <- getcodelistinfo(
    df = df,
    domain = domain,
    APT = apt,
    adam_meta_loc = get_path(am_in),
    sdtm_meta_loc = get_path(dm_in)
  )

  code_decode$DOMAIN <- domain
  code_decode <- code_decode |>
    relocate("DOMAIN")

  code_decode
}

### get all available ADaM domains from a_in folder (PREPROD/PROD)
domains <- list.files(get_path(a_in)[[a_in_fixed_search]], pattern = "^ad.*\\.sas7bdat$")

code_decode_all <- lapply(domains, getcodelist_singledomain, envsetup_environ = a_in_fixed_search)

# combine into single dataframe
code_decode_df <- bind_rows(code_decode_all)

code_decode_df <- code_decode_df |>
  mutate(across(where(is.character), ~ na_if(.x, "")))

# investigate whether/why there are records with CODEVAL missing
code_decode_df_check <- code_decode_df |>
  filter(is.na(CODEVAL))

# exclude these for further processing
code_decode_df <- code_decode_df |>
  filter(!is.na(CODEVAL))

# save information in xlsx file in dpspath
create_codelistinfo_xlsx(code_decode_df, file = write_path(dpspath, "codelistinfo.xlsx"))


# STEP 2/ check that retrieval process when using codelistinfo.xlsx should be safe  ----

message("------------------------------------------------")
message("STEP2: Check for codelistinfo retrieval process ")
message("------------------------------------------------")

### Here perform the check that the code_decode_df dataframe can be reconstructed from the codelistinfo.xlsx file
# confirmation that with these datasets the required information from combined can be recreated

# keep only relevant information
#### !!!!!!!!!!! PARAMCD is needed in the processing of convert2factor
#### need to keep !!!
target_code_decode <- unique(code_decode_df |> select(-c("defined", "observed"))) |>
  arrange("DOMAIN", "VARNAME", "PARAMCD", "CODELST", "RNK")

reconstructed_code_decode <- get_codelistinfo_xlsx(
  domain = NULL,
  file = read_path(dpspath, "codelistinfo.xlsx")
) |>
  select(-c("defined", "observed"))

# compare_df requires that the dataset contain the same columns!!!
stopifnot(identical(
  sort(names(reconstructed_code_decode)),
  sort(names(target_code_decode))
))

comparison_output <- compare_df_ext(
  df_new = reconstructed_code_decode,
  df_old = target_code_decode,
  group_col = c("DOMAIN", "VARNAME", "PARAMCD", "CODELST", "RNK", "DECOD"),
  exclude = NULL,
  tolerance = 0,
  tolerance_type = "ratio",
  keep_unchanged_rows = FALSE,
  keep_unchanged_cols = TRUE,
  change_markers = c("+", "-", "="),
  round_output_to = 3
)

if (comparison_output$compareDF_error == TRUE) {
  message("!!!!!!!!!!ERROR: Please investigate why compareDF resulted in error")
} else if (sum(comparison_output$change_summary[c("changes", "additions", "removals")]) > 0) {
  create_output_table(
    comparison_output,
    output_type = "xlsx",
    file_name = write_path(dpspath, "getcodelist_issues.xlsx"),
    limit = 50,
    color_scheme = c(
      addition = "#52854C",
      removal = "#FC4E07",
      unchanged_cell = "#999999",
      unchanged_row = "#293352"
    ),
    headers = NULL,
    change_col_name = "chng_type",
    group_col_name = "grp"
  )

  message(
    "!!!!!!!!!!ERROR: All required information from codelistinfo could NOT be regenerated
          from excel file codelistinfo.xlsx for further processing."
  )
} else {
  message(
    "Note: All required information from codelistinfo could be regenerated
          from excel file codelistinfo.xlsx for further processing."
  )
}
