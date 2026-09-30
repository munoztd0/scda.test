###############################################################################
## Original Reporting Effort: Standards
## Program Name:              sas2rdsadagocmq.r
## R version:                 4.2.1
## Short Description:         Create adagocmq.rds
## Author:                    Technology Solutions
## Date:                      2025-09-30
## Input:                     ADAGOCMQ, ADSL, SDTM Metadata, ADaM Metadata
## Output:                    adagocmq.rds
## Remarks:                   Contains OCMQ data with fields:
##                            USUBJID, ATERMN, ATERM
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

source(read_path(cl, "utils_jjcs_internal.r"))
library(envsetup)
library(dplyr)

domain <- "adagocmq"

## Read in dataset
df <- haven::read_sas(envsetup::read_path(a_in, paste0(domain, ".sas7bdat")))
df_orig <- df

adsl <- readRDS(envsetup::read_path(a_in, "adsl.rds"))

vars <- setdiff(names(df), names(adsl))
vars <- c("USUBJID", vars)
adslvars <- intersect(names(adsl), names(df))
adslvarsdecodelst <- paste0(adslvars, "_DECODE")
adslvarsdecode <- intersect(names(adsl), adslvarsdecodelst)
adslvarsall <- c(adslvars, adslvarsdecode)
adsl <- adsl[adslvarsall]

df <- df %>%
  select(all_of(vars))

## retrieve a code_decode dataframe using get_codelistinfo_method and df
# note : method_sas2rds should be defined in your .Rprofile
code_decode <- get_codelistinfo_method(
  df = df,
  domain = domain,
  method = method_sas2rds,
  file = read_path(dpspath, "codelistinfo.xlsx"),
  APT = apt,
  adam_meta_loc = get_path(am_in),
  sdtm_meta_loc = get_path(dm_in),
  define_loc = get_path(sadpath)
)

### Actual conversion
## step 1 -- convert vars into factors and add decode versions (if needed)
df <- convert2factor(df, code_decode = code_decode)

## step 2 -- get vars that have a numeric counterpart and were not handled in step1
vars_init <- names(df)
vars_done <- unique(code_decode$VARNAME)

vars_init2 <- setdiff(vars_init, vars_done)
sortvars <- paste0(vars_init2, "N")
selvars <- sortvars %in% vars_init2
sortvars <- sortvars[selvars]
vars <- vars_init2[selvars]

# Ensure ATERM is properly handled
df <- mutate(df, ATERM = case_when(ATERM == "" ~ "Uncoded", .default = ATERM))

# Convert remaining variables to factors
df <- convert2factor(df, vars = vars, sortvars = sortvars)

### replace blank with NA for remaining character variables - and convert to factor
df <- df_na(df, char_as_factor = TRUE)

# Add variables from adsl
df <- df %>%
  left_join(., adsl, by = "USUBJID")

# Restore any dropped labels from the original dataset
df <- restore_labels(df, df_orig)

# Have variables in same order as in the original dataset, then new vars, then the _decode vars
df <- set_varorder(df, df_orig)

# Store the dataset in rds format in the analysis (a_out) folder
saveRDS(df, file = envsetup::write_path(a_out, paste0(domain, ".rds")))
