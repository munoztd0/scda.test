###############################################################################
## Original Reporting Effort: Standards
## Program Name:              sas2rdsadexsum.r
## R version:                 4.2.1
## Short Description:         Create adexsum.rds
## Author:                    Technology Solutions
## Date:                      2024-10-24
## Input:                     ADEXSUM, ADSL, SDTM Metadata, ADaM Metadata
## Output:                    adexsum.rds
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

source(read_path(cl, "utils_jjcs_internal.r"))
library(envsetup)
library(dplyr)


domain <- "adexsum"


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

### Part 0 :
#### assign alternative labels -- to be used on displays
### to draft a proposed excel file + function to perform this step

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

# manually add some extra --- harder to automate - leave a manual process for now

# vars <- c(vars,"CQ01NAM")
# sortvars <- c(sortvars,"CQ01NAM")

# Need to add all CQxxNAM vars into this I think, these will also
# need to be factors.
# Not sure why we aren't converting all char vars to factors
# Need to set missing to NA_inter_ or NA_character in this process?

# this part is not working so commenting out for now
df <- convert2factor(df, vars = vars, sortvars = sortvars)

### replace blank with NA for remaining character variables - and convert to factor?
df <- df_na(df, char_as_factor = TRUE)


# Part 3: add variables from adsl
df <- df %>%
  left_join(., adsl, by = "USUBJID")

# Part 4:  restore any dropped labels from the original dataset
df <- restore_labels(df, df_orig)

# Part 5: have variables in same order as in the original dataset, then new vars, then the _decode vars
df <- set_varorder(df, df_orig)

# Part 6:  store  the dataset in rds format in the analysis (a_out) folder
saveRDS(df, file = envsetup::write_path(a_out, paste0(domain, ".rds")))
