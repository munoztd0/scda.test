###############################################################################
## Original Reporting Effort: Standards
## Program Name:              sas2rdsadsl.r
## R version:                 4.2.1
## Short Description:         Create adsl.rds
## Author:                    Technology Solutions
## Date:                      2024-10-24
## Input:                     ADSL, SDTM Metadata, ADaM Metadata
## Output:                    adsl.rds
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


domain <- "adsl"


## Read in dataset
df <- haven::read_sas(envsetup::read_path(a_in, paste0(domain, ".sas7bdat")))
df_orig <- df

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
### drafted a study/domain specific proposal for alternative labels  + function to perform this step

df <- set_alternative_labels(df, domain)


### updates to code_decode to re-order / drop levels for some variables
### this is not yet a well-defined process and requires manual judgement on which to re-order, which vars to contain only observed levels
### some of these would to be used across all sas2rds programs
code_decode_adjusted <- code_decode %>%
  mutate(
    RNK = case_when(
      VARNAME == "SEX" & CODEVAL == "M" ~ "0001",
      VARNAME == "SEX" & CODEVAL == "F" ~ "0002",
      TRUE ~ RNK
    )
  )


code_decode_adjusted <- code_decode_adjusted %>%
  ### might be study specific
  ### get rid of unobserved levels for some variables (not all)
  ### there is no clear list of which variables to include in the following
  filter(!(VARNAME %in% c("RACE", "COUNTRY", "REGION1") & is.na(observed))) %>%
  filter(
    !(VARNAME == "SEX" &
      is.na(observed) &
      CODEVAL %in% c("UNDIFFERENTIATED", "U"))
  )


### Actual conversion
## step 1 -- convert vars into factors and add decode versions if needed
df <- convert2factor(df, code_decode = code_decode_adjusted)

## step 2 -- get vars that have a numeric counterpart and were not handled in step1
vars_init <- names(df)
vars_done <- unique(code_decode$VARNAME)

vars_init2 <- setdiff(vars_init, vars_done)
sortvars <- paste0(vars_init2, "N")
indices <- sortvars %in% vars_init2
sortvars <- sortvars[indices]
vars <- vars_init2[indices]

df <- convert2factor(df, vars = vars, sortvars = sortvars)

code_decode2 <- tribble(
  ~VARNAME                                    ,
  ~CODELST                                    ,
  ~CODEVAL                                    ,
  ~defined                                    ,
  ~PARAMCD                                    ,
  ~source                                     ,
  ~observed                                   ,
  ~DECOD                                      ,
  ~RNK                                        ,
  "ETHNIC"                                    ,
  "ETHNIC"                                    ,
  "HISPANIC OR LATINO"                        ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Hispanic or Latino"                        ,
  "00001"                                     ,
  "ETHNIC"                                    ,
  "ETHNIC"                                    ,
  "NOT HISPANIC OR LATINO"                    ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Not Hispanic or Latino"                    ,
  "00002"                                     ,
  "ETHNIC"                                    ,
  "ETHNIC"                                    ,
  "NOT REPORTED"                              ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Not reported"                              ,
  "00003"                                     ,
  "ETHNIC"                                    ,
  "ETHNIC"                                    ,
  "UNKNOWN"                                   ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Unknown"                                   ,
  "00004"                                     ,
  "SEX"                                       ,
  "SEX"                                       ,
  "M"                                         ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Male"                                      ,
  "00001"                                     ,
  "SEX"                                       ,
  "SEX"                                       ,
  "F"                                         ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Female"                                    ,
  "00002"                                     ,
  "SEX"                                       ,
  "SEX"                                       ,
  "INTERSEX"                                  ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Intersex"                                  ,
  "00003"                                     ,
  "SEX"                                       ,
  "SEX"                                       ,
  "U"                                         ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Unknown"                                   ,
  "00004"                                     ,
  "RACE"                                      ,
  "RACE"                                      ,
  "AMERICAN INDIAN OR ALASKA NATIVE"          ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "American Indian or Alaska Native"          ,
  "00001"                                     ,
  "RACE"                                      ,
  "RACE"                                      ,
  "ASIAN"                                     ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Asian"                                     ,
  "00002"                                     ,
  "RACE"                                      ,
  "RACE"                                      ,
  "BLACK OR AFRICAN AMERICAN"                 ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Black or African American"                 ,
  "00003"                                     ,
  "RACE"                                      ,
  "RACE"                                      ,
  "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Native Hawaiian or other Pacific Islander" ,
  "00004"                                     ,
  "RACE"                                      ,
  "RACE"                                      ,
  "WHITE"                                     ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "White"                                     ,
  "00005"                                     ,
  "RACE"                                      ,
  "RACE"                                      ,
  "MULTIPLE"                                  ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Multiple"                                  ,
  "00006"                                     ,
  "RACE"                                      ,
  "RACE"                                      ,
  "NOT REPORTED"                              ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Not reported"                              ,
  "00007"                                     ,
  "RACE"                                      ,
  "RACE"                                      ,
  "UNKNOWN"                                   ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Unknown"                                   ,
  "00008"                                     ,
  "RACE"                                      ,
  "RACE"                                      ,
  "OTHER"                                     ,
  "Y"                                         ,
  NA                                          ,
  "Table Shell"                               ,
  "Y"                                         ,
  "Other"                                     ,
  "00009"                                     ,
)

df <- convert2factor(df, code_decode = code_decode2)

### replace blank with NA for remaining character variables
df <- df_na(df)

# Part 4:  restore any dropped labels from the original dataset
df <- restore_labels(df, df_orig)

# Part 5: have variables in same order as in the original dataset, then new vars, then the _decode vars
df <- set_varorder(df, df_orig)

# Part 6:  store  the dataset in rds format in the analysis (a_out) folder
saveRDS(df, file = envsetup::write_path(a_out, paste0(domain, ".rds")))
