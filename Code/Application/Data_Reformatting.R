### Load Packages for Data Manipulation
library(tidyverse)
library(foreign)

### Set path to data sets in Stata format
data_path <- 'D:/Data/jtpa_national_evaluation/application_data/'
subdirectory_paths <- list.dirs(path = data_path)[-1]
main_directory_paths <- list.files(path = data_path, pattern = '*.dta')

### Load Data Files
main_directory_data <- purrr::map(.x = main_directory_paths, 
                                  .f = ~read.dta(paste0(data_path, .x)))

subdirectory_files <- purrr::map(.x = subdirectory_paths, 
                                 .f = ~paste0(.x, '/', list.files(path = .x, pattern = '*.dta')))

subdirectory_data <- purrr::map(subdirectory_files, 
                                function(x) purrr::map(x, function(y) read.dta(file = y)))
                                                       

### Merge files
mdir_data_fulljoin <- main_directory_data |> reduce(full_join, by = 'recid')
IRS_data_fulljoin <- subdirectory_data[[4]] |> reduce(full_join, by = 'recid')
monthly_data_fulljoin <- subdirectory_data[[3]] |> reduce(full_join, by = 'recid')

### Merge 
data_merged <- list(mdir_data_fulljoin, IRS_data_fulljoin, monthly_data_fulljoin) |> 
  reduce(full_join, by = 'recid')

### Save merged data
saveRDS(object = data_merged, file = 'D:/Data/jtpa_national_evaluation/Formatted_Data/full_merge.RDS')

