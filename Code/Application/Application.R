### Load Packages
library(tidyverse)

### Load data
load('D:/Data/jtpa_national_evaluation/Formatted_Data/jtpa.RData')
raw <- readRDS('D:/Data/jtpa_national_evaluation/Formatted_Data/full_merge.RDS') |> 
  mutate(recid = as.integer(recid))

abadie_et_al <- x[1:9]
names(abadie_et_al) <- c('recid', '30_month_earnings', 'eligible', 'treated', 'male', 
              'high-school', 'black', 'hispanic', 'married') 
abadie_et_al <- abadie_et_al |> 
  mutate(recid = as.integer(recid))

data <- abadie_et_al |> inner_join(y = raw, by = 'recid') |> 
  arrange(by = 'recid')

### Remove Duplicated Columns
data_test <- data[!duplicated(as.list(data))] |> 
  # Remove Observations with non-matching ages
  filter(xor((age.x == age.y), (age.y == ''))) |>  
  # Filter out Observations with incorrect education codes
  filter(educ07 != 9) |> 
  filter(educ08 != 9) |> 
  filter(educ09 != 9) |> 
  filter(educ10 != 9) |> 
  filter(educ11 != 9) |> 
  filter(educ12 != 9) |> 
  filter(educ13 != 9) |> 
  filter(educ14 != 9) |> 
  filter(educ15 != 9) |> 
  filter(educ16 != 9) |> 
  filter(educ17 != 9) |> 
  filter(educ18 != 9) |> 
  # Remove Observations with incorrect number of education dummies
  mutate(across(starts_with('educ'), ~ as.integer(.x))) |> 
  filter(educ07 + educ08 + educ09 + educ10 + educ11 + educ12  + educ13 + educ14 + educ15 + educ16 + educ17 + educ18 == 1) |> 
  # filter out incorrect binaries
  filter(married %in% c(0,1), `high-school` %in% c(0,1)) |> 
  # Reformat Education
  mutate(educ = 7*educ07 + 8*educ08 + 9*educ09 + 10*educ10 
         + 11*educ11 + 12*educ12 + 13*educ13 + 14*educ14 
         + 15*educ15 + 16*educ16 + 17*educ17 + 18*educ18) |> 
  select(recid, `30_month_earnings`, eligible, treated, male.x, 
         `high-school`, black.x, hispanic.x, married, age.x, educ)
