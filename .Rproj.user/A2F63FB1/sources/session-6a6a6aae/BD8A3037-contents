library(dplyr)
library(tidyverse)
library(stringr)

Urine_Raw_Corrected <- readxl::read_excel("data_raw/Urine samples thermograms raw 4.23.2023 corrected.xlsx", 
                                            sheet = "Data")

Check.na.cols <- Urine_Raw_Corrected %>% slice(1) %>% is.na() %>% which()
colnames(Urine_Raw_Corrected)[Check.na.cols]
Urine_Raw_Corrected_Full <- Urine_Raw_Corrected %>% select(-all_of(Check.na.cols))

Urine.Working <- NULL
Total.Samples <- ncol(Urine_Raw_Corrected_Full)/2
for(j in 1:Total.Samples)
{
  lwr <- 2*j - 1
  upr <- 2*j
  if(j %% 20 == 0) cat(j, ' of ', Total.Samples, 'completed. \n')
  temp.col <- Urine_Raw_Corrected_Full %>% select(lwr:upr)
  temp.col <- temp.col %>% mutate(SampleID = colnames(temp.col)[2])
  colnames(temp.col)[1:2] <- c('Temperature', 'dCp')
  Urine.Working <- Urine.Working %>% rbind(temp.col)
}
### Create variables for patient identification and tracking
### Remove all NA rows (non-aligned temperatures)
Urine.Working <- Urine.Working %>% 
  filter(!is.na(Temperature)) %>%
  mutate(SampleNumber = factor(str_extract(SampleID, '\\d+')),
         SampleIteration = factor(str_extract(SampleID, '[a-f]'))) %>%
  mutate(SampleID = factor(SampleID)) %>%
  relocate(SampleID)

rm(list = setdiff(ls(), c('Urine.Working')))
save.image('data/Urine_Working.RData')
