Plasma.Working.Corrected <- NULL
Total.Samples <- ncol(plasma_corrected_data)/2
for(j in 1:Total.Samples)
{
  lwr <- 2*j - 1
  upr <- 2*j
  if(j %% 20 == 0) cat(j, ' of ', Total.Samples, 'completed. \n')
  temp.col <- plasma_corrected_data %>% select(lwr:upr)
  temp.col <- temp.col %>% mutate(SampleID = colnames(temp.col)[2])
  colnames(temp.col)[1:2] <- c('Temperature', 'dCp')
  Plasma.Working.Corrected <- Plasma.Working.Corrected %>% rbind(temp.col)
}

Plasma.Working.Corrected <- Plasma.Working.Corrected %>% 
  filter(!is.na(Temperature)) %>%
  mutate(SampleNumber = factor(str_extract(SampleID, '\\d+')),
         SampleIteration = factor(str_extract(SampleID, '[a-f]$'))) %>%
  mutate(SampleID = factor(SampleID)) %>%
  relocate(SampleID)


thing2 <- Plasma.Working.Corrected %>% filter(SampleID == All.IDs[10]) %>% filter(between(Temperature, 45, 90))
g2 <- ggplot(thing2, aes(x = Temperature, y = dCp)) + geom_line()
