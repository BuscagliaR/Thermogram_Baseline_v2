library(tidyverse)
library(knitr)
library(forecast)
load('data/Urine_Working.RData')

final.sample.interpolate <- function(x, grid.temp, plot.on = TRUE)
{
  spline.fit <- smooth.spline(x$Temperature, x$dCp, cv = TRUE)
  interpolated.sample.pred <- predict(spline.fit, grid.temp)
  interpolated.sample <- data.frame(Temperature = grid.temp, 
                                    dCp = interpolated.sample.pred$y)
  
  if(plot.on)
  {
    g.out <- ggplot(interpolated.sample, aes(x = Temperature, y = dCp)) + geom_point() +
      labs(title = 'Interpolated Result')
    print(g.out)
  }
  
  return(interpolated.sample)
}

grid.temp <- seq(45, 90, 0.1)
final.sample <- data.frame(Temperature = grid.temp)
for(i in 1:length(All_Urine_ID)){
  sample <- Urine.Working.Final %>% filter(SampleID ==All_Urine_ID[i])
  single.sample <- final.sample.interpolate(
    x = sample, 
    grid.temp = seq(45, 90, 0.1), plot.on = FALSE)
  final.sample <- final.sample %>% cbind(out = single.sample$dCp)
  print(i)
}

All.IDs <- Urine.Working.Final %>% pull(SampleID) %>% unique() %>% as.vector()
colnames(final.sample)[-1] <- All.IDs

rm(list = setdiff(ls(), c('final.sample')))
write.csv(x = final.sample, file = 'generated_output/Interpolated_truncated_data.csv')
