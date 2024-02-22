library(tidyverse)
library(knitr)
library(forecast)
library(ggplot2)
load('data/Urine_Working.RData')

head(table.innermost.smaller)

table.innermost.smaller <- table.innermost.smaller %>% mutate(SampleNumber = factor(str_extract(SampleID, '\\d+')))

Endpoint_Table <- read.csv("~/GitHub/Thermogram_Baseline_v2/generated_output/Endpoint_Table.csv")

Endpoint_Table <- Endpoint_Table %>% mutate(SampleNumber = factor(str_extract(SampleID, '\\d+')))

Endpoint_Table <- Endpoint_Table %>% mutate(New_lower = filter)           

Inner_Endpoint_Table <- Endpoint_Table %>% filter(method == "innermost")

max <- 393

new_frame <- data.frame()
for(i in 1:max){
  table <- table.innermost.smaller %>% filter(i == SampleNumber)
  new_lower <- max(table$lower)
  table <- table %>% mutate(new_lower = new_lower)
  new_upper <- min(table$upper)
  table <- table %>% mutate(new_upper = new_upper)
  new_frame <- new_frame %>% rbind(table)
}

grid.temp <- seq(54, 86, 0.1)

### Setup frame to store all data generated
Urine.Final.Data <- data.frame(Temperature = grid.temp)

Urine.Final.Data <- data.frame(Temperature = grid.temp)
for(i in 1:length(All_Urine_ID)){
  working.sample <- Urine.Working.Final %>% filter(SampleID == All_Urine_ID[i])
  endpoints <- new_frame %>% filter(SampleID == All_Urine_ID[i]) %>% select(new_lower,new_upper)
  baseline.output <- baseline.subtraction.byhand(
    x = working.sample, 
    lwr.temp = endpoints$new_lower, 
    upr.temp = endpoints$new_upper,
    plot.on = FALSE)
  final.sample <- final.sample.interpolate(
    x = baseline.output, 
    grid.temp = seq(54, 86, 0.1),
    plot.on = FALSE)
  Urine.Final.Data <- Urine.Final.Data %>% cbind(out = final.sample$dCp)
}

colnames(Urine.Final.Data)[-1] <- All_Urine_ID


write.csv(x = Urine.Final.Data, file = 'generated_output/Same_Inner_Endpoint_Graphs.csv')

#rm(list=setdiff(ls(),c("All_Urine_ID", "Inner_Endpoint_Table", "Urine.Final.Data")))

#save.image('generated_output/Inner_Average.RData')

pdf('generated_output/Same_Endpoint_Graphs.pdf')
{
  gs <- NULL
  for(i in 1:n.samples)
  {
    g1 <- Urine.Final.Data %>% select(Temperature, All_Urine_ID[i])
    gs <- ggplot(g1, aes(x = Temperature,y = g1[,2]))+
      geom_point(size = 0.5)+
      labs(title = paste0('Final Automated Sample for ', str_sub(All_Urine_ID[i], 1)))
    print(gs)
  }
}
dev.off()

Urine.Final.Data <- Urine.Final.Data %>% mutate(SampleNumber = factor(str_extract(SampleID, '\\d+')))

t <- colnames(Urine.Final.Data)[-1]
numbers <- str_extract(t, '\\d+')
numbers <- as.integer(numbers)

new_frame <- Urine.Final.Data

colnames(new_frame)[-1] <- numbers


urine2 <- Urine.Final.Data %>% pivot_longer(names_to = "Sample_ID", values_to = "dCp", 2:782)

urine2 <- urine2 %>% mutate(sample_number = str_extract(Sample_ID, '\\d+'), sample_rep = str_extract(Sample_ID, '\\D+'))


numbers <- urine2 %>% pull(sample_number) %>% unique()

pdf('generated_output/Same_Inner_Endpoint_Graphs_Double.pdf')
{
  gs <- NULL
  for(i in 1:length(numbers))
  {
    thing <- urine2 %>% filter(sample_number == numbers[i])
    gs <- ggplot(thing, aes(x = Temperature,y = dCp))+
      geom_line(linewidth = 1.5, aes(color = sample_rep))+
      labs(title = paste0('Final Automated Sample for sample number ', numbers[i]))
    print(gs)
  }
}
dev.off()


###Functions
run.indicators <- function(x) return(ifelse(x <= 0, -1, 1))

updown.count <- function(x){
  require(tidyverse)
  j = 1
  n.its <- length(x)
  df.out <- data.frame()
  count <- 1
  while(j <= n.its){
    if(x[j] == x[j+1] & j != n.its)
    {
      j <- j+1
      count <- count+1
    } else {
      out <- c(count, run.indicators(x[j]))
      df.out <- df.out %>% rbind(out)
      j <- j+1
      count <- 1
    }
  }
  colnames(df.out) <- c('RunSize', 'Up/Down')
  df.out$`Up/Down` <- ifelse(df.out$`Up/Down` < 0, 'Down', 'Up')
  return(df.out)
}

baseline.subtraction.byhand <- function(x, lwr.temp, upr.temp, plot.on = TRUE)
{
  ### check-conditions of boundaries - this effects automation
  if(lwr.temp < min(x$Temperature)+1) lwr.temp = lwr.temp + 1
  if(upr.temp > max(x$Temperature)-1) upr.temp = upr.temp - 1
  ### Extract the baseline regions
  work.lower <- x %>% filter(Temperature < lwr.temp)
  work.upper <- x %>% filter(Temperature > upr.temp)
  ### Splines for lower/upper regions
  spline.lower <- smooth.spline(work.lower$Temperature, work.lower$dCp, cv = TRUE)
  spline.upper <- smooth.spline(work.upper$Temperature, work.upper$dCp, cv = TRUE)
  ### Store data for graphing
  spline.lower.fit <- data.frame(Temperature = work.lower$Temperature, fit = spline.lower$y)
  spline.upper.fit <- data.frame(Temperature = work.upper$Temperature, fit = spline.upper$y)
  
  ### store middle (signal) region
  work.mid <- x %>% filter(between(Temperature, lwr.temp, upr.temp))
  ### find endpoints of splines
  spline.connect.points <- rbind(
    spline.lower.fit %>% filter(Temperature == max(Temperature)),
    spline.upper.fit %>% filter(Temperature == min(Temperature)))
  ### connect endpoints and store
  spline.connect.lm <- lm(fit ~ Temperature, data = spline.connect.points)
  spline.connect.fit <- data.frame(
    Temperature = work.mid$Temperature, 
    fit = predict(spline.connect.lm, data.frame(Temperature = work.mid$Temperature)))
  
  ### store baseline as one unit
  working.baseline.final <- rbind(spline.lower.fit, spline.connect.fit, spline.upper.fit)
  ### join for tidyverse simplification
  baseline.join <- full_join(x, working.baseline.final, by = 'Temperature')
  ### final sample!
  baseline.sample <- baseline.join %>% mutate(final.dcp = dCp - fit) %>%
    select(Temperature, final.dcp) %>% rename(dCp = final.dcp)
  
  if(plot.on)
  {
    ### graph of raw with spline
    g.spline <- working.sample %>% ggplot(aes(x = Temperature, y = dCp)) + geom_point() + 
      geom_line(data = spline.lower.fit, aes(x = Temperature, y = fit), color = 'red') + 
      geom_line(data = spline.upper.fit, aes(x = Temperature, y = fit), color = 'red') +
      geom_line(data = spline.connect.fit, aes(x = Temperature, y = fit), color = 'red') +
      labs(title = 'Raw Curve with Spline Overlay')
    ### final baseline subtracted sample
    g.final <- baseline.sample %>% ggplot(aes(x = Temperature, y = dCp)) + geom_point() +
      labs(title = 'Baseline Subtracted Sample')
    ### overlaid output
    print(cowplot::plot_grid(g.final, g.spline, nrow=2))
  }
  
  return(baseline.sample)
}

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

moving.window <- function(x, w = 90, exclusion.lwr = 60, exclusion.upr = 80, point.selection = "outmost")
{
  if(point.selection == "outmost"){
    k <- w
  }
  if(point.selection == "innermost"){
    k <- 0
  }
  if(point.selection == "mid"){
    k <- round(w/2)
  }
  
  ### fit a CV spline
  full.spline.fit <- smooth.spline(x$Temperature, x$dCp, cv=TRUE)
  ### spline residuals into a data.frame with ids for tracking.
  r <- resid(full.spline.fit)
  r.df <- data.frame(Temperature = working.sample$Temperature, 
                     r=r, 
                     id = 1:length(working.sample$Temperature))
  ### scan through lower region calculating variance within a window size of w
  cat('Scanning Lower. \n')
  i=0
  df.var <- data.frame()
  ### how far do we need to scan?
  points.in.lower <- nrow(working.sample %>% filter(Temperature < exclusion.lwr))
  ### scan and calculate variance for each window
  while ((w+i) < points.in.lower){
    rout <- r.df %>% slice((1+i):(w+i)) %>% summarise(temp.stop=(w+i),mean = mean(r), sd=sd(r))
    df.var <- rbind(df.var, rout)
    i=i+1
  }
  ### where was the minimum variance?
  low <- df.var %>% filter(sd == min(sd))
  lower <- working.sample$Temperature[low$temp.stop-k+1]
  
  ### total size of thermogram
  cat('Scanning Upper. \n')
  l <- length(working.sample$Temperature)
  j=0
  df.var.upper <- data.frame()
  ### how far from upper endpoint do we need to scan?
  smallest.point.in.upr <- nrow(working.sample) - nrow(working.sample %>% filter(Temperature > exclusion.upr))
  ### scan from highest point to region of exclusion
  while ((l-w-j) > smallest.point.in.upr) {
    rout <- r.df %>% slice((l-w-j):(l-j)) %>% summarise(i=(l-w-j),mean = mean(r), sd=sd(r))
    df.var.upper <- rbind(df.var.upper, rout)
    j=j+1
  }
  ### where was the minimum variance?
  up <- df.var.upper %>% filter(sd == min(sd))
  upper <- working.sample$Temperature[up$i+k]
  
  output <- data.frame(lower = lower, upper = upper, method = point.selection)
  return(output)
}

auto.baseline <- function(x, w = 90, exclusion.lwr = 60, exclusion.upr = 80,
                          grid.temp = seq(45, 90, 0.1), plot.on = FALSE, 
                          point = "outmost")
{
  ### automate selection of endpoints
  endpoints <- moving.window(
    x = working.sample,
    w = 90,
    exclusion.lwr = 60,
    exclusion.upr = 80,
    point.selection = point)
  ### baseline subtraction with auto-selected upr/lwr points
  baseline.output <- baseline.subtraction.byhand(
    x = working.sample, 
    lwr.temp = endpoints$lower, 
    upr.temp = endpoints$upper,
    plot.on = plot.on)
  ### generate a final sample on chosen grid!
  final.sample <- final.sample.interpolate(
    x = baseline.output, 
    grid.temp = seq(45, 90, 0.1),
    plot.on = plot.on)
  ### return the interpolated baseline-subtracted result
  return(final.sample)
}