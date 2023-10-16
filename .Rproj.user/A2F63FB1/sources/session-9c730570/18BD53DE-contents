library(tidyverse)
library(knitr)
library(forecast)
library(ggplot2)
load('data/Urine_Working.RData')

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



All.IDs <- Urine.Working.Final %>% pull(SampleID) %>% unique() %>% as.vector()
n.samples <- length(All.IDs)

### Final temperature grid
grid.temp <- seq(45, 90, 0.1)

### Setup frame to store all data generated
Urine.Final.Data <- data.frame(Temperature = grid.temp)

### loop over all samples and store the dCp with chosen temperature grid
for(j in 1:n.samples)
{
  cat('Working on Sample ', All.IDs[j], 'element ', j,' of', n.samples,' \n')
  ### select a sample
  working.sample <- Urine.Working.Final %>% 
    filter(SampleID == All.IDs[j]) %>% 
    select(Temperature, dCp)
  ### get a baseline-subtracted and interpolated final result!
  auto.output <- auto.baseline(x = working.sample, grid.temp = grid.temp)
  Urine.Final.Data <- Urine.Final.Data %>% cbind(out = auto.output$dCp)
  cat("\014")
}

table.outmost <- data.frame()
for(i in 1:n.samples){
  working.sample <- Urine.Working.Final %>% 
    filter(SampleID == All.IDs[i]) %>% 
    select(Temperature, dCp)
  thing <- moving.window(x = working.sample)
  thing <- thing %>% cbind(SampleID = All.IDs[i])
  table.outmost <- table.outmost %>% rbind(thing)
  print(i)
}


table.innermost <- data.frame()
for(i in 1:n.samples){
  working.sample <- Urine.Working.Final %>% 
    filter(SampleID == All.IDs[i]) %>% 
    select(Temperature, dCp)
  thing <- moving.window(x = working.sample, point.selection = "innermost")
  thing <- thing %>% cbind(SampleID = All.IDs[i])
  table.innermost <- table.innermost %>% rbind(thing)
  print(i)
}

table.mid <- data.frame()
for(i in 1:n.samples){
  working.sample <- Urine.Working.Final %>% 
    filter(SampleID == All.IDs[i]) %>% 
    select(Temperature, dCp)
  thing <- moving.window(x = working.sample, point.selection = "mid")
  thing <- thing %>% cbind(SampleID = All.IDs[i])
  table.mid <- table.mid %>% rbind(thing)
  print(i)
}

head(table.mid)

final.table <- data.frame()
final.table <- final.table %>% rbind(table.innermost) %>% rbind(table.mid) %>% rbind(table.outmost)
copy1 <- final.table


grid.temp <- seq(45, 90, 0.1)

### Setup frame to store all data generated
thing1 <- data.frame(Temperature = grid.temp)

for(i in 1:n.samples){
    working.sample <- Urine.Working.Final %>% 
      filter(SampleID == All.IDs[i]) %>% 
      select(Temperature, dCp)
    endpoints <- copy1 %>% filter(SampleID == All.IDs[i]) %>% filter(method == "innermost")
    baseline.output <- baseline.subtraction.byhand(
      x = working.sample, 
      lwr.temp = endpoints$lower, 
      upr.temp = endpoints$upper,
      plot.on = FALSE)
    thing <- final.sample.interpolate(
      x = baseline.output, 
      grid.temp = seq(45, 90, 0.1),
      plot.on = FALSE)
    thing1 <- thing1 %>% cbind(out = thing$dCp)
}

thing2 <- data.frame(Temperature = grid.temp)

for(i in 1:n.samples){
  working.sample <- Urine.Working.Final %>% 
    filter(SampleID == All.IDs[i]) %>% 
    select(Temperature, dCp)
  endpoints <- copy1 %>% filter(SampleID == All.IDs[i]) %>% filter(method == "outmost")
  baseline.output <- baseline.subtraction.byhand(
    x = working.sample, 
    lwr.temp = endpoints$lower, 
    upr.temp = endpoints$upper,
    plot.on = FALSE)
  thing <- final.sample.interpolate(
    x = baseline.output, 
    grid.temp = seq(45, 90, 0.1),
    plot.on = FALSE)
  thing2 <- thing2 %>% cbind(out = thing$dCp)
}

thing3 <- data.frame(Temperature = grid.temp)

for(i in 1:n.samples){
  working.sample <- Urine.Working.Final %>% 
    filter(SampleID == All.IDs[i]) %>% 
    select(Temperature, dCp)
  endpoints <- copy1 %>% filter(SampleID == All.IDs[i]) %>% filter(method == "mid")
  baseline.output <- baseline.subtraction.byhand(
    x = working.sample, 
    lwr.temp = endpoints$lower, 
    upr.temp = endpoints$upper,
    plot.on = FALSE)
  thing <- final.sample.interpolate(
    x = baseline.output, 
    grid.temp = seq(45, 90, 0.1),
    plot.on = FALSE)
  thing3 <- thing3 %>% cbind(out = thing$dCp)
}

colnames(thing1)[-1] <- All.IDs
Sample_Names <- colnames(thing1)[-1]

colnames(thing2)[-1] <- All.IDs
Sample_Names <- colnames(thing2)[-1]

colnames(thing3)[-1] <- All.IDs
Sample_Names <- colnames(thing3)[-1]

one.sample <- thing1 %>% select(Temperature,Sample_Names[1] )

g1 <- thing1 %>% select(Temperature, Sample_Names[1])
g1 <- g1 %>% mutate(method = "Inner")
g2 <- thing2 %>% select(Temperature, Sample_Names[1])
g2<- g2 %>% mutate(method = "Outer")
g3 <- thing3 %>% select(Temperature, Sample_Names[1])
g3<- g3 %>% mutate(method = "Mid")

totalframe <- rbind(g1,g2,g3)

ggplot(totalframe, aes(x = Temperature,y = totalframe[,2], color=method))+
  geom_line()


pdf('generated_output/Plots_of_Varying_point_selection_with_names.pdf')
{
gs <- NULL
for(i in 1:n.samples)
{
  g1 <- thing1 %>% select(Temperature, Sample_Names[i])
  g1 <- g1 %>% mutate(method = "Inner")
  g2 <- thing2 %>% select(Temperature, Sample_Names[i])
  g2<- g2 %>% mutate(method = "Outer")
  g3 <- thing3 %>% select(Temperature, Sample_Names[i])
  g3<- g3 %>% mutate(method = "Mid")
  
  totalframe <- rbind(g1,g2,g3)
  
  gs <- ggplot(totalframe, aes(x = Temperature,y = totalframe[,2], color=method))+
    geom_point(size = 0.5)+
    labs(title = paste0('Final Automated Sample for ', str_sub(All_Urine_ID[i], 1)))
  print(gs)
}
}
dev.off()



write.csv(x = totalframe, file = 'generated_output/graphable_Table.csv')

totalframe <- Endpoint_Table

pdf('generated_output/Plots_of_Varying_point_selection_with_names.pdf')
{
  gs <- NULL
  for(i in 1:n.samples){
    gs <- ggplot(totalframe, aes(x = Temperature,y = totalframe[,2], color=method))+
      geom_point(size = 0.5)+
      labs(title = paste0('Final Automated Sample for ', str_sub(All_Urine_ID[i], 1)))
    print(gs)
  }
}
dev.off()
