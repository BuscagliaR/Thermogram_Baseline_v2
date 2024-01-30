library(tidyverse)
library(knitr)
library(forecast)
library(gridExtra)

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
    x = x,
    w = w,
    exclusion.lwr = exclusion.lwr,
    exclusion.upr = exclusion.upr,
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

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

### Run the data 

plasma_corrected_data <- readxl::read_excel("data_raw/Baseline.xlsx", 
                                          sheet = "Sheet1")
plasma_raw_data <- readxl::read_excel("data_raw/Normalized data for Robert 12.11.2023.xlsx", 
                                            sheet = "Sheet1")

Plasma.Working <- NULL
Total.Samples <- ncol(plasma_raw_data)/2
for(j in 1:Total.Samples)
{
  lwr <- 2*j - 1
  upr <- 2*j
  if(j %% 20 == 0) cat(j, ' of ', Total.Samples, 'completed. \n')
  temp.col <- plasma_raw_data %>% select(lwr:upr)
  temp.col <- temp.col %>% mutate(SampleID = colnames(temp.col)[2])
  colnames(temp.col)[1:2] <- c('Temperature', 'dCp')
  Plasma.Working <- Plasma.Working %>% rbind(temp.col)
}

Plasma.Working <- Plasma.Working %>% 
  filter(!is.na(Temperature)) %>%
  mutate(SampleNumber = factor(str_extract(SampleID, '\\d+')),
         SampleIteration = factor(str_extract(SampleID, '[a-f]'))) %>%
  mutate(SampleID = factor(SampleID)) %>%
  relocate(SampleID)

final.plasma <- Plasma.Working %>% filter(between(Temperature, 45, 90))

All.IDs <- final.plasma %>% pull(SampleID) %>% unique() %>% as.vector()

n.samples <- length(All.IDs)

### Final temperature grid
grid.temp <- seq(45, 90, 0.1)

### Setup frame to store all data generated
Plasma.Final.Data <- data.frame(Temperature = grid.temp)

### loop over all samples and store the dCp with chosen temperature grid
for(j in 1:n.samples)
{
  cat('Working on Sample ', All.IDs[j], 'element ', j,' of', n.samples,' \n')
  ### select a sample
  working.sample <- final.plasma %>% 
    filter(SampleID == All.IDs[j]) %>% 
    select(Temperature, dCp)
  ### get a baseline-subtracted and interpolated final result!
  auto.output <- auto.baseline(x = working.sample, grid.temp = grid.temp,exclusion.lwr = 48)
  Plasma.Final.Data <- Plasma.Final.Data %>% cbind(out = auto.output$dCp)
  cat("\014")
}

colnames(Plasma.Final.Data)[-1] <- All.IDs


g1 <- Plasma.Final.Data %>% select(Temperature, All.IDs[10]) %>% 
  ggplot(aes(x = Temperature, y = .[,2])) + 
  geom_line() + 
  labs(title = paste0('Final Automated Sample for ', str_sub(All.IDs[10], 2)))
print(g1)

multiplot(g1,g2)

thing <- Plasma.Final.Data %>% pivot_longer(names_to="SampleID", values_to = "dCp", 2:73)

thing <- thing %>% mutate(SampleNumber = str_extract(SampleID, '\\d+'), SampleIteration = str_extract(SampleID, '\\D$'))

numbers <- thing %>% pull(SampleNumber) %>% unique()

other.thing <- Plasma.Working.Corrected %>% filter(between(Temperature, 45, 90))

pdf('generated_output/Plasma_Final_Auto_Graphics_Eight.pdf')
{
g<-NULL
for(i in 1:length(numbers)){
  sample1 <- thing %>% filter(SampleNumber == numbers[i]) %>% mutate(method = str_c("Auto",SampleIteration,sep=" ")) %>% mutate(type = "Auto")
  sample2 <- other.thing %>% filter(SampleNumber == numbers[i]) %>% mutate(method = str_c("Tech",SampleIteration,sep=" ")) %>% mutate(type = "Not Auto")
  
  combinded <- sample1 %>% rbind(sample2)
  
  g[[i]]<- ggplot(combinded, aes(Temperature, dCp, color = method, linetype = type))+
    geom_line()+
    scale_color_manual(breaks = c("Auto a", "Auto b", "Tech a", "Tech b"),
                       values=c("green1", "red1", "green4", "red4"))+
    labs(title = paste0('Final Automated Sample for ', str_sub(numbers[i], 1)))+
    theme(legend.position = "top",
          legend.key.height = unit(0.10, 'cm'),
          legend.key.width = unit(0.30, 'cm'),
          legend.box.spacing = unit(0.20, 'cm'),
          legend.margin = margin(1, 1, 1, 1),
          legend.title = element_text(size=5),
          legend.text=element_text(size=5),
          title = element_text(size = 8))
}
  marrangeGrob(g,nrow=4,ncol=2)
}
dev.off()


save.image('data/Plasma_Working.RData')
load('data/Plasma_Working.RData')

new_thing <- thing
other_thing2 <- other.thing

new_thing <- new_thing %>% mutate(new_title = str_c(SampleID,"Auto",sep="_"))
new_thing2 <- other_thing2 %>% mutate(new_title = str_c(SampleID,"Tech",sep="_"))

new_thing_wide <- new_thing[c("Temperature", "new_title", "dCp")] %>% pivot_wider(names_from = new_title, values_from = dCp)
new_thing_wide_2 <- new_thing2[c("new_title", "dCp")] %>% pivot_wider(names_from = new_title, values_from = dCp)

data.set <- data.frame(Temperature = grid.temp)

for(i in 1:n.samples){
  working.sample <- other.thing %>% 
    filter(SampleID == All.IDs[i]) %>% 
    select(Temperature, dCp)
  
  inter_tech <- final.sample.interpolate(working.sample,grid.temp,plot.on = FALSE)
  
  data.set <- data.set %>% cbind(inter_tech[2])
  print(i)
}

colnames(data.set)[-1] <- str_c(All.IDs,"Tech",sep="_")


final.data.table <- inner_join(data.set,new_thing_wide,by="Temperature")

write.csv(x = final.data.table, file = 'generated_output/Plasma_Data_comparison.csv')

endpoints <- data.frame()
for(j in 1:n.samples)
{
  ### select a sample
  working.sample <- final.plasma %>% 
    filter(SampleID == All.IDs[j]) %>% 
    select(Temperature, dCp)
  ### get a baseline-subtracted and interpolated final result!
  auto.output <- moving.window(x = working.sample)
  endpoints <- endpoints %>% rbind(data.frame(sample = All.IDs[j],lower = auto.output[1], upper=auto.output[2]))
}

write.csv(x = endpoints, file = 'generated_output/Plasma_Endpoints.csv')
