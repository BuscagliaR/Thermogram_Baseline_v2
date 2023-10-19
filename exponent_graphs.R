library(tidyverse)
library(knitr)
library(forecast)
library(ggplot2)
load('data/Urine_Working.RData')

#working <- Urine.Working.Final %>% filter(SampleID==All_Urine_ID[1])

#Our cutting parameter
t <- 4

pdf('generated_output/Urine_Raw_Data_Exponent_Cut_trial_2.pdf')
{
for(i in 1:length(All_Urine_ID)){
  #set up data
  working <- Urine.Working.Final %>% filter(SampleID==All_Urine_ID[i])
  if(min(working$dCp) < 0){
    shifted_dCp = working$dCp + abs(min(working$dCp))
  }
  else{
    shifted_dCp = working$dCp - abs(min(working$dCp))
  }
  working <- working %>% mutate(shifted_dCp = shifted_dCp)
  working <- working %>% mutate(edCp = exp(shifted_dCp))
  number_points <- working %>% summarise(count = n())
  t <- 4
  
  #Find our potential endpoints
  x = 0
  while(x == 0){
  temps <- working %>% filter(edCp >= max(edCp)/t) %>% summarise(min = min(Temperature),
                                                                 max = max(Temperature)) %>% as.vector()
  
  working <- working %>% mutate(edCp = ifelse(between(Temperature,temps$min[1], temps$max[1]), 0,edCp))
  
  #possible_points <- working %>% filter(edCp >= max(edCp)/t) %>%  summarise(min = min(Temperature),
                                                                                #max = max(Temperature))%>% as.vector()
  leftover_points <- working %>% filter(edCp >= max(edCp)/t) %>% summarise(count = n())
  
  non_zero <- working %>% filter(edCp != 0) %>% summarise(count = n())
  
  if(leftover_points$count[1] < 0.20*non_zero$count[1]){ #0.25*number_points$count[1]
    t = t #t*2
  }
  else{
    x = 1
  }
}
  thing <- working %>% filter(edCp == 0) %>% summarise(min = min(Temperature),
                                                                max = max(Temperature))
  g <- ggplot(working, aes(x = Temperature, y = edCp))+
    geom_line()+
    labs(title = paste0('Exponentiated with Zero Replaced of Sample ', str_sub(All_Urine_ID[i], 1)))
  g1 <- ggplot(working, aes(x = Temperature, y= shifted_dCp))+
    geom_line()
  multiplot(g,g1)
}
}
dev.off() 
  
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

head(working)

working <- working %>% filter(edCp >= max(edCp)/t) %>%  mutate(edCp = exp(dCp))

head(working)

ggplot(working, aes(x = Temperature, y = edCp))+
  geom_line()

half_max <- working %>% filter(edCp == max(edCp)) %>% pull(edCp)/2

temps <- working %>% filter(edCp >= max(edCp)/2) %>% summarise(min = min(Temperature),
                                                      max = max(Temperature)) %>% as.vector()

working.2 <- working %>% mutate(edCp = ifelse(between(Temperature,temps$min[1], temps$max[1]), 0,edCp))


ggplot(working.2, aes(x = Temperature, y = edCp))+
  geom_line()

temps.2 <- working.2 %>% filter(edCp >= max(edCp)/2) %>% summarise(min = min(Temperature),
                                                               max = max(Temperature)) %>% as.vector()

working.3 <- working.2 %>% filter(!between(Temperature,temps.2$min[1], temps.2$max[1]))

ggplot(working.3, aes(x = Temperature, y = edCp))+
  geom_line()

temps.3 <- working.3 %>% filter(edCp >= max(edCp)/2) %>% summarise(min = min(Temperature),
                                                                   max = max(Temperature)) %>% as.vector()

working.4 <- working.3 %>% filter(!between(Temperature,temps.3$min[1], temps.3$max[1]))

ggplot(working.4, aes(x = Temperature, y = edCp))+
  geom_line()


temps.4 <- working.4 %>% filter(edCp >= max(edCp)/2) %>% summarise(min = min(Temperature),
                                                                   max = max(Temperature)) %>% as.vector()

working.5 <- working.4 %>% filter(!between(Temperature,temps.4$min[1], temps.4$max[1]))

ggplot(working.5, aes(x = Temperature, y = edCp))+
  geom_line()


temps.5 <- working.5 %>% filter(edCp >= max(edCp)/2) %>% summarise(min = min(Temperature),
                                                                   max = max(Temperature)) %>% as.vector()

working.6 <- working.5 %>% filter(!between(Temperature,temps.5$min[1], temps.5$max[1]))

ggplot(working.6, aes(x = Temperature, y = edCp))+
  geom_line()



temps.6 <- working.6 %>% filter(edCp >= max(edCp)/2) %>% summarise(min = min(Temperature),
                                                                   max = max(Temperature)) %>% as.vector()

working.7 <- working.6 %>% filter(!between(Temperature,temps.6$min[1], temps.6$max[1]))

ggplot(working.7, aes(x = Temperature, y = edCp))+
  geom_line()


temps.7 <- working.7 %>% filter(edCp >= max(edCp)/2) %>% summarise(min = min(Temperature),
                                                                   max = max(Temperature)) %>% as.vector()

working.8 <- working.7 %>% filter(!between(Temperature,temps.7$min[1], temps.7$max[1]))

ggplot(working.8, aes(x = Temperature, y = edCp))+
  geom_line()


temps.8 <- working.8 %>% filter(edCp >= max(edCp)/2) %>% summarise(min = min(Temperature),
                                                                   max = max(Temperature)) %>% as.vector()

working.9 <- working.8 %>% filter(!between(Temperature,temps.8$min[1], temps.8$max[1]))

ggplot(working.9, aes(x = Temperature, y = edCp))+
  geom_line()

temps.9 <- working.9 %>% filter(edCp >= max(edCp)/2) %>% summarise(min = min(Temperature),
                                                                   max = max(Temperature)) %>% as.vector()

working.10 <- working.9 %>% filter(!between(Temperature,temps.9$min[1], temps.9$max[1]))

ggplot(working.10, aes(x = Temperature, y = edCp))+
  geom_line()


temps.10 <- working.10 %>% filter(edCp >= max(edCp)/1000) %>% summarise(min = min(Temperature),
                                                                   max = max(Temperature)) %>% as.vector()

working.11 <- working.10 %>% filter(!between(Temperature,temps.10$min[1], temps.10$max[1]))

ggplot(working.11, aes(x = Temperature, y = edCp))+
  geom_line()


temps.10 <- working.10 %>% filter(edCp >= max(edCp)/5) %>% summarise(min = min(Temperature),
                                                                     max = max(Temperature)) %>% as.vector()

working.11 <- working.10 %>% filter(!between(Temperature,temps.10$min[1], temps.10$max[1]))

ggplot(working.11, aes(x = Temperature, y = edCp))+
  geom_line()
