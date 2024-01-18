library(tidyverse)
library(knitr)
library(forecast)
library(ggplot2)
library(mosaic)
load('data/Urine_Working.RData')

first_sample <- Urine.Working.Final %>% filter(SampleID == All_Urine_ID[781])

thing1 <- smooth.spline(first_sample[,2:3],cv=TRUE)

plot(thing1)



poss <- predict(thing1, deriv = 1)
plot(poss)
abline(h=0,type=2)

#poss2 <- predict(thing1, deriv = 2)
#plot(poss2)

der_spline <- smooth.spline(poss, cv=TRUE)
plot(der_spline$fit)

findZeros(der_spline$fit~0,near = 60)

funspline <- splinefun(poss)

findZeros(funspline(x)~x,near = 60, nearest = 20)






pdf('generated_output/graph_deravitives.pdf')
{
for(i in 1:length(All_Urine_ID)){
  sample <- Urine.Working.Final %>% filter(SampleID==All_Urine_ID[i])
  thing <- sample %>% select(Temperature,dCp)
  spline1 <- smooth.spline(thing,cv=TRUE)
  poss <- predict(spline1, deriv = 1)
  g1 <- ggplot(thing, aes(x = Temperature, y = dCp))+
    geom_point()+
    labs(title = paste0('Sample with derivative ', str_sub(All_Urine_ID[i], 1)))
  g2 <- ggplot(data.frame(poss), aes(x = x, y = y))+
    geom_point()
  multiplot(g1,g2)
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
