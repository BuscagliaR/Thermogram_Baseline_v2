#' Determines the two endpoints for baseline subtraction
#'
#' 
#'
#' @param x
#' @param w The number of points in the window
#' @param exclusion.lwr The lower bound of the exclusion window
#' @param exclusion.upr The upper bound of the exclusion window
#' @param point.selection The method of selecting the endpoint. Default is "outmost"
#' @return data frame containing the two endpoints for selection
#' @export
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