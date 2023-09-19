#' moving window
#' 
#' @use
#'  Takes a raw thermogram and tries to estimate where the upper and lower
#' cutoffs for baseline should be. Not an easy task! This function can
#' be slow when scanning over large windows or using a small window size (w)
#' 
#' @require dplyr
#'
#' @param x  raw thermogram as data.frame with Temperature and dCp columns.
#' @param w  window size (default 90, try others!)
#' @param exclusion.lwr lower temperature exclusion point
#' @param exclusion.upr upper temperature exclusion point
#' @return an upper and a lower endpoints in data.frame format
#' @export

moving.window <- function(x, w = 90, exclusion.lwr = 60, exclusion.upr = 80)
{
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
  points.in.lower <- nrow(working.sample %>% filter(Temperature < exclude.lwr))
  ### scan and calculate variance for each window
  while ((w+i) < points.in.lower){
    rout <- r.df %>% slice((1+i):(w+i)) %>% summarise(temp.stop=(w+i),mean = mean(r), sd=sd(r))
    df.var <- rbind(df.var, rout)
    i=i+1
  }
  ### where was the minimum variance?
  low <- df.var %>% filter(sd == min(sd))
  lower <- working.sample$Temperature[low$temp.stop-w]
  ### total size of thermogram
  cat('Scanning Upper. \n')
  l <- length(working.sample$Temperature)
  j=0
  df.var.upper <- data.frame()
  ### how far from upper endpoint do we need to scan?
  smallest.point.in.upr <- nrow(working.sample) - nrow(working.sample %>% filter(Temperature > exclude.upr))
  ### scan from highest point to region of exclusion
  while ((l-w-j) > smallest.point.in.upr) {
    rout <- r.df %>% slice((l-w-j):(l-j)) %>% summarise(i=(l-w-j),mean = mean(r), sd=sd(r))
    df.var.upper <- rbind(df.var.upper, rout)
    j=j+1
  }
  ### where was the minimum variance?
  up <- df.var.upper %>% filter(sd == min(sd))
  upper <- working.sample$Temperature[up$i+w]
  output <- data.frame(lower = lower, upper = upper)
  return(output)
}