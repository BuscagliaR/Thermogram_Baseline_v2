#' Function that counts runs of ups/downs
#' Relies run.indicators() to count how many consecutive ups/downs
#' 
#' returns two column data frame:
#' RunSize  with run of same monotonic changes (Ex. 1, 1, 1, -1 is 3) 
#' Up/Down whether the current run of monotonic changes is positive or negative
#' (Up being positive)
#'
#' @param x 
#' @return data frame with 2 columns
#' @export

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