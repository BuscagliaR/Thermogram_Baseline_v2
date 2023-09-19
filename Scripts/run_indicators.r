#' Function that determines if the signal moves up or down
#' Requires being provided first difference curve (i.e. the deltas)
#' 
#'
#' @param x The first difference curve
#' @return returns -1 when x<=0. Returns 1 when x>0
#' @export
run.indicators <- function(x) return(ifelse(x <= 0, -1, 1))