#' automated baseline function
#' 
#' @use
#' Automatically determines a baseline for RAW thermograms and returns
#' baseline subtracted set on a chosen temperature grid!
#' 
#' @require dplyr
#'
#' @param x  raw thermogram as data.frame with Temperature and dCp columns.
#' @param w  window size (default 90, try others!)
#' @param exclusion.lwr lower temperature exclusion point
#' @param exclusion.upr upper temperature exclusion point
#' @param grid.temp  chosen temperature grid for final data
#' 
#' @return return 
#' @export

auto.baseline <- function(x, w = 90, exclusion.lwr = 60, exclusion.upr = 80,
                          grid.temp = seq(45, 90, 0.1), plot.on = FALSE)
{
  ### automate selection of endpoints
  endpoints <- moving.window(
    x = working.sample,
    w = 90,
    exclusion.lwr = 60,
    exclusion.upr = 80)
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