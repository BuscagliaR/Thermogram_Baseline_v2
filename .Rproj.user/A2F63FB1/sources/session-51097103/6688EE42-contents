#' Does the full algorithm
#'
#' Determines the two endpoints for baseline subtraction, then does the baseline subtraction,
#' Then interpolates to the specified temperature grid
#'
#' @param x
#' @param w The number of points in the window
#' @param exclusion.lwr The lower bound of the exclusion window
#' @param exclusion.upr The upper bound of the exclusion window
#' @param grid.temp The set of temperatures to interpolate too
#' @param plot.on Should graphs of the samples be generated. Default is FALSE
#' @param point The method of selecting the endpoint. Default is "outmost"
#' @return data frame with the sample interpolated and the baseline subtracted
#' @export
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