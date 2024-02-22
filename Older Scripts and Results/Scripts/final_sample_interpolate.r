#' final sample interpolate
#' 
#' @use
#' Takes a baseline subtracted sampled and produces an interpolated result
#' on a chosen temperature grid. Must use baseline.subtraction.byhand() for x
#' 
#' @require dplyr
#' @require cowplot
#'
#' @param x  baseline subtracted data.frame (use baseline function above!)
#' @param grid.temp  grid of desired temperatures
#' @param upr.temp  upper cutoff temperature
#' @param plot.on  outputs a graphic of interpolated sample
#' @return a tibble with 2 columns: Temperature and dCp
#' @export

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