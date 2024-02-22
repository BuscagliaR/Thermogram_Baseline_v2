#' Interpolates each sample on a fixed grid of temperature
#'
#' Fits a spline to the baseline subtracted data and interpolates the thermogram on the given temperature grid.
#'
#' @param x Baseline subtracted thermogram data
#' @param grid.temp The grid of temperatures the sample is interpolated onto
#' @param plot.on logical: should the output be graphed?
#' @return Data frame of temperature and dCp
#' @export
final.sample.interpolate <- function(x, grid.temp, plot.on = TRUE)
{
  require(ggplot2)
  require(dplyr)
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
