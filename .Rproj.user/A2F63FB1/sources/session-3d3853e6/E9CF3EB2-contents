#' Interpolates each sample on a fixed grid of temperature
#'
#' 
#'
#' @param x
#' @param grid.temp The grid of temperatures the sample is interpolated onto
#' @param plot.on
#' @return data frame of temperature and dCp with the sample interpolated to the grid temperatureds
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