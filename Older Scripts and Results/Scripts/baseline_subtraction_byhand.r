#' baseline subtraction by hand
#' 
#' @use
#' Function for by-hand baseline subtraction
#' This can be used to introduce some automation to the temperature choices
#' but does provide a function in R for doing baseline subtraction!
#' 
#' @require dplyr
#' @require cowplot
#'
#' @param x data.frame with two columns: Temperature and dCp
#' @param lwr.temp lower cutoff temperature
#' @param upr.temp  upper cutoff temperature
#' @param plot.on  outputs a graphic of final sample + baseline procedure when TRUE
#' @return a tibble with 2 columns: Temperature and dCp but with baseline subtracted
#' @export

baseline.subtraction.byhand <- function(x, lwr.temp, upr.temp, plot.on = TRUE)
{
  ### check-conditions of boundaires - this effects automation
  if(lwr.temp < min(x$Temperature)+1) lwr.temp = lwr.temp + 1
  if(upr.temp > max(x$Temperature)-1) upr.temp = upr.temp - 1
  ### Extract the baseline regions
  work.lower <- x %>% filter(Temperature < lwr.temp)
  work.upper <- x %>% filter(Temperature > upr.temp)
  ### Splines for lower/upper regions
  spline.lower <- smooth.spline(work.lower$Temperature, work.lower$dCp, cv = TRUE)
  spline.upper <- smooth.spline(work.upper$Temperature, work.upper$dCp, cv = TRUE)
  ### Store data for graphing
  spline.lower.fit <- data.frame(Temperature = work.lower$Temperature, fit = spline.lower$y)
  spline.upper.fit <- data.frame(Temperature = work.upper$Temperature, fit = spline.upper$y)
  ### store middle (signal) region
  work.mid <- x %>% filter(between(Temperature, lwr.temp, upr.temp))
  ### find endpoints of splines
  spline.connect.points <- rbind(
    spline.lower.fit %>% filter(Temperature == max(Temperature)),
    spline.upper.fit %>% filter(Temperature == min(Temperature)))
  ### connect endpoints and store
  spline.connect.lm <- lm(fit ~ Temperature, data = spline.connect.points)
  spline.connect.fit <- data.frame(
    Temperature = work.mid$Temperature,
    fit = predict(spline.connect.lm, data.frame(Temperature = work.mid$Temperature)))
  ### store baseline as one unit
  working.baseline.final <- rbind(spline.lower.fit, spline.connect.fit, spline.upper.fit)
  ### join for tidyverse simplification
  baseline.join <- full_join(x, working.baseline.final, by = 'Temperature')
  ### final sample!
  baseline.sample <- baseline.join %>% mutate(final.dcp = dCp - fit) %>%
    select(Temperature, final.dcp) %>% rename(dCp = final.dcp)
  if(plot.on)
  {
    ### graph of raw with spline
    g.spline <- working.sample %>% ggplot(aes(x = Temperature, y = dCp)) + geom_point() +
      geom_line(data = spline.lower.fit, aes(x = Temperature, y = fit), color = 'red') +
      geom_line(data = spline.upper.fit, aes(x = Temperature, y = fit), color = 'red') +
      geom_line(data = spline.connect.fit, aes(x = Temperature, y = fit), color = 'red') +
      labs(title = 'Raw Curve with Spline Overlay')
    ### final baseline subtracted sample
    g.final <- baseline.sample %>% ggplot(aes(x = Temperature, y = dCp)) + geom_point() +
      labs(title = 'Baseline Subtracted Sample')
    ### overlaid output
    print(cowplot::plot_grid(g.final, g.spline, nrow=2))
  }
  return(baseline.sample)
}
