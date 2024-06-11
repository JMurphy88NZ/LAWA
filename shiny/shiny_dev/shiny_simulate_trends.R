


# This script consists of  functions that simulate new series to replace existing trend components
library(dplyr)
library(tsibble)
library(imputeTS)
library(fable)
library(tibble)
library(purrr)
library(feasts)



#' Generate a Cosine Series with Decay
#'
#' This function generates a cosine series with an exponential decay in amplitude.
#'
#' @param total_length Integer. The total length of the time series.
#' @param initial_amplitude Numeric. The initial amplitude of the cosine wave.
#' @param decay_rate Numeric. The rate at which the amplitude decays.
#' @param num_peaks Integer. The number of peaks in the cosine wave.
#' @param phase_shift Numeric. The phase shift of the cosine wave.
#'
#' @return A numeric vector representing the cosine series.
#' @examples
#' generate_cosine_series(total_length = 100, initial_amplitude = 1, decay_rate = 0.01, num_peaks = 5, phase_shift = 0)
generate_cosine_series <- function(total_length, initial_amplitude, decay_rate, num_peaks, phase_shift) {
  
  # Time vector
  time <- 1:total_length
  
  # Angular frequency
  angular_frequency <- 2 * pi * num_peaks / total_length
  
  # Amplitude decay function
  amplitude <- initial_amplitude * exp(-decay_rate * time)
  
  # Generate the cosine wave with phase shift
  cosine_wave <- amplitude * cos(angular_frequency * time + phase_shift)
  
  return(cosine_wave)
}

#' Generate a Linear Trend
#'
#' This function generates a linear trend.
#'
#' @param total_length Integer. The total length of the time series.
#' @param slope Numeric. The slope of the linear trend.
#' @param intercept Numeric. The intercept of the linear trend.
#'
#' @return A numeric vector representing the linear trend.
#' @examples
#' generate_linear_trend(total_length = 100, slope = 0.5, intercept = 10)
generate_linear_trend <- function(total_length, slope, intercept) {
  # Time vector
  time <- 1:total_length
  
  # Linear trend
  linear_trend <- slope * time + intercept
  
  return(linear_trend)
}

#' Generate a Level Shift with Ramp
#'
#' This function generates a level shift with a ramp in the time series.
#'
#' @param total_length Integer. The total length of the time series.
#' @param baseline_amp Numeric. The baseline amplitude before the shift. Default is 1.
#' @param amp_change Numeric. The change in amplitude during the ramp, expressed as a percentage of the baseline amplitude. Default is 0.2.
#' @param ramp_start Numeric. The start point of the ramp, expressed as a percentage of the total length. Default is 0.5.
#' @param ramp_length Numeric. The length of the ramp, expressed as a percentage of the total length. Default is 0.05.
#' @param steepness Numeric. The steepness of the ramp. Default is 1.
#'
#' @return A numeric vector representing the time series with a level shift and ramp.
#' @examples
#' level_shift_with_ramp(total_length = 100, baseline_amp = 1, amp_change = 0.2, ramp_start = 0.3, ramp_length = 0.4, steepness = .8)
generate_level_shift_with_ramp <- function(total_length, baseline_amp = 1, amp_change = 0.2, 
                                  ramp_start = 0.5, ramp_length = 0.05, steepness = 1){
  
  # Convert percentages into actual values
  amp_change <- amp_change * baseline_amp  
  ramp_start <- round(ramp_start * total_length, digits = 0)
  ramp_length <- round(ramp_length * total_length, digits = 0)
  
  # Define the time series
  time <- 1:total_length
  
  # Create a placeholder for the entire series
  values <- rep(baseline_amp, total_length)
  
  # Set default values if not provided
  if (is.null(ramp_start)) {
    ramp_start <- floor(total_length / 2) # Default to start the ramp at the midpoint
  }
  if (is.null(ramp_length)) {
    ramp_length <- max(1, floor(total_length * 0.1)) # Default to 10% of total length, with a minimum of 1
  }
  
  # Define the ramp
  ramp_time <- ramp_start:(ramp_start + ramp_length - 1)
  ramp_end_level <- baseline_amp + amp_change
  midpoint <- ramp_length / 2
  scale <- amp_change
  
  ramp_values <- baseline_amp + scale / (1 + exp(-steepness * (ramp_time - ramp_start - midpoint + 1)))
  
  # Place the ramp values into the series
  values[ramp_start:(ramp_start + ramp_length - 1)] <- ramp_values
  
  # Ensure the series ends at the final value after the ramp
  if ((ramp_start + ramp_length) <= total_length) {
    values[(ramp_start + ramp_length):total_length] <- ramp_end_level
  }
  
  return(values)
}

#level_shift_with_ramp(total_length = 50)


#' Generate a gaussian time series
#' 
#' Put the centre near one the end of the timeseries to reproduce the examples from lawa
#'
#' @param total_length Integer. The total length of the time series.
#' @param centre Numeric.Point on which the gaussian is centred
#' @param width Numeric. The width (standard deviation) of the gaussian
#' @param amplitude Numeric. Maximum value of the returned time series
#'
#' @return A numeric vector representing the time series with a gaussian
#' @export
#'
#' @examples
#' generate_gaussian(100, 1, 3, 3)
#' generate_gaussian(100, 1, 97, 3)
generate_gaussian <- function(total_length, amplitude, centre, width) {
  # Time vector
  time <- 1:total_length
  y <- dnorm(time, centre, width)
  # Scale the function, for ease of plotting
  max_y <- max(y)
  min_y <- min(y)
  y <- amplitude*(y - min_y)/(max_y - min_y)
}

#' Title
#'
#' @param total_length Integer. The total length of the time series.
#' @param x_min Numeric Position of curve minimum within the time series
#' @param scale Numeric Scale parameter - for larger values, the rhs doesn't rise much. for small values, it does
#' @param amplitude The maximum value of the returned series
#'
#' @return
#' @export
#'
#' @examples
#' generate_linex(100, 1, 20, 0.15)
#' generate_linex(100, 1, 20, 0.2)
generate_linex <- function(total_length, amplitude, x_min, scale){
  # Time vector
  time <- 1:total_length
  # Linex function
  linex <- exp(- scale * (time - x_min)) + scale*(time - x_min) - 1
  # Scale the function, for ease of plotting
  max_l <- max(linex)
  min_l <- min(linex)
  linex <- amplitude*(linex - min_l)/(max_l - min_l)
}
