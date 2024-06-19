
#load LandwaterPeople (LWP) functions for trend estimation 
#source("For_Andrew/LWPTrends_v2102.R")
source(here::here("LWPTrends_v2102.R"))

#library(tidyverse)
library(ggplot2)
library(withr)
library(feasts)
library(tsibble)

#' Get STL Decomposition
#'
#' This function performs STL decomposition on the input time series data, handling gaps and imputing missing values if necessary.
#'
#' @param data A data frame containing the time series data. The data frame should include columns `yearmon` and `RawValue`.
#' @return A list containing the STL decomposition results and original data.
#' @importFrom tsibble as_tsibble has_gaps fill_gaps
#' @importFrom imputeTS na_interpolation
#' @importFrom fable model STL
#' @export
Get_STL <- function(data){
  
  # Convert data to tsibble for convenience functions, including STL
  data_tsibble <- data %>% 
    tsibble::as_tsibble(index = yearmon, regular = TRUE)
  
  start <- data$myDate[1] # Assumes arranged by date already
  
  # Test for gaps/implicit missingness
  has_gaps_check <- data_tsibble %>%  
    tsibble::has_gaps()
  
  # If gaps present, fill them in with NAs
  if (has_gaps_check == TRUE){ 
    data_tsibble <- data_tsibble %>%  
      tsibble::fill_gaps()
    
    # Convert to ts data
    ts_data <- ts(data_tsibble$RawValue, start = start, frequency = 12) 
    ts_data <- imputeTS::na_interpolation(ts_data) 
    
    data_tsibble$final_series <- ts_data
    
    data_stl <- data_tsibble %>% 
      model(stl = STL(final_series))
    
    # Add tag whether imputed values included in series or not
    data_stl$imputed <- "YES"
    
  } else { # If no gaps, no need for imputation
    
    ts_data <- ts(data_tsibble$RawValue, start = start, frequency = 12) 
    
    data_tsibble$final_series <- ts_data
    
    data_stl <- data_tsibble %>% 
      model(stl = STL(final_series))
    
    # Add tag whether imputed values included in series or not
    data_stl$imputed <- "NO"
  } 
  
  # Add original columns from data as needed for LWP functions
  data_stl$orig_data <- list(data_tsibble)
  #data_stl$orig_data <- list(data %>%  select(lawa_site_id, CenType, Censored, yearmon, Season, Year, myDate, RawValue))
  
  return(data_stl)
}

#test 
#test_stl <- Get_STL(N03N_filtered$`GW-00004`)



#TODO add in option to test for seasonality if not already known
#' Scale STL Noise
#'
#' This function scales the noise in the STL decomposition by a series of lambda values and performs trend analysis.
#'
#' @param stl_data The STL decomposition data from `Get_STL`.
#' @param lambda A numeric vector of lambda values to scale the noise.
#' @param is_seasonal A logical indicating whether to perform seasonal or nonseasonal trend analysis
#' @param analysis_params A list of additional parameters for trend analysis functions.
#' @return A data frame with the estimated results for each lambda value.
#' @importFrom tsibble as_tibble
#' @importFrom dplyr left_join mutate select
#' @importFrom purrr imap_dfr
#' @export
scale_stl_noise <- function(stl_data, lambda = seq(from = -1, to = 3, by = 0.5), 
                            is_seasonal = TRUE, analysis_params = list()) {
  
  # Container lists
  complist <- list()
  senslope_res_list <- list()
  
  comp_df <- components(stl_data) %>% as_tibble()
  
  # Needs extra columns for LWP functions to work
  orig_data <- stl_data$orig_data[[1]] %>% 
    select(lawa_site_id, CenType, Censored, 
           yearmon, Season, Year, myDate, RawValue)
  
  # Add columns so LWP functions work
  comp_df <- comp_df %>% left_join(orig_data, by = "yearmon")
  
  ## Scale noise by lambda here in a for loop (estimating slope each time)
  for (lambda_value in lambda) {
    comp_df <- comp_df %>% 
      mutate(simulated = final_series + (remainder * lambda_value),
             RawValue = simulated)  # Only so named because LWP functions expect that name
    
    if (is_seasonal == TRUE) {
      senslope_res <- do.call(SeasonalTrendAnalysis, c(list(as.data.frame(comp_df)), analysis_params))  # Pass additional arguments
    } else {
      senslope_res <- do.call(NonSeasonalTrendAnalysis, c(list(as.data.frame(comp_df)), analysis_params))  # Pass additional arguments
    }
    
    slope_df <- tibble(est_slope = senslope_res$AnnualSenSlope / 12,
                       lci = senslope_res$Sen_Lci / 12,
                       uci = senslope_res$Sen_Uci / 12,
                       CI_width = uci - lci)
    
    senslope_res_list[[as.character(lambda_value)]] <- slope_df
  }
  
  # Combine each estimate into same tidy dataframe
  senslope_res_list <- imap_dfr(senslope_res_list, ~ tibble(lambda = .y, .x))
  
  estimate_results <- senslope_res_list %>% 
    mutate(CI_width = uci - lci)
  
  return(estimate_results)
}


#test_noise <- scale_stl_noise(test_stl)



#' Analyze Trend with Noise
#'
#' This wrapper function performs STL decomposition, optionally modifies the trend component, scales the noise, and performs trend analysis.
#'
#' @param data A data frame containing the time series data.
#' @param lambda A numeric vector of lambda values to scale the noise.
#' @param is_seasonal A logical indicating whether to perform seasonal or nonseasonal trend analysis.
#' @param trend_params A list of parameters for modifying the trend component, including `modify_trend` and `initial_amplitude` (see simulate_trends_functions.R) .
#' @param analysis_params A list of additional parameters for trend analysis functions.
#' @param ... Additional arguments passed to the trend modification function.
#' @return A list containing the STL decomposition data and the estimated results.
#' @export
analyze_trend_with_noise <- function(data, lambda = seq(from = -1, to  = 3, by =  0.5), 
                                     is_seasonal = TRUE, 
                                     trend_params = list(modify_trend = NULL, initial_amplitude = NULL), 
                                     analysis_params = list(), mod_fun, ...){
  # Step 1: Get STL decomposition
  stl_data <- Get_STL(data)
  
  # Join with metadata
  comp_df <- stl_data$stl[[1]]$fit$decomposition
  
  # Step 2: Optionally modify the trend component
  if (is.null(trend_params$modify_trend)) {
    message("Keeping original trend component")
    estimate_results <- scale_stl_noise(stl_data, lambda = lambda, is_seasonal = is_seasonal, analysis_params = analysis_params)
  } else {
    trend_params$total_length <- nrow(comp_df)  
    message(paste("modify_trend option:", trend_params$modify_trend))
    
    # Remove the modify_trend element from trend_params
    trend_params <- trend_params[setdiff(names(trend_params), "modify_trend")]
    
    # Modify trend component
    comp_df$trend <- do.call(mod_fun, trend_params)
    
    # Update STL data
    comp_df$final_series <- comp_df$trend + comp_df$remainder
    comp_df$season_adjust <- comp_df$final_series - comp_df$season_year
    stl_data$stl[[1]]$fit$decomposition <- comp_df
    
    estimate_results <- scale_stl_noise(stl_data, lambda = lambda, is_seasonal = is_seasonal, analysis_params = analysis_params)
  }
  
  # Step 4: Return the results
  return(list(stl_data = stl_data, estimate_results = estimate_results))
}





analyze_trend_rolling <- function(data, lambda = seq(from = -1, to  = 3, by =  0.5), 
                                     is_seasonal = TRUE, 
                                     trend_params = list(modify_trend = NULL, initial_amplitude = NULL), 
                                     analysis_params = list(), mod_fun, years = 5, ...){
  # Step 1: Get STL decomposition
  stl_data <- Get_STL(data)
  
  # Join with metadata
  comp_df <- stl_data$stl[[1]]$fit$decomposition
  
  # Step 2: Optionally modify the trend component
  if (is.null(trend_params$modify_trend)) {
    message("Keeping original trend component")
    estimate_results <- scale_stl_noise(stl_data, lambda = lambda, is_seasonal = is_seasonal, analysis_params = analysis_params)
  } else {
    trend_params$total_length <- nrow(comp_df)  
    message(paste("modify_trend option:", trend_params$modify_trend))
    
    # Remove the modify_trend element from trend_params
    trend_params <- trend_params[setdiff(names(trend_params), "modify_trend")]
    
    # Modify trend component
    comp_df$trend <- do.call(mod_fun, trend_params)
    
    # Update STL data
    comp_df$final_series <- comp_df$trend + comp_df$remainder
    comp_df$season_adjust <- comp_df$final_series - comp_df$season_year
    
    ###rolling window for   analysis
    
    #add param for full length
    full_length <- nrow(comp_df)/12  
    
    comp_df <- filter_time_series <- function(comp_df, years = full_length)
    

    stl_data$stl[[1]]$fit$decomposition <- comp_df
    
  
    
    estimate_results <- scale_stl_noise(stl_data, lambda = lambda, is_seasonal = is_seasonal, analysis_params = analysis_params)
  }
  
  # Step 4: Return the results
  return(list(stl_data = stl_data, estimate_results = estimate_results))
}





#' Filter Time Series Data for Custom Periods
#'
#' This function filters a time series data frame for custom periods.
#'
#' @param ts_data A data frame containing the time series data.
#' @param period A numeric vector of length 2 specifying the start and end years. The default is c(5, 0).
#'
#' @return A filtered data frame containing the data for the specified period.
#' @export
#'
#' @examples
#' # Example time series data
#' ts_data <- tibble(
#'   yearmon = seq(as.Date("2000-01-01"), by = "month", length.out = 240),
#'   value = rnorm(240)
#' )
#'
#' # Valid period
#' filtered_data <- filter_custom_period(ts_data, c(5, 3))
#' print(filtered_data)
#'
#' # Invalid period: start period longer than total length
#' # This will stop with an error
#' try(filter_custom_period(ts_data, c(25, 5)))
#'
filter_custom_period <- function(ts_data, period = c(5, 0)) {
  
  # Check if period is a numeric vector of length 2
  if (!is.numeric(period) || length(period) != 2) {
    stop("The 'period' argument must be a numeric vector of length 2 specifying the start and end years.")
  }
  
  total_years <- nrow(ts_data) / 12
  start_years <- period[1]
  end_years <- period[2]
  
  # Check if requested start of period is outside time series
  if (start_years > total_years) {
    stop("The start period is greater than the total length of the data.")
  }
  
  start_months <- round(start_years * 12)
  end_months <- round(end_years * 12)
  
  total_months <- nrow(ts_data)
  
  # Calculate the start and end rows for slicing
  start_row <- max(total_months - start_months + 1, 1)
  end_row <- total_months - end_months
  
  # Filter for the given range
  ts_data_filt <- ts_data %>% 
    slice(start_row:end_row)
  
  return(ts_data_filt)
}



#' Plot Period Comparison
#'
#' This function creates a plot comparing the original time series data with the estimated slopes for different periods.
#'
#' @param full_data A data frame containing the full time series data with columns `yearmon` and `final_series`.
#' @param period_df A data frame containing the period-specific slope estimates and nested data for each period. The data frame should have columns `period`, `est_slope`, and `data` (where `data` is a nested data frame containing `yearmon`).
#'
#' @return A ggplot object showing the original time series data and the estimated slopes for different periods.
#' @export
#'
#' @examples
#' # Example usage
#' full_data <- tibble(
#'   yearmon = seq(as.Date("2000-01-01"), by = "month", length.out = 240),
#'   final_series = rnorm(240)
#' )
#'
#' period_df <- tibble(
#'   period = c("period_10_0", "period_10_5", "period_5_0"),
#'   est_slope = c(-0.0003, 0.00167, -0.00252),
#'   data = list(
#'     tibble(yearmon = seq(as.Date("2000-01-01"), by = "month", length.out = 10)),
#'     tibble(yearmon = seq(as.Date("2005-01-01"), by = "month", length.out = 10)),
#'     tibble(yearmon = seq(as.Date("2010-01-01"), by = "month", length.out = 10))
#'   )
#' )
#'
#' plot_period_comparison(full_data, period_df)
plot_period_comparison <- function(full_data, period_df){
  
  orig_plot <- full_data %>% 
    ggplot(aes(x = yearmon, y = final_series))+
    geom_line()
  
  
  slope_df <- period_df %>% 
    tidyr::hoist(data, "yearmon") %>% 
    tidyr::unnest(yearmon)
  
  
  slope_df <- slope_df %>% 
    mutate(period = as.factor(period)) %>% 
    group_by(period ) %>% 
    mutate(rowid = dplyr::row_number(),
           slope_line = est_slope*rowid)
  
  orig_plot <- orig_plot+
    geom_line(data = slope_df,
              aes(x = yearmon, y = slope_line, colour = period, ),show.legend = FALSE)
  
  return(orig_plot)
  
}



#' Rolling Trend Analysis
#'
#' This function performs rolling trend analysis on STL decomposed time series data.
#'
#' @param stl_data A list containing the STL decomposed time series data.
#' @param periods A list specifying the periods for which the trend analysis should be performed. Can include "full_length" and numeric vectors specifying start and end years.
#' @param is_seasonal A logical value indicating whether the data is seasonal. Defaults to TRUE.
#' @param analysis_params A list of additional parameters to be passed to the trend analysis functions.
#'
#' @return A data frame containing the estimated slopes and confidence intervals for each period.
#' @export
#'
#' @examples
#' # Example usage
#' results <- rolling_trend(stl_data, periods = list("full_length", c(5, 0)), is_seasonal = TRUE, analysis_params = list())
#' print(results)

rolling_trend <- function(stl_data, periods = list("full_length", c(5, 0)), 
                          is_seasonal = TRUE, analysis_params = list()) {
  
  # Container lists
  complist <- list()
  senslope_res_list <- list()
  
  comp_df <- components(stl_data) %>% as_tibble()
  
  # Needs extra columns for LWP functions to work
  orig_data <- stl_data$orig_data[[1]] %>% 
    select(lawa_site_id, CenType, Censored, 
           yearmon, Season, Year, myDate, RawValue)
  
  # Add columns so LWP functions work
  comp_df <- comp_df %>% left_join(orig_data, by = "yearmon")
  
  # Iterate over each period and perform trend analysis
  for (period in periods) {
    
    if (is.character(period) && period == "full_length") {
      # Full length analysis without filtering
      if (is_seasonal) {
        senslope_res <- do.call(SeasonalTrendAnalysis, c(list(as.data.frame(comp_df)), analysis_params))
      } else {
        senslope_res <- do.call(NonSeasonalTrendAnalysis, c(list(as.data.frame(comp_df)), analysis_params))
      } 
        
        # Create a data frame with the results
        slope_df <- tibble(est_slope = senslope_res$AnnualSenSlope / 12,
                           lci = senslope_res$Sen_Lci / 12,
                           uci = senslope_res$Sen_Uci / 12,
                           CI_width = uci - lci,
                           data = list(comp_df))
        
        senslope_res_list[[paste0("period_", paste(period, collapse = "_"))]] <- slope_df  
        
        

    } else if (is.numeric(period) && length(period) == 2) {
      # Filter the data for the specified period
      comp_df_filt <- filter_custom_period(comp_df, period)
      
      if (is_seasonal) {
        senslope_res <- do.call(SeasonalTrendAnalysis, c(list(as.data.frame(comp_df_filt)), analysis_params))
      } else {
        senslope_res <- do.call(NonSeasonalTrendAnalysis, c(list(as.data.frame(comp_df_filt)), analysis_params))
      }
      
      # Create a data frame with the results
      slope_df <- tibble(est_slope = senslope_res$AnnualSenSlope / 12,
                         lci = senslope_res$Sen_Lci / 12,
                         uci = senslope_res$Sen_Uci / 12,
                         CI_width = uci - lci,
                         data = list(comp_df_filt))
      
      senslope_res_list[[paste0("period_", paste(period, collapse = "_"))]] <- slope_df
      
      
      
    } else {
      stop("Invalid period specified. Period should be 'full_length' or a numeric vector of length 2.")
    }
    

  }
  
  # Combine each estimate into the same tidy data frame
  senslope_res_list <- imap_dfr(senslope_res_list, ~ tibble(period = .y, .x))
  
  estimate_results <- senslope_res_list %>% 
    mutate(CI_width = uci - lci)
  
  
  
  period_plot <- plot_period_comparison(comp_df, period_df =  estimate_results)
  
  return( list(estimate_results,period_plot ))
}


#undebug(rolling_trend)
#test_rolling <- rolling_trend(test_stl, periods = list("full_length", c(15,0), c(15,10), c(10,5), c(5,0)))

testdf <- rolling_trend(test_stl, periods = list("full_length"))

testdf[[1]]$data

#plot_period_comparison(full_data = test_stl$orig_data[[1]],period_df = test_rolling[[1]])

                     