
#load LandwaterPeople (LWP) functions for trend estimation 
source("For_Andrew/LWPTrends_v2102.R")

library(ggplot2)

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
