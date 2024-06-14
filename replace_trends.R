

source(here::here("simulate_trends_functions.R"))




#source("For_Andrew/simulate_trends_functions.R")
#load LandwaterPeople (LWP) functions for trend estimation 

source(here::here("LWPTrends_v2102.R"))
#source("For_Andrew/LWPTrends_v2102.R")

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
#autoplot(components(test_stl))



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
analyze_trend_with_noise <- function(data, lambda = seq(from = -1, to = 3, by = 0.5), 
                                     is_seasonal = TRUE, 
                                     trend_params = list(modify_trend = NULL,initial_amplitude = NULL), analysis_params = list(), ...){
  
  # Step 1: Get STL decomposition
  stl_data <- Get_STL(data)
  
  
  # join with metadata
  comp_df <- stl_data$stl[[1]]$fit$decomposition
  
  
  
  # Step 2: Optionally modify the trend component
  if (is.null(trend_params$modify_trend)) {
    message("keeping original trend component")
    estimate_results <- scale_stl_noise(stl_data, lambda = lambda, is_seasonal = is_seasonal, analysis_params = analysis_params)
  } else {
    trend_params$total_length <- nrow(comp_df)  
    
    if (trend_params$modify_trend == "cosine") {
      message("replacing original trend component with simulated cosine series")
      mod_fun <- generate_cosine_series
      
      if (is.null(trend_params$initial_amplitude)) {
        trend_params$initial_amplitude <- mean(comp_df$trend, na.rm = TRUE)
        message("setting amplitude as mean of original trend component")
      }
    } else if (trend_params$modify_trend == "level shift") {
      mod_fun <- level_shift_with_ramp
      message("replacing original trend component with simulated level-shift series")
    }  else if (trend_params$modify_trend == "linear") {
      mod_fun <- generate_linear_trend
      message("replacing original trend component with simulated linear trend series")
    } else if (trend_params$modify_trend == "Linex") {
      mod_fun <- generate_linex
      message("replacing original trend component with simulated linex trend series")
    }else {
      stop("Invalid modify_trend option")
    }
    
    # Remove the modify_trend component from trend_params (is not needed for  mod_fun)
    trend_params <- trend_params[-which(names(trend_params) == "modify_trend")]
    
    comp_df$trend <- do.call(mod_fun, trend_params)
    
    comp_df$final_series <-   comp_df$trend + comp_df$remainder
    comp_df$season_adjust <-  comp_df$final_series - comp_df$season_year
    
    
    stl_data$stl[[1]]$fit$decomposition <- comp_df
    
    estimate_results <- scale_stl_noise(stl_data, lambda = lambda, is_seasonal = is_seasonal, analysis_params = analysis_params)
  }
  
  # Step 4: Return the results
  return(list(stl_data = stl_data, estimate_results = estimate_results))
}


# Examples usage of the wrapper function

#linear example
linear_trend_params <- list( modify_trend = "linear",
                             slope = .02, 
                             intercept = 1)

linear_result <- analyze_trend_with_noise(N03N_filtered$`GW-00004`, 
                                          lambda = seq(-1, 3, 0.5), 
                                          is_seasonal = TRUE, 
                                          trend_params = linear_trend_params)


autoplot(components(linear_result$stl_data))

linear_result$estimate_results %>% 
  ggplot(aes(x = as.numeric(lambda), y = est_slope))+
  geom_point()+
  geom_linerange(aes(ymin = lci, ymax = uci))

# Ramp Example
Level_shift_params <- list( modify_trend = "level shift",
                            baseline_amp = .2, 
                            amp_change = .3, 
                            steepness = .5 )

ramp_result <- analyze_trend_with_noise(N03N_filtered$`GW-00004`, 
                                        lambda = seq(-1, 3, 0.5), 
                                        is_seasonal = TRUE, 
                                        trend_params = Level_shift_params)


autoplot(components(ramp_result$stl_data))

ramp_result$estimate_results %>% 
  ggplot(aes(x = as.numeric(lambda), y = est_slope))+
  geom_point()+
  geom_linerange(aes(ymin = lci, ymax = uci))

# Cosine Example
cosine_params <-  list(modify_trend = "cosine",
                       decay_rate = 0.01,      
                       num_peaks = 3,           
                       phase_shift = 0          
)


cos_result <- analyze_trend_with_noise(N03N_filtered$`GW-00002`, 
                                       lambda = seq(-1, 3, 0.5), 
                                       is_seasonal = TRUE, 
                                       trend_params = cosine_params)


autoplot(components(cos_result$stl_data))


cos_result$estimate_results %>% 
  ggplot(aes(x = as.numeric(lambda), y = est_slope))+
  geom_point()+
  geom_linerange(aes(ymin = lci, ymax = uci))



# plot original data from results
cos_result$stl_data$orig_data[[1]] %>% #TODO missing Y value
  ggplot(aes(x = yearmon, y= RawValue))+
  geom_line()+
  geom_point(aes(y = 0, x = yearmon, colour = Censored, alpha = Censored, size = 2))+
  #scale_color_manual(values = c("FALSE" = "white", "TRUE" = "blue"))+
  scale_alpha_manual(values = c("FALSE" = 0, "TRUE" = 1))



ramp_result$stl_data$orig_data[[1]] %>% #TODO missing Y value
  ggplot(aes(x = yearmon, y= RawValue))+
  geom_line()+
  geom_point(aes(y = 0, x = yearmon, colour = Censored, alpha = Censored, size = 2))+
  #scale_color_manual(values = c("FALSE" = "white", "TRUE" = "blue"))+
  scale_alpha_manual(values = c("FALSE" = 0, "TRUE" = 1))



#Linex example
linex_trend_params <- list( modify_trend = "Linex",
                             amplitude = 2, 
                             x_min = 20,
                             scale = 0.5)




linex_result <- analyze_trend_with_noise(N03N_filtered$`GW-00004`, 
                                          lambda = seq(-1, 3, 0.5), 
                                          is_seasonal = TRUE, 
                                          trend_params = linex_trend_params)


autoplot(components(linex_result$stl_data))

linex_result$estimate_results %>% 
  ggplot(aes(x = as.numeric(lambda), y = est_slope))+
  geom_point()+
  geom_linerange(aes(ymin = lci, ymax = uci))
##################This can be added into wrapper also, too much??
# 
# 
# get_adj_slope_trends <- function(data, is_seasonal = NULL, 
#                                  slopes = seq(-0.1, 0.1, by = 0.025)) {
#   #container list
#   estimate_results <- list()
#   simulated_data <- data
#   
#   #test seasonality if not aleady indicated
#   if (is.null(is_seasonal)){
#     seasonal_test <- SeasonalityTest(as.data.frame(data))
#     # indicator if seasonal or not
#     is_seasonal <-  seasonal_test$pvalue <= .05
#   }
#   
#   
#   if ("rowid" %in% colnames(simulated_data) == FALSE){
#     simulated_data <- simulated_data %>% rowid_to_column()
#   }
#   
#   simulated_data <- simulated_data %>% 
#     mutate(Original = RawValue)
#   
#   for (slope in slopes) {
#     # Modify the RawValue based on the slope
#     simulated_data <- simulated_data %>% mutate(RawValue = Original + (rowid * slope))
#     
#     
#     if (is_seasonal == TRUE){
#       senslope_res <- SeasonalTrendAnalysis(as.data.frame(simulated_data)) # doPlot = TRUE 
#     } else {
#       senslope_res <-  NonSeasonalTrendAnalysis(as.data.frame(simulated_data))    
#       
#     }
#     
#     slope_estimate_df <-  tibble(slope_estimate = senslope_res$AnnualSenSlope / 12,
#                                  lci = senslope_res$Sen_Lci/12,
#                                  uci =  senslope_res$Sen_Uci/12,
#                                  CI_width =uci - lci)
#     
#     # Store results
#     estimate_results[[as.character(slope)]] <-  slope_estimate_df 
#     
#   }
#   #tidy output into dataframe
#   
#   estimate_results <- imap_dfr(estimate_results,~ tibble(slope_added = as.numeric(.y), .x, ))
#   
#   
#   estimate_results <- estimate_results %>% 
#     mutate(slope_change = slope_added - lag(slope_added),
#            est_change = slope_estimate - lag(slope_estimate))
#   
#   return(estimate_results)
# }
# 
# 
# ########
# 
# #TODO potenitally modify funciton so defaults are scaled relative to SD of noise component
# 
# analyze_trend_with_noise <- function(data, lambda = seq(from = -1, to = 3, by = 0.5), 
#                                      is_seasonal = TRUE, 
#                                      trend_params = list(modify_trend = NULL,initial_amplitude = NULL), analysis_params = list(), ...){
#   
#   # Step 1: Get STL decomposition
#   stl_data <- Get_STL(data)
#   
#   
#   # join with metadata
#   comp_df <- stl_data$stl[[1]]$fit$decomposition
#   
#   
#   
#   # Step 2: Optionally modify the trend component
#   if (is.null(trend_params$modify_trend)) {
#     message("keeping original trend component")
#     estimate_results <- scale_stl_noise(stl_data, lambda = lambda, is_seasonal = is_seasonal, analysis_params = analysis_params)
#   } else {
#     trend_params$total_length <- nrow(comp_df)  
#     
#     if (trend_params$modify_trend == "cosine") {
#       message("replacing original trend component with simulated cosine series")
#       mod_fun <- generate_cosine_series
#       
#       if (is.null(trend_params$initial_amplitude)) {
#         trend_params$initial_amplitude <- mean(comp_df$trend, na.rm = TRUE)* scale_factor * noise_scale
#         message("setting amplitude as mean of original trend component")
#       }else {
#         trend_params$initial_amplitude <- trend_params$initial_amplitude * scale_factor * noise_scale
#       }
#       
#       
#     } else if (trend_params$modify_trend == "level shift") {
#       mod_fun <- level_shift_with_ramp
#       message("replacing original trend component with simulated level-shift series")
#     }  else if (trend_params$modify_trend == "linear") {
#       mod_fun <- generate_linear_trend
#       message("replacing original trend component with simulated linear trend series")
#     } else if (trend_params$modify_trend == "Linex") {
#       mod_fun <- generate_linex
#       message("replacing original trend component with simulated linex trend series")
#     }else {
#       stop("Invalid modify_trend option")
#     }
#     
#     # Remove the modify_trend component from trend_params (is not needed for  mod_fun)
#     trend_params <- trend_params[-which(names(trend_params) == "modify_trend")]
#     
#     comp_df$trend <- do.call(mod_fun, trend_params)
#     
#     comp_df$final_series <-   comp_df$trend + comp_df$remainder
#     comp_df$season_adjust <-  comp_df$final_series - comp_df$season_year
#     
    
    stl_data$stl[[1]]$fit$decomposition <- comp_df
    
    estimate_results <- scale_stl_noise(stl_data, lambda = lambda, is_seasonal = is_seasonal, analysis_params = analysis_params)
  }
  
  # Step 4: Return the results
  return(list(stl_data = stl_data, estimate_results = estimate_results))
}
    
    # Remove the modify_trend component from trend_params (is not needed for mod_fun)
    trend_params <- trend_params[-which(names(trend_params) == "modify_trend")]
    
    comp_df$trend <- do.call(mod_fun, trend_params)
    
    comp_df$final_series <- comp_df$trend + comp_df$remainder
    comp_df$season_adjust <- comp_df$final_series - comp_df$season_year
    
    stl_data$stl[[1]]$fit$decomposition <- comp_df
    
    estimate_results <- scale_stl_noise(stl_data, lambda = lambda, is_seasonal = is_seasonal, analysis_params = analysis_params)
  }
  
  # Step 4: Return the results
  return(list(stl_data = stl_data, estimate_results = estimate_results))
}

# Ramp Example
Level_shift_params <- list( modify_trend = "level shift",
                            baseline_amp = .1, 
                            amp_change = 100, 
                            steepness = .8 )

ramp_result <- analyze_trend_with_noise(N03N_filtered$`GW-00004`, 
                                        lambda = seq(-1, 3, 0.5), 
                                        is_seasonal = TRUE, 
                                        trend_params = Level_shift_params,scale_factor = 1)


# analyze_trend_with_noise(data, 
#                          lambda = seq(-1, 3, 0.5), 
#                          is_seasonal = TRUE, 
#                          trend_params = list(modify_trend = "cosine", decay_rate = 0.1, num_peaks = 5, phase_shift = 0), 
#                          analysis_params = list(doPlot = TRUE, some_other_param = "example"),
#                          scale_factor = 1.2)



autoplot(components(ramp_result$stl_data))

ramp_result$estimate_results %>% 
  ggplot(aes(x = as.numeric(lambda), y = est_slope))+
  geom_point()+
  geom_linerange(aes(ymin = lci, ymax = uci))
