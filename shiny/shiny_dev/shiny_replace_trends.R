
#load LandwaterPeople (LWP) functions for trend estimation 

#source(here::here("LWPTrends_v2102.R"))

#library(tidyverse)
# library(ggplot2)
# library(withr)
# library(feasts)
# library(tsibble)

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
      fabletools::model(stl = STL(final_series))
    
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
#' @param trend_params A list of parameters for modifying the trend component, including `initial_amplitude` (see simulate_trends_functions.R) .
#' @param analysis_params A list of additional parameters for trend analysis functions.
#' @param ... Additional arguments passed to the trend modification function.
#' @return A list containing the STL decomposition data and the estimated results.
#' @export
analyze_trend_with_noise <- function(data, lambda = seq(from = -1, to  = 3, by =  0.5), 
                                     is_seasonal = TRUE, 
                                     trend_params = list(initial_amplitude = NULL), 
                                     analysis_params = list(), mod_fun = NULL, ...){
  # Step 1: Get STL decomposition
  stl_data <- Get_STL(data)
  
  # Join with metadata
  comp_df <- stl_data$stl[[1]]$fit$decomposition
  
  # Step 2: Optionally modify the trend component

  if (is.null(mod_fun)) {
    message("Keeping original trend component")
    estimate_results <- scale_stl_noise(stl_data, lambda = lambda, is_seasonal = is_seasonal, analysis_params = analysis_params)

  } else {
    trend_params$total_length <- nrow(comp_df)  
    message("simulating trend using  function provided from 'mod_fun' parameter")
  
    
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




# Set cosine_params to match  generate_cosine_series (total_length, initial_amplitude, decay_rate, num_peaks, phase_shift)

# cosine_params <-  list(
#                        decay_rate = 0.01,
#                        initial_amplitude = 5,
#                        num_peaks = 5,
#                        phase_shift = 0
# )
# # 
# # 
# test_output_lamba <- analyze_trend_with_noise(N03N_filtered$`GW-00002`,
#                                               lambda = seq(-1, 3, 0.5),
#                                               is_seasonal = TRUE,
#                                               trend_params = cosine_params,
#                                               mod_fun = generate_cosine_series)

# analyze_trend_with_noise(N03N_filtered$`GW-00002`,
#                                               lambda = seq(-1, 3, 0.5),
#                                               is_seasonal = TRUE,
#                                               trend_params = cosine_params,
#                                               mod_fun = testNULL)
# 
# testNULL <- as.null()

#
# # 
#  autoplot(components(test_output_lamba$stl_data)) 
# 

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
# plot_period_comparison <- function(full_data, period_df){
#   
#   orig_plot <- full_data %>% 
#     ggplot(aes(x = yearmon)) +
#     geom_line(aes(y = final_series))+
#     theme_bw()
#   
#   
#   slope_df <- period_df %>% 
#     tidyr::hoist(data, "yearmon") %>% 
#     tidyr::unnest(yearmon)
#   
#   
#   slope_df <- slope_df %>% 
#     mutate(period = as.factor(period)) %>% 
#     dplyr::group_by(period ) %>% 
#     dplyr::mutate(rowid = dplyr::row_number(),
#            slope_line = est_slope*rowid,
#            slope_line_uci = uci*rowid,
#            slope_line_lci = lci*rowid,
#     )
#   
#   # Adding the slope line
#   orig_plot <- orig_plot +
#     geom_line(data = slope_df, aes(x = yearmon, y = slope_line, group = period), 
#               show.legend = FALSE,
#               size = 1.5,
#               colour = "grey")
#   
#   # Adding the slope CI 
#   orig_plot <- orig_plot +
#     geom_ribbon(data = slope_df, 
#                 aes(x = yearmon, ymin = slope_line_lci, ymax = slope_line_uci, group = period), 
#                 alpha = 0.15, color = NA, show.legend = FALSE, fill = "grey")
#   
#   return(orig_plot)
#   
# }

plot_period_comparison <- function(full_data, period_df){
  
  
  
  ##
  MK_data <- period_df$MK
  names(MK_data) <- period_df$period
  MK_data <- imap_dfr(MK_data, ~tibble(period = .y,.x))
  #add conf_int
  MK_data <- get_ConfCat_from_MK_data(MK_data)
  
  period_df <- left_join(period_df, select(MK_data,ConfCat, period ))
  
  
  
  orig_plot <- full_data %>% 
    ggplot(aes(x = yearmon)) +
    geom_line(aes(y = final_series))+
    theme_bw()
  
  
  slope_df <- period_df %>% 
    tidyr::hoist(data, "yearmon") %>% 
    tidyr::unnest(yearmon)
  
  
  slope_df <- slope_df %>% 
    mutate(period = as.factor(period)) %>% 
    dplyr::group_by(period ) %>% 
    dplyr::mutate(rowid = dplyr::row_number(),
                  slope_line = est_slope*rowid,
                  slope_line_uci = uci*rowid,
                  slope_line_lci = lci*rowid,
    )
  
  # Adding the slope line
  orig_plot <- orig_plot +
    geom_line(data = slope_df, aes(x = yearmon, y = slope_line, group = period,
                                   colour = ConfCat), 
              show.legend = FALSE,
              size = 1)
  
  
  
  
  # Adding the slope CI 
  orig_plot <- orig_plot +
    geom_ribbon(data = slope_df, 
                aes(x = yearmon, ymin = slope_line_lci, 
                    ymax = slope_line_uci,
                    group = period,
                    fill = ConfCat), 
                alpha = 0.1, show.legend = FALSE, )+
    scale_fill_manual(values = color_mapping, drop = FALSE)+
    scale_color_manual(values = color_mapping, drop = FALSE)
  
  
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
           yearmon, Season, Year, myDate)
  
  # Add columns so LWP functions work
  comp_df <- comp_df %>% left_join(orig_data, by = "yearmon")
  #needs this name for LWP functions
  comp_df$RawValue <- comp_df$final_series
  
  # Iterate over each period and perform trend analysis
  for (period in periods) {
    
    if (is.character(period) && period == "full_length") {
      
      
      #get date range in string
      date_range <- range(comp_df$yearmon)
      date_range <- paste(date_range[1],   date_range[2], sep = ":")
      
      # Full length analysis without filtering
      if (is_seasonal) {
        senslope_res <- do.call(SeasonalTrendAnalysis, c(list(as.data.frame(comp_df)), analysis_params))
      } else {
        senslope_res <- do.call(NonSeasonalTrendAnalysis, c(list(as.data.frame(comp_df)), analysis_params))
      } 
        
        # Create a data frame with the results
        slope_df <- tibble(date_range = date_range,
                           est_slope = senslope_res$AnnualSenSlope / 12,
                           lci = senslope_res$Sen_Lci / 12,
                           uci = senslope_res$Sen_Uci / 12,
                           CI_width = uci - lci,
                           data = list(comp_df),
                           MK = list(senslope_res))
        
        senslope_res_list[[paste0("period_", paste(period, collapse = "_"))]] <- slope_df  
        
        

    } else if (is.numeric(period) && length(period) == 2) {
      # Filter the data for the specified period
      comp_df_filt <- filter_custom_period(comp_df, period)
      #get date range in string
     date_range <- range(comp_df_filt$yearmon)
      date_range <- paste(date_range[1],   date_range[2], sep = ":")
      
      
      if (is_seasonal) {
        senslope_res <- do.call(SeasonalTrendAnalysis, c(list(as.data.frame(comp_df_filt)), analysis_params))
      } else {
        senslope_res <- do.call(NonSeasonalTrendAnalysis, c(list(as.data.frame(comp_df_filt)), analysis_params))
      }
      
      # Create a data frame with the results
      slope_df <- tibble(date_range = date_range,
                         est_slope = senslope_res$AnnualSenSlope / 12,
                         lci = senslope_res$Sen_Lci / 12,
                         uci = senslope_res$Sen_Uci / 12,
                         CI_width = uci - lci,
                         data = list(comp_df_filt),
                         MK = list(senslope_res))
      
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



#test_rolling <- rolling_trend(test_stl, periods = list("full_length"))

#debug(rolling_trend)
############     


analyze_trend_rolling <- function(data, 
                                  is_seasonal = TRUE, 
                                  trend_params = list(initial_amplitude = NULL), 
                                  analysis_params = list(), mod_fun = NULL, periods = list("full_length"), ...){
  # Step 1: Get STL decomposition
  stl_data <- Get_STL(data)
  
  # Join with metadata
  comp_df <- stl_data$stl[[1]]$fit$decomposition
  
  # Step 2: Optionally modify the trend component
  if (is.null(mod_fun)) {
    message("Keeping original trend component")
    estimate_results <- rolling_trend(stl_data, periods = periods, is_seasonal = is_seasonal, analysis_params = analysis_params) 
    
  } else {
    trend_params$total_length <- nrow(comp_df)  
    message("simulating trend using  function provided from 'mod_fun' parameter")
    
    
    # Modify trend component
    comp_df$trend <- do.call(mod_fun, trend_params)
    
    # Update STL data
    comp_df$final_series <- comp_df$trend + comp_df$remainder
    comp_df$season_adjust <- comp_df$final_series - comp_df$season_year
    
    ###rolling window for   analysis
    

    stl_data$stl[[1]]$fit$decomposition <- comp_df
    
    #estimate_results <- scale_stl_noise(stl_data, lambda = lambda, is_seasonal = is_seasonal, analysis_params = analysis_params)
     estimate_results <- rolling_trend(stl_data, periods = periods, is_seasonal = is_seasonal, analysis_params = analysis_params) 
    
    
  }
  
  # Step 4: Return the results
  return(list(stl_data = stl_data, estimate_results = estimate_results))
}



##
get_GAM_trend <- function(stl_data,  GAM_params = list()){
  
  
  comp_df <- components(stl_data) %>% as_tibble()
  
  # Needs extra columns for LWP functions to work
  orig_data <- stl_data$orig_data[[1]] %>% 
    select(lawa_site_id, CenType, Censored, 
           yearmon, Season, Year, myDate)
  
  # Add columns so LWP functions work
  comp_df <- comp_df %>% left_join(orig_data, by = "yearmon")
  #needs this name for LWP functions
  comp_df$RawValue <- comp_df$final_series
  
  
  GAM_results <-  screeningmodeling(.data = comp_df , datevar = myDate,values = RawValue )
  
  # Create a data frame with the results
  GAM_plot <- plot_individual_trend(GAM_results)
  
  return(list(GAM_results,GAM_plot) )
}




analyze_GAM_wrapper <- function(data, 
                                trend_params = list(initial_amplitude = NULL), 
                                GAM_params= list(),
                                mod_fun = NULL, ...) {
  # Step 1: Get STL decomposition
  stl_data <- Get_STL(data)
  
  # Join with metadata
  comp_df <- stl_data$stl[[1]]$fit$decomposition
  
  # Step 2: Optionally modify the trend component
  if (is.null(mod_fun)) {
    message("Keeping original trend component")
    GAM_results <- get_GAM_trend(stl_data = stl_data, GAM_params = GAM_params)
    
  } else {
    trend_params$total_length <- nrow(comp_df)  
    message("simulating trend using  function provided from 'mod_fun' parameter")
    
    
    # Modify trend component
    comp_df$trend <- do.call(mod_fun, trend_params)
    
    # Update STL data
    comp_df$final_series <- comp_df$trend + comp_df$remainder
    comp_df$season_adjust <- comp_df$final_series - comp_df$season_year
    
    
    stl_data$stl[[1]]$fit$decomposition <- comp_df
    
    
    GAM_results <- get_GAM_trend(stl_data = stl_data, GAM_params= GAM_params )
  }
  
  
  
  # Step 4: Return the results
  return(GAM_results)
}


# cosine_params <-  list(
#                        decay_rate = 0.01,
#                        initial_amplitude = 5,
#                        num_peaks = 5,
#                        phase_shift = 0
# )

# test_rol_td <- analyze_trend_rolling(N03N_filtered$`GW-00002`,
#                                periods = list("full_length", c(10,0), c(8,3), c(6,1), c(5,0)),
#                                is_seasonal = TRUE,
#                                trend_params = cosine_params,
#                                mod_fun = generate_cosine_series)
# 
# 
# debug(analyze_trend_rolling)


# test_rol_td <- analyze_trend_rolling(N03N_filtered$`GW-00002`,
#                                periods = list("full_length", c(10,0), c(8,3), c(6,1), c(5,0)),
#                                is_seasonal = TRUE,
#                                trend_params = cosine_params,
#                                mod_fun = NULL)

# test_rol_td <- analyze_trend_rolling(N03N_filtered$`GW-00002`,
#                                      #periods = list("full_length", c(10,0), c(8,3), c(6,1), c(5,0)),
#                                      is_seasonal = TRUE,
#                                      trend_params = cosine_params,
#                                      mod_fun = NULL)

###




color_mapping <- c(
  "Very likely improving" = "green4",
  "Likely improving" = "seagreen2",
  "Indeterminate" = "blue",
  "Likely degrading" = "orange",
  "Very likely degrading" = "red"
)
# from cawthron  
# add confidence classes
###modify so take input direct from rolling results
get_ConfCat <-  function(rolling_results){
  
  #get MK data stored in results from 'analyze_trend_rolling'
  MK_data <- rolling_results$estimate_results[[1]]$MK
  names( MK_data) <- rolling_results$estimate_results[[1]]$date_range 
  MK_data <- imap_dfr(MK_data, ~tibble(period = .y,.x))
  
  conf_cat_df <- MK_data %>%
    mutate(
      ConfCat = cut(Cd,
                    breaks = c(-0.1, 0.1, 0.33, 0.67, 0.90, 1.1),
                    labels = c(
                      "Very likely improving", "Likely improving", "Indeterminate",
                      "Likely degrading", "Very likely degrading"
                    )
      ),
      ConfCat = factor(ConfCat,
                       levels = rev(c(
                         "Very likely improving",
                         "Likely improving",
                         "Indeterminate",
                         "Likely degrading",
                         "Very likely degrading"
                       ))
      ),
      TrendScore = as.numeric(ConfCat) - 3,
      TrendScore = ifelse(is.na(TrendScore), NA, TrendScore))
  
  return(conf_cat_df )
}
## 

# plot MK estimates
get_MK_plot <- function(MKdata, equiv_zone = c(-.015, .015)){
  
  MK_plot <- MKdata  %>% 
    ggplot(aes(colour = ConfCat))+
    geom_point(aes(y = AnnualSenSlope, x = period,size = 1),show.legend = F)+
    geom_linerange(aes(ymin = Sen_Lci, ymax = Sen_Uci, x = period), size = 1.5)+
    geom_hline(yintercept = 0, linetype = "dashed", size = .3)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = -90, hjust = 1))+
    # Rectangle for the equivalence zone
    geom_rect(aes(ymin = min(equiv_zone), ymax = max(equiv_zone), xmin = -Inf, xmax = Inf),
              fill = "grey80", alpha = 0.12, colour = NA) +
    scale_color_manual(values = color_mapping, drop = FALSE)
  
  return( MK_plot)   
}



# scale components
scale_components<- function(stl_data, scaling_factor = list(c(1,1)), 
                            is_seasonal = TRUE, analysis_params = list()) {
  
  # Container lists
  complist <- list()
  senslope_res_list <- list()
  
  #comp_df <- components(stl_data) #%>% as_tibble()
  comp_df <- stl_data$stl[[1]]$fit$decomposition
  # Needs extra columns for LWP functions to work
  orig_data <- stl_data$orig_data[[1]] %>% 
    select(lawa_site_id, CenType, Censored, 
           yearmon, Season, Year, myDate, RawValue)
  
  # Add columns so LWP functions work
  comp_df <- comp_df %>% left_join(orig_data, by = "yearmon")
  
  ## Scale noise by lambda here in a for loop (estimating slope each time)
  for (i in 1:length(scaling_factor)) {
    
    # comp_df <- comp_df %>% 
    #    mutate(final_series = trend +  (season_year * scaling_factor[[i]][1]) + (remainder * scaling_factor[[i]][2]),
    #           RawValue = final_series)  # Only so named because LWP functions expect that name
    
    
    # Update STL data
    
    comp_df$remainder <- comp_df$remainder* scaling_factor[[i]][2]
    comp_df$season_year <-comp_df$season_year*scaling_factor[[i]][1]
    
    comp_df$final_series <- comp_df$trend + comp_df$season_year + comp_df$remainder
    comp_df$season_adjust <- comp_df$trend - comp_df$season_year
    
    #need this value name for LWP function
    comp_df$RawValue <- comp_df$final_series     
    
    if (is_seasonal == TRUE) {
      senslope_res <- do.call(SeasonalTrendAnalysis, c(list(as.data.frame(comp_df)), analysis_params))  # Pass additional arguments
    } else {
      senslope_res <- do.call(NonSeasonalTrendAnalysis, c(list(as.data.frame(comp_df)), analysis_params))  # Pass additional arguments
    }
    
    #replace STL components
    stl_data_mod <- stl_data
    stl_data_mod$stl[[1]]$fit$decomposition <- comp_df
    
    #components(stl_data_mod) %>% autoplot()
    
    slope_df <- tibble(est_slope = senslope_res$AnnualSenSlope / 12,
                       lci = senslope_res$Sen_Lci / 12,
                       uci = senslope_res$Sen_Uci / 12,
                       CI_width = uci - lci,
                       seasonal_sf = scaling_factor[[i]][1],
                       noise_sf = scaling_factor[[i]][2],
                       data = list(stl_data_mod),
                       seasonal = is_seasonal,
                       MKdata_res = senslope_res)
    
    senslope_res_list[[i]] <- slope_df
  }
  
  # Combine each estimate into same tidy dataframe
  senslope_res_list <- map_df(senslope_res_list, ~ .x)
  
  estimate_results <- senslope_res_list %>% 
    mutate(CI_width = uci - lci)
  
  return(estimate_results)
}



# new wrapper that allows scaling and setting periods ---------------------



rolling_trend_alt <- function(stl_data, periods = list("full_length", c(5, 0)), 
                              is_seasonal = TRUE, analysis_params = list(),
                              ...) {
  
  # Container lists
  complist <- list()
  senslope_res_list <- list()
  
  comp_df <- components(stl_data) %>% as_tibble()
  

  # Iterate over each period and perform trend analysis
  for (period in periods) {
    
    if (is.character(period) && period == "full_length") {
      
      
      #get date range in string
      date_range <- range(comp_df$yearmon)
      date_range <- paste(date_range[1],   date_range[2], sep = ":")
      
      # Full length analysis without filtering
      if (is_seasonal) {
        senslope_res <- do.call(SeasonalTrendAnalysis, c(list(as.data.frame(comp_df)), analysis_params))
      } else {
        senslope_res <- do.call(NonSeasonalTrendAnalysis, c(list(as.data.frame(comp_df)), analysis_params))
      } 
      
      # Create a data frame with the results
      slope_df <- tibble(date_range = date_range,
                         est_slope = senslope_res$AnnualSenSlope / 12,
                         lci = senslope_res$Sen_Lci / 12,
                         uci = senslope_res$Sen_Uci / 12,
                         CI_width = uci - lci,
                         data = list(comp_df),
                         MK = list(senslope_res))
      
      senslope_res_list[[paste0("period_", paste(period, collapse = "_"))]] <- slope_df  
      
      
      
    } else if (is.numeric(period) && length(period) == 2) {
      # Filter the data for the specified period
      comp_df_filt <- filter_custom_period(comp_df, period)
      #get date range in string
      date_range <- range(comp_df_filt$yearmon)
      date_range <- paste(date_range[1],   date_range[2], sep = ":")
      
      
      if (is_seasonal) {
        senslope_res <- do.call(SeasonalTrendAnalysis, c(list(as.data.frame(comp_df_filt)), analysis_params))
      } else {
        senslope_res <- do.call(NonSeasonalTrendAnalysis, c(list(as.data.frame(comp_df_filt)), analysis_params))
      }
      
      # Create a data frame with the results
      slope_df <- tibble(date_range = date_range,
                         est_slope = senslope_res$AnnualSenSlope / 12,
                         lci = senslope_res$Sen_Lci / 12,
                         uci = senslope_res$Sen_Uci / 12,
                         CI_width = uci - lci,
                         data = list(comp_df_filt),
                         MK = list(senslope_res))
      
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

# scale components
scale_components<- function(stl_data, 
                            scaling_factor = list(c(1,1)), 
                            periods = list("full_length", c(5, 0)), 
                            is_seasonal = TRUE, 
                            analysis_params = list(),
                            ...) {
  
  # Container lists
  complist <- list()
  senslope_res_list <- list()
  
  #comp_df <- components(stl_data) #%>% as_tibble()
  comp_df <- stl_data$stl[[1]]$fit$decomposition
  # Needs extra columns for LWP functions to work
  orig_data <- stl_data$orig_data[[1]] %>% 
    select(lawa_site_id, CenType, Censored, 
           yearmon, Season, Year, myDate, RawValue)
  
  # Add columns so LWP functions work
  comp_df <- comp_df %>% left_join(orig_data, by = "yearmon")
  
  ## Scale noise by lambda here in a for loop (estimating slope each time)
  for (i in 1:length(scaling_factor)) {
    
    #set to original
    comp_df_mod <- comp_df
    # Update STL data
    
    comp_df_mod$remainder <- comp_df_mod$remainder* scaling_factor[[i]][2]
    comp_df_mod$season_year <-comp_df_mod$season_year*scaling_factor[[i]][1]
    
    comp_df_mod$final_series <- comp_df_mod$trend + comp_df_mod$season_year + comp_df_mod$remainder
    comp_df_mod$season_adjust <- comp_df_mod$trend - comp_df_mod$season_year
    
    #need this value name for LWP function
    comp_df_mod$RawValue <- comp_df_mod$final_series   
    
    # #replace STL components
    stl_data_mod <- stl_data
    stl_data_mod$stl[[1]]$fit$decomposition <- comp_df_mod
    
    rol_est <- rolling_trend_alt(stl_data_mod,
                                 is_seasonal = is_seasonal, 
                                 periods = periods,
                                 analysis_params = analysis_params)
    
    rol_est[[1]]$seasonal_sf <- scaling_factor[[i]][1]
    rol_est[[1]]$noise_sf = scaling_factor[[i]][2]
    
    rol_est[[1]] <- rol_est[[1]] %>% relocate(period,seasonal_sf,noise_sf)
    
    senslope_res_list[[i]] <- rol_est
  }
  return(senslope_res_list)
}  


analyze_trend_wrapper <- function(data, 
                                      scaling_factor = list(c(1,1)),
                                      is_seasonal = TRUE, 
                                      trend_params = list(initial_amplitude = NULL), 
                                      analysis_params = list(), mod_fun = NULL, 
                                      periods = list("full_length"), 
                                      ...){
  # Step 1: Get STL decomposition
  stl_data <- Get_STL(data)
  
  # Join with metadata
  comp_df <- stl_data$stl[[1]]$fit$decomposition
  
  # Step 2: Optionally modify the trend component
  if (is.null(mod_fun)) {
    message("Keeping original trend component")
    #estimate_results <- rolling_trend_alt(stl_data, periods = periods, is_seasonal = is_seasonal, analysis_params = analysis_params) 
    estimate_results <- scale_components(stl_data,  periods = periods, 
                                         scaling_factor = scaling_factor,
                                         is_seasonal = is_seasonal, 
                                         analysis_params = analysis_params)
  } else {
    trend_params$total_length <- nrow(comp_df)  
    message("simulating trend using  function provided from 'mod_fun' parameter")
    
    
    # Modify trend component
    comp_df$trend <- do.call(mod_fun, trend_params)
    
    # Update STL data
    comp_df$final_series <- comp_df$trend + comp_df$remainder
    comp_df$season_adjust <- comp_df$final_series - comp_df$season_year
    
    ###rolling window for   analysis
    
    
    stl_data$stl[[1]]$fit$decomposition <- comp_df
    
 
    estimate_results <- scale_components(stl_data,  
                                         periods = periods, 
                                         scaling_factor = scaling_factor,
                                         is_seasonal = is_seasonal, 
                                         analysis_params = analysis_params)
    
  }
  
  # Step 4: Return the results
  #return(list(stl_data = stl_data, estimate_results = estimate_results))
  return(estimate_results)
}

##Example
# cosine_params <-  list(
#   decay_rate = 0.01,
#   initial_amplitude = 2,
#   num_peaks = 2,
#   phase_shift = 0
# )
# 
# # first element = seasonal, second element = remainder
# # e.g.1, c(1,1): original seasonal and remanider component
# #e.g.,2  c(0,2): remove seasonal, double  remainder component
# scaling_factor <- list(c(1,1), c(0,1), c(2,1), c(2,2))
# 
         # MK_result<-        analyze_trend_wrapper(site_data$`GW-00002`,
         #                                 scaling_factor = scaling_factor,
         #                                 periods = list("full_length", c(10,5)),
         #                                 is_seasonal = TRUE,
         #                                 trend_params = cosine_params,
         #                                 mod_fun = generate_cosine_series)
         # 
       
         #   output$summary1 <- renderPlot({
         #     
             # Extract seasonal_sf for the title
         #     seasonal_sf <- unique(result[[1]][[1]]$seasonal_sf)
         #     noise_sf <- unique(result[[1]][[1]]$noise_sf)
         # 
         #   result[[1]][[2]]+
         #     labs(title = paste("Seasonal SF:", seasonal_sf, ",", " Noise SF: ",noise_sf ))
         # 
#          # 
#            plot_list <- map(MK_result, ~ {seasonal_sf <- unique(.x[[1]]$seasonal_sf)
#                           noise_sf <- unique(.x[[1]]$noise_sf)
# 
#                           plot <- .x[[2]]+
#                           labs(title = paste("Seasonal SF:", seasonal_sf, ",", " Noise SF: ",noise_sf ))
# 
#                           return(plot)
#                           })
#            
#            
#            MK_results_list <- map(MK_result, ~ {seasonal_sf <- unique(.x[[1]]$seasonal_sf)
#            noise_sf <- unique(.x[[1]]$noise_sf)
#            
#            MK <- .x[[1]]
#      
#            
#            return(plot)
#            })  
#            
#            MK_list <- map(MK_result, ~ select(.[[1]], 1:4, MK) %>% unnest(MK) )
#                  
#            MK_df <-  map_df(MK_list, ~ tibble(.))
#         
#            
#            
# MK_data_list <- map(MK_result, ~ select(.[[1]], 1:4, data) %>% unnest(data) )

#MK_data_df <-  map_df(MK_data_list, ~ tibble(.))




         # 
         #   # Dynamically wrap plots
         #   wrapped_plots <- do.call(patchwork::wrap_plots, c(plot_list, ncol = 1))
         # 
         # })

# result

  
get_ConfCat_from_MK_data <-  function(MK_data ){
  
  conf_cat_df <- MK_data %>%
    mutate(
      ConfCat = cut(Cd,
                    breaks = c(-0.1, 0.1, 0.33, 0.67, 0.90, 1.1),
                    labels = c(
                      "Very likely improving", "Likely improving", "Indeterminate",
                      "Likely degrading", "Very likely degrading"
                    )
      ),
      ConfCat = factor(ConfCat,
                       levels = rev(c(
                         "Very likely improving",
                         "Likely improving",
                         "Indeterminate",
                         "Likely degrading",
                         "Very likely degrading"
                       ))
      ),
      TrendScore = as.numeric(ConfCat) - 3,
      TrendScore = ifelse(is.na(TrendScore), NA, TrendScore))
  
  return(conf_cat_df )
}

#QR functions


get_QR_est <- function(data, trm = "Month", 
                       form = formula(paste0("RawValue~Year+", trm)),
                       analysis_params = list()) {
  
  data$RawValue[data$RawValue == 0] <- 0.00001  
  
  data$myDate <-   as.Date(data$yearmon) 
  data$Year <- lubridate::year( data$myDate )
  data$Month <- as.factor(lubridate::month(data$myDate ))
  
  
  
  
  qr_med <-  quantreg::rq(formula = form,
                          data =  data,
                          tau =  .5, #median
                          contrasts = list(Month = contr.sum),
                          method = "fn")
  
  ##
  
  sjplot <- qr_med %>%  sjPlot::plot_model(terms = "Year") 
  conf_data <- sjplot$data
  
  conf_cat_df <- conf_data  %>%
    mutate(Cd = if_else(estimate < 0, 1 - p.value/2, p.value/2),
           ConfCat = cut(Cd,
                         breaks = c(-0.1, 0.1, 0.33, 0.67, 0.90, 1.1),
                         labels = c(
                           "Very likely improving", "Likely improving", "Indeterminate",
                           "Likely degrading", "Very likely degrading"
                         )
           ),
           ConfCat = factor(ConfCat,
                            levels = rev(c(
                              "Very likely improving",
                              "Likely improving",
                              "Indeterminate",
                              "Likely degrading",
                              "Very likely degrading"
                            ))
           ),
           TrendScore = as.numeric(ConfCat) - 3,
           TrendScore = ifelse(is.na(TrendScore), NA, TrendScore))
  
  return(conf_cat_df)
  
}


rolling_trend_QR <- function(stl_data, periods = list("full_length", c(5, 0)), 
                             is_seasonal = TRUE, analysis_params = list(),
                             ...) {
  
  # Container lists
  complist <- list()
  senslope_res_list <- list()
  
  comp_df <- components(stl_data) %>% as_tibble()
  
  #Removed getting of orig data, as already from other function.
  
  #needs this name for LWP functions
  comp_df$RawValue <- comp_df$final_series
  
  
  # Iterate over each period and perform trend analysis
  for (period in periods) {
    
    if (is.character(period) && period == "full_length") {
      
      
      #get date range in string
      date_range <- range(comp_df$yearmon)
      date_range <- paste(date_range[1],   date_range[2], sep = ":")
      
      # Full length analysis without filtering
      senslope_res <- do.call(get_QR_est, c(list(as.data.frame(comp_df)), analysis_params))
      
      
      senslope_res_list[[paste0("period_", paste(period, collapse = "_"))]] <- senslope_res 
      
      
      
    } else if (is.numeric(period) && length(period) == 2) {
      
      # Filter the data for the specified period
      comp_df_filt <- filter_custom_period(comp_df, period)
      #get date range in string
      date_range <- range(comp_df_filt$yearmon)
      date_range <- paste(date_range[1],   date_range[2], sep = ":")
      
      
      senslope_res <- do.call(get_QR_est, c(list(as.data.frame(comp_df_filt)), analysis_params))
      
      
      senslope_res_list[[paste0("period_", paste(period, collapse = "_"))]] <- senslope_res 
      
      
    } else {
      stop("Invalid period specified. Period should be 'full_length' or a numeric vector of length 2.")
    }
    
    
  }
  
  
  return( senslope_res_list)
}
# scale components
scale_components_QR<- function(stl_data, 
                               scaling_factor = list(c(1,1)), 
                               periods = list("full_length", c(5, 0)), 
                               is_seasonal = TRUE, 
                               analysis_params = list(),
                               ...) {
  
  # Container lists
  complist <- list()
  QR_res_list <- list()
  
  #comp_df <- components(stl_data) #%>% as_tibble()
  comp_df <- stl_data$stl[[1]]$fit$decomposition
  # Needs extra columns for LWP functions to work
  orig_data <- stl_data$orig_data[[1]] %>% 
    select(lawa_site_id, CenType, Censored, 
           yearmon, Season, Year, myDate, RawValue)
  
  # Add columns so LWP functions work
  comp_df <- comp_df %>% left_join(orig_data, by = "yearmon")
  
  ## Scale noise by lambda here in a for loop (estimating slope each time)
  for (i in 1:length(scaling_factor)) {
    
    #set to original
    comp_df_mod <- comp_df
    # Update STL data
    
    comp_df_mod$remainder <- comp_df_mod$remainder* scaling_factor[[i]][2]
    comp_df_mod$season_year <-comp_df_mod$season_year*scaling_factor[[i]][1]
    
    comp_df_mod$final_series <- comp_df_mod$trend + comp_df_mod$season_year + comp_df_mod$remainder
    comp_df_mod$season_adjust <- comp_df_mod$trend - comp_df_mod$season_year
    
    #need this value name for LWP function
    comp_df_mod$RawValue <- comp_df_mod$final_series   
    
    # #replace STL components
    stl_data_mod <- stl_data
    stl_data_mod$stl[[1]]$fit$decomposition <- comp_df_mod
    
    
    
    QR_est <- rolling_trend_QR(stl_data_mod,is_seasonal = is_seasonal, 
                               periods = periods,
                               analysis_params = analysis_params)
    
    # rol_est <- rolling_trend_alt(stl_data_mod,
    #                              is_seasonal = is_seasonal, 
    #                              periods = periods,
    #                              analysis_params = analysis_params)
    
    
    
    QR_est <- map(QR_est, ~ {.$seasonal_sf <- scaling_factor[[i]][1]
    .$noise_sf <-  scaling_factor[[i]][2] 
    return(.) } )            
    
    
    
    QR_res_list[[i]] <- QR_est
  }
  return(QR_res_list)
}  


###
analyze_trend_wrapper_QR <- function(data, 
                                     scaling_factor = list(c(1,1)),
                                     is_seasonal = TRUE, 
                                     trend_params = list(initial_amplitude = NULL), 
                                     analysis_params = list(), mod_fun = NULL, 
                                     periods = list("full_length"), 
                                     ...){
  # Step 1: Get STL decomposition
  stl_data <- Get_STL(data)
  
  # Join with metadata
  comp_df <- stl_data$stl[[1]]$fit$decomposition
  
  # Step 2: Optionally modify the trend component
  if (is.null(mod_fun)) {
    message("Keeping original trend component")
    #estimate_results <- rolling_trend_alt(stl_data, periods = periods, is_seasonal = is_seasonal, analysis_params = analysis_params) 
    estimate_results <- scale_components_QR(stl_data,  periods = periods, 
                                            scaling_factor = scaling_factor,
                                            is_seasonal = is_seasonal, 
                                            analysis_params = analysis_params)
  } else {
    trend_params$total_length <- nrow(comp_df)  
    message("simulating trend using  function provided from 'mod_fun' parameter")
    
    
    # Modify trend component
    comp_df$trend <- do.call(mod_fun, trend_params)
    
    # Update STL data
    comp_df$final_series <- comp_df$trend + comp_df$remainder
    comp_df$season_adjust <- comp_df$final_series - comp_df$season_year
    
    ###rolling window for   analysis
    
    
    stl_data$stl[[1]]$fit$decomposition <- comp_df
    
 
    estimate_results <- scale_components_QR(stl_data,  
                                            periods = periods, 
                                            scaling_factor = scaling_factor,
                                            is_seasonal = is_seasonal, 
                                            analysis_params = analysis_params)
    
  }
  
  # Step 4: Return the results
  #return(list(stl_data = stl_data, estimate_results = estimate_results))
  return(estimate_results)
}



# result_QR <- analyze_trend_wrapper_QR(site_data$`GW-00002`,
#                                       scaling_factor = scaling_factor,
#                                       periods = list("full_length", c(10,0),c(7,0)),
#                                       is_seasonal = TRUE,
#                                       trend_params = cosine_params,
#                                       mod_fun = generate_cosine_series)
# 
# 
# # Flatten the list and combine into a single tibble
# test_QR_df <-  result_QR %>%
#   map_dfr(~ {
#     bind_rows(.x, .id = "period")
#   }, .id = "scaling_factor") 
# 
# test_QR_df %>% 
#   ggplot(aes(x = period, y = estimate, colour = ConfCat))+
#   geom_point(aes(size = 1),show.legend = F)+
#   geom_linerange(aes(ymin = conf.low, ymax = conf.high, x = period), size = 1.5)+
#   geom_hline(yintercept = 0, linetype = "dashed", size = .3)+
#   theme_bw()+
#   scale_color_manual(values = color_mapping) +  # Use the custom color mapping
#   labs(title = "Estimates with Confidence Intervals",
#        x = "Period",
#        y = "Estimate",
#        color = "Confidence Category")+
#   facet_wrap(~scaling_factor)

