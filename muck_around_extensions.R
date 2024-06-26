
library(quantreg)

 cosine_params <-  list(
                        decay_rate = 0.01,
                        initial_amplitude = 2,
                        num_peaks = 2,
                        phase_shift = 0
 )

 test_rol_td <- analyze_trend_rolling(site_data$`GW-00002`,
                                periods = list("full_length",c(15,10), c(15,5), c(15,0), c(10,5), c(10,0), c(8,3), c(6,1), c(5,0)),
                                is_seasonal = TRUE,
                                trend_params = cosine_params,
                                mod_fun = generate_cosine_series)
 
 result <-  test_rol_td
 
 
 lin_params <- list(slope = 0.0025, intercept = 1)
 
 test_rol_td_lin <- analyze_trend_rolling(site_data$`GW-00002`,
                                      periods = list("full_length",c(15,10), c(15,5), c(15,0), c(10,5), c(10,0), c(8,3), c(6,1), c(5,0)),
                                      is_seasonal = TRUE,
                                      trend_params = lin_params,
                                      mod_fun = generate_linear_trend)
 
 
 test_rol_td_lin <- analyze_trend_rolling(site_data$`GW-00002`,
                                          periods = list("full_length",c(15,10), c(15,5), c(15,0), c(10,5), c(10,0), c(8,3), c(6,1), c(5,0)),
                                          is_seasonal = TRUE,
                                          trend_params = lin_params,
                                          mod_fun = NULL)

result <-  test_rol_td_lin

#########

 MK <-  result$estimate_results[[1]]$MK
 names(MK) <- result$estimate_results[[1]]$period
 MK_df <- imap_dfr(MK, ~tibble(period = .y,.x))
 
 
# from cawthron  
 # add confidence classes
 trend_tab <-  MK_df %>%
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
       TrendScore = ifelse(is.na(TrendScore), NA, TrendScore)
    )
## 
 
 # Define the colors for each category
 color_mapping <- c(
    "Very likely improving" = "green4",
    "Likely improving" = "seagreen2",
    "Indeterminate" = "blue",
    "Likely degrading" = "orange",
    "Very likely degrading" = "red"
 )
 
 
 equiv_zone <-  c(-.015, .015)
 
 
 get_MK_plot <- function(MKdata, equiv_zone = c(-.015, .015)){
    
 MK_plot <- MKdata  %>% 
   ggplot(aes(colour = ConfCat))+
 geom_point(aes(y = AnnualSenSlope, x = period,size = 1),show.legend = F)+
    geom_linerange(aes(ymin = Sen_Lci, ymax = Sen_Uci, x = period), size = 1.5)+
    geom_hline(yintercept = 0, linetype = "dashed", size = .3)+
   theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    # Rectangle for the equivalence zone
    geom_rect(aes(ymin = min(equiv_zone), ymax = max(equiv_zone), xmin = -Inf, xmax = Inf),
              fill = "grey80", alpha = 0.05, colour = NA) +
    scale_color_manual(values = color_mapping, drop = FALSE)

 return( MK_plot)   
}
 
 #Add site and measurement...
 get_MK_plot(trend_tab)
## 
 
 

#equivalence
example_data <- site_data$`GW-00002`

example_result <- NonSeasonalTrendAnalysis(as.data.frame(example_data))

slope_df <-  tibble(est_slope = example_result$AnnualSenSlope / 12,
                    lci = example_result$Sen_Lci/12,
                    uci =  example_result$Sen_Uci/12)


CI_range <- c(slope_df$lci, slope_df$uci)
equiv_zone <- c(-.01,.01)



# Function to check the position of the CI relative to the equivalence region
check_equivalence <- function(CI, equiv) {
   lower_CI <- CI[1]
   upper_CI <- CI[2]
   lower_equiv <- equiv[1]
   upper_equiv <- equiv[2]
   
   if (lower_CI >= lower_equiv && upper_CI <= upper_equiv) {
      return("Fully within the equivalence region")
   } else if (upper_CI < lower_equiv || lower_CI > upper_equiv) {
      return("Fully outside the equivalence region")
   } else {
      return("Partly within the equivalence region")
   }
}

# Check the CI against the equivalence region
result <- check_equivalence(CI_range, equiv_zone)

#  PLOT
slope_df %>% 
   ggplot() +
   geom_point(aes(x = est_slope, y =.6))+
   geom_linerange(aes(xmin = lci, xmax = uci, y = .6))+
   # Rectangle for the equivalence zone
   geom_rect(aes(xmin = min(equiv_zone), xmax = max(equiv_zone), ymin = 0, ymax = 1),
             fill = "grey80", alpha = 0.5) +
   theme_classic()+
   coord_cartesian(xlim = c(-.02,.02))
# Horizontal error bar for the CI


#####

trend_tab %>% 
   ggplot() +
   geom_point(aes(x = est_slope, y =.6))+
   geom_linerange(aes(xmin = lci, xmax = uci, y = .6))+
   # Rectangle for the equivalence zone
   geom_rect(aes(xmin = min(equiv_zone), xmax = max(equiv_zone), ymin = 0, ymax = 1),
             fill = "grey80", alpha = 0.5) +
   theme_classic()+
   coord_cartesian(xlim = c(-.02,.02))

########
test_stl <- Get_STL(site_data$`GW-00004`)
test_rolling <- rolling_trend(test_stl, periods = list("full_length", c(5,0)))

test_rol_wrap <- analyze_trend_rolling(site_data$`GW-00004`,periods = list("full_length", c(5,0))) 

undebug(rolling_trend)

yearrange

# Concatenate the elements
yearrange_ccat <- paste(yearrange[1], yearrange[2], sep = ":")


#########

#MK results
MK <- result$estimate_results[[1]]$MK
names(MK) <- result$estimate_results[[1]]$period
MK_df <- imap_dfr(MK, ~tibble(period = .y,.x))

MK_df <-  get_ConfCat(MK_df)
get_MK_plot(MK_df)


rolling_results <- test_rol_wrap


test_cc <-   get_ConfCat(test_rol_wrap) 
get_MK_plot(test_cc) 


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

######
scaling_factor = list(c(1,1), c(0,2))
scale_comp_test <- scale_components(stl_data, scaling_factor = scaling_factor )

components(scale_comp_test$data[[1]]) %>% autoplot()
components(scale_comp_test$data[[2]]) %>% autoplot()


stl_data <- Get_STL(site_data$`GW-00002`)
stl_data <- test_stl
scaling_factor = list(c(3,2), c(4,0))
lambda_value <-scaling_factor[[1]]
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









senslope_res_list[[1]]$data[[1]] %>% 
   components() %>% autoplot()


####
get_QR_est <- function(data, trm = "Month", form = formula(paste0("RawValue~Year+", trm))) {
   
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
   
   
   
   qr_med_sum <- summary(qr_med, se = "boot", covariance = TRUE, R = 300) 
   
   # Extract the coefficients table from the summary
   coef_table <- qr_med_sum$coefficients
   
   # Extract the lower and upper bounds of the confidence intervals
   coef_table_df <- data.frame(
      intercept = coef_table[1, 1],
      lci_int = coef_table[1,1] - 1.96 * coef_table[1, 2],
      uci_int = coef_table[1, 1] + 1.96 * coef_table[1, 2],
      trend = coef_table[2, 1],
      trend_CI = 1.96 * coef_table[2, 2],
      lci_trend = (coef_table[2,1] - 1.96 * coef_table[2, 2]),
      uci_trend = (coef_table[2, 1] + 1.96 * coef_table[2, 2])
   )
   
   
   # Example data creation for plotting (modify according to your actual data)
   years <- seq(min(data$Year), max(data$Year), by = 1)
   predicted_df <- data.frame(Year = years)
   
   
   # Calculate predicted values using the coefficients from your model output
   predicted_df$Predicted_Value <- coef_table_df$intercept  + (coef_table_df$trend * predicted_df$Year)
   predicted_df$trend_CI  <- coef_table_df$trend_CI 
   
   
   
   
   data$fitted <-  predict(qr_med, newdata = data)
   
   QR_plot <-  ggplot( data, aes(x = Year)) +
      geom_point(aes(y = RawValue, color = "Data Points"), alpha = 0.6) +  
      geom_point(aes(y = fitted, color = "Fitted"), alpha = 0.6, ) +  
      geom_line(data = predicted_df,
                aes(y = Predicted_Value, 
                    x = Year, 
                    group = 1,
                    color = "Quantile Regression"), size = 1) +
      geom_ribbon(data = predicted_df,aes(
         ymin = Predicted_Value - trend_CI, 
         ymax = Predicted_Value + trend_CI),
         alpha = .4)
   
   
   
   return(list(qr_med,coef_table_df,QR_plot ))
   
}


scale_stl_noise_QR(stl_data) 

scale_stl_noise_QR <- function(stl_data, lambda = seq(from = -1, to = 3, by = 0.5), is_seasonal = TRUE, analysis_params = list()) {
   
   # Container lists
   QR_results <- list()
   complist <- list()
   senslope_res_list <- list()
   
   comp_df <- components(stl_data) %>% as_tibble()
   
   #needs extra columns for LWP functions to work
   #   orig_data <- stl_data$orig_data[[1]] %>% 
   #     select(lawa_site_id, CenType,Censored, 
   #            yearmon,Season, Year, myDate,RawValue)
   # 
   # comp_df <- comp_df %>% left_join(orig_data, by = "yearmon")
   
   
   
   ## Simulate and replace trend here
   # comp_df$trend
   ##
   
   for (lambda_value in lambda) {
      
      
      comp_df <- comp_df %>% 
         mutate(  simulated = final_series + (remainder * lambda_value),
                  RawValue = simulated  ) 
      
      get_QR_est <- safely(get_QR_est)
      
      QR_est <- get_QR_est(comp_df)
      #comp_df$RawValue <-  comp_df$final_series
      
      
      QR_results[[as.character(lambda_value)]] <- list(QR_est)
      
      # if (is_seasonal == TRUE){
      #   senslope_res <- do.call(SeasonalTrendAnalysis, c(list(as.data.frame(comp_df)), analysis_params))  # Pass additional arguments
      # } else {
      #   senslope_res <- do.call(NonSeasonalTrendAnalysis, c(list(as.data.frame(comp_df)), analysis_params))  # Pass additional arguments
      # }
      # 
      # slope_df <- tibble(est_slope = senslope_res$AnnualSenSlope / 12,
      #                    lci = senslope_res$Sen_Lci / 12,
      #                    uci = senslope_res$Sen_Uci / 12,
      #                    CI_width =uci - lci)
      
      #senslope_res_list[[as.character(lambda_value)]] <- slope_df
   }
   
   #senslope_res_list <- imap_dfr(senslope_res_list, ~ tibble(lambda = .y, .x))
   
   #estimate_results <- senslope_res_list %>% 
   #  mutate(CI_width = uci - lci)
   
   return(QR_results)
}


#incorporate scaline with rolling


undebug(analyze_trend_rolling_alt)

analyze_trend_rolling_alt(data = site_data$`GW-00002`)




undebug(scale_components)
undebug(rolling_trend_alt)


test_sc_rol <- scale_components(stl_data),scaling_factor = list(c(1,1), c(1,2), c(2,1)))


test_sc_rol[[3]]

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
      
      #replace STL components
      stl_data_mod <- stl_data
      stl_data_mod$stl[[1]]$fit$decomposition <- comp_df
      
      rol_est <- rolling_trend_alt(stl_data_mod,is_seasonal = is_seasonal)
      
      rol_est[[1]]$seasonal_sf <- scaling_factor[[i]][1]
      rol_est[[1]]$noise_sf = scaling_factor[[i]][2]
   
   
      senslope_res_list[[i]] <- rol_est
   }
    return(senslope_res_list)
   }  


#######


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

test_rt <- rolling_trend(stl_data) 

test_rt[[1]]$seasonal_sf <- scaling_factor[[i]][1]
test_rt[[1]]$noise_sf = scaling_factor[[i]][2]




rolling_trend_alt <- function(stl_data, periods = list("full_length", c(5, 0)), 
                          is_seasonal = TRUE, analysis_params = list(),
                          ...) {
   
   # Container lists
   complist <- list()
   senslope_res_list <- list()
   
   comp_df <- components(stl_data) %>% as_tibble()
   
   #Removed getting of orig data, as already from other function.
   
   #needs this name for LWP functions
   #comp_df$RawValue <- comp_df$final_series
   
   
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

#
analyze_trend_rolling_alt <- function(data, 
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
      estimate_results <- scale_components(stl_data,  periods = periods, is_seasonal = is_seasonal, analysis_params = analysis_params)
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
      #estimate_results <- rolling_trend_alt(stl_data, periods = periods, is_seasonal = is_seasonal, analysis_params = analysis_params) 
      estimate_results <- scale_components(stl_data,  periods = periods, is_seasonal = is_seasonal, analysis_params = analysis_params)
      
   }
   
   # Step 4: Return the results
   #return(list(stl_data = stl_data, estimate_results = estimate_results))
   return(estimate_results)
}

