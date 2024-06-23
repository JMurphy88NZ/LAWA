if (!("pacman" %in% rownames(installed.packages()))) {install.packages("pacman")}
pacman::p_load(pacman,
               tidyverse,
               readxl,
               lubridate,
               mgcv,
               gratia,
               furrr,
               tictoc,
               magrittr,
               beepr,
               RcppRoll,
               scales)

##
source("Gam_functions.R")

N03N_filtered <- readRDS("~/LAWA/LAWAgit/N03N_filtered.rds")


test_df <- N03N_filtered$`GW-00002`

#test_df$yearmon <- tsibble::yearmonth(test_df$yearmon)

tscm <- screeningmodeling(.data = N03N_filtered$`GW-00002` , datevar = myDate,values = RawValue )
#tscm2 <- screeningmodeling(.data = N03N_filtered$`GW-00002` , datevar = yearmon,values = RawValue )

plot_individual_trend(tscm)

dplyr::between()
plot_individual_trend(tscm2)

undebug(plot_individual_trend)
##




testdf <- N03N_filtered[[62]]

#testdf$yearmon

testdf$ymon <-  tsibble::yearmonth(testdf$myDate)
test_gam <- screeningmodeling(.data = testdf , datevar = ymon,values = RawValue )

debug(screeningmodeling)

test_gam <- screeningmodeling(.data = testdf , datevar = myDate,values = RawValue )

plot_individual_trend(test_gam)

plot_individual_trend_log(test_gam)
###


N03N_filtered_df <- map_df(N03N_filtered,~.)

N03N_filtered_df %>%
  group_by(lawa_site_id, measurement) %>% 
  #filter(measurement == "bd") %>% 
  #drop_na( value) %>% #TRY WITH NAs too
  mutate(n=length(yearmon )) %>% 
  screeningmodeling(values=  RawValue ,
                    datevar = myDate , 
                    link = "identity", 
                    conf.type = "conf",
                    conf.level=0.95,
                    beep = TRUE, 
                    tdist = F,
                    autocor = TRUE,
                    lawa_site_id, measurement) -> 
  LAWA_screenmodel_out




screen_plot <- LAWA_screenmodel_out %>%
  #arrange(desc(SWEREFID))%>%
  plot_screeningtrends(sorting =lawa_site_id,
                       y_id = lawa_site_id,
                       wrappingvar = measurement)

screen_plot$data
#plotly::ggplotly(screen_plot)

LAWA_screenmodel_out$fderiv[[1]] %>% tail()



LAWA_screenmodel_out$fderiv[[1]] %>%
  mutate(deriv = .derivative, 
         deriv_se = .se, 
         lower=.lower_ci, 
         upper=.upper_ci) %>%
  as_tibble %>%
  rowwise %>%
  mutate(signif = !between(0, lower, upper),
         sign=sign(deriv),
         signif_sign = signif*sign) 


#get info on significance across series

get_signif <- function(mod){
  
  signif_ts <- mod$fderiv[[1]] %>%
    mutate(deriv = .derivative, 
           deriv_se = .se, 
           lower=.lower_ci, 
           upper=.upper_ci) %>%
    as_tibble %>%
    rowwise %>%
    mutate(signif = !between(0, lower, upper),
           sign=sign(deriv),
           signif_sign = signif*sign) 
  return(signif_ts)
}


sig_ts <- get_signif(LAWA_screenmodel_out[1,])

 
sig_ts %>% 
  group_by(signif,sign) %>% 
  
##################

#test GAM with simulate...

###
#get_GAM_trend(stl_data = stl_data)

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
  

testGAM <-   analyze_GAM_wrapper(data =  N03N_filtered$`GW-00002`)

cosine_params <-  list(
                       decay_rate = 0.01,
                       initial_amplitude = 5,
                       num_peaks = 2,
                       phase_shift = 0
)

# test_rol_td <- analyze_trend_rolling(N03N_filtered$`GW-00002`,
#                                periods = list("full_length", c(10,0), c(8,3), c(6,1), c(5,0)),
#                                is_seasonal = TRUE,
#                                trend_params = cosine_params,
#                                mod_fun = generate_cosine_series)

testGAMcos <-   analyze_GAM_wrapper(data =  N03N_filtered$`GW-00002`,
                                     trend_params = cosine_params,
                                    mod_fun = generate_cosine_series)

