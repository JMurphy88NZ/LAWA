#source("For_Andrew/shiny/shiny_dev/shiny_simulate_trends.R")
#source("For_Andrew/shiny/shiny_dev/shiny_replace_trends.R")

source(here::here("LWPTrends_v2102.R"))

source(here::here("shiny/shiny_dev/shiny_simulate_trends.R"))
source(here::here("shiny/shiny_dev/shiny_replace_trends.R"))




library(shiny)
library(ggplot2)
library(withr)
library(tsibble)
library(withr)
library(imputeTS)
library(tibble)
library(purrr)
library(dplyr)
library(fable)
library(feasts)
#library(shiny)
#library(tsibble)
#library(dplyr)
#library(fable)
#library(imputeTS)
#library(feasts)
#library(ggplot2)

# Load data


server <- function(input, output, session) {
  updateSelectInput(session, "site", choices = names(site_data))
  
  data <- reactive({
    req(input$site)
    site_data[[input$site]]
  })
  
  # Perform STL decomposition on original data and plot components
  output$originalDataPlot <- renderPlot({
    req(data())
    
    input_data <- data()
    
    # Ensure necessary columns are present
    if (!all(c("yearmon", "RawValue", "myDate") %in% names(input_data))) {
      stop("Data must contain columns 'yearmon', 'RawValue', and 'myDate'")
    }
    
    # Perform STL decomposition
    stl_data <- Get_STL(input_data)
    
    # Plot STL components
    autoplot(components(stl_data$stl[[1]]$fit))
  })
  
  # Function to generate simulated trend based on selected parameters
  simulate_trend <- reactive({
    req(data())
    req(input$series_type)
    
    if (input$series_type == "NULL") {
      #mod_fun <- as.null
      return(NULL)
    }
    
    input_data <- data()
    total_length <- nrow(input_data)
    
    trend_params <- list()
    
    if (input$series_type == "Cosine") {
      trend_params <- list(
        total_length = total_length,
        initial_amplitude = input$initial_amplitude,
        decay_rate = input$decay_rate,
        num_peaks = input$num_peaks,
        phase_shift = input$phase_shift
      )
      trend_function <- generate_cosine_series
    } else if (input$series_type == "Linear Trend") {
      trend_params <- list(
        total_length = total_length,
        slope = input$slope,
        intercept = input$intercept
      )
      trend_function <- generate_linear_trend
    } else if (input$series_type == "Level Shift With Ramp") {
      trend_params <- list(
        total_length = total_length,
        baseline_amp = input$baseline_amp,
        amp_change = input$amp_change,
        ramp_start = input$ramp_start,
        ramp_length = input$ramp_length,
        steepness = input$steepness
      )
      trend_function <- generate_level_shift_with_ramp
    } else if (input$series_type == "Linex") {
      trend_params <- list(
        total_length = total_length,
        amplitude = input$amplitude,
        x_min = input$x_min,
        scale = input$scale
      )
      trend_function <- generate_linex
    } else {
      stop("Invalid modify_trend option")
    }
    
    simulated_trend <- do.call(trend_function, trend_params)
    return(list(input_data = input_data, simulated_trend = simulated_trend))
  })
  
  output$simulatedTrendPlot <- renderPlot({
    trend_data <- simulate_trend()
    if (is.null(trend_data)) {
      ggplot() + 
        labs(title = "No Trend Modification", x = "Time", y = "Value")
    } else {
      ggplot(trend_data$input_data, aes(x = yearmon)) +
        geom_line(aes(y = trend_data$simulated_trend), color = "blue") +
        labs(title = paste("Simulated Trend:", input$series_type),
             x = "Time",
             y = "Value")
    }
  })
  
  observeEvent(input$estimate_noise_btn, {
    trend_data <- simulate_trend()
    lambda <- seq(input$lambda_min, input$lambda_max, by = input$lambda_step)
    
    if (is.null(trend_data)) {
      input_data <- data()
      mod_fun <- as.null()
    } else {
      input_data <- trend_data$input_data
      
      #trend_params <- list(modify_trend = input$series_type)
      
      if (input$series_type == "Cosine") {
        trend_params <- list()
        trend_params$initial_amplitude <- input$initial_amplitude
        trend_params$decay_rate <- input$decay_rate
        trend_params$num_peaks <- input$num_peaks
        trend_params$phase_shift <- input$phase_shift
      } else if (input$series_type == "Linear Trend") {
        trend_params <- list()
        trend_params$slope <- input$slope
        trend_params$intercept <- input$intercept
      } else if (input$series_type == "Level Shift With Ramp") {
        trend_params <- list()
        trend_params$baseline_amp <- input$baseline_amp
        trend_params$amp_change <- input$amp_change
        trend_params$ramp_start <- input$ramp_start
        trend_params$ramp_length <- input$ramp_length
        trend_params$steepness <- input$steepness
      } else if (input$series_type == "Linex") {
        trend_params <- list()
        trend_params$amplitude <- input$amplitude
        trend_params$x_min <- input$x_min
        trend_params$scale <- input$scale
      }
      
      mod_fun <- switch(input$series_type,
                        "Cosine" = generate_cosine_series,
                        "Linear Trend" = generate_linear_trend,
                        "Level Shift With Ramp" = generate_level_shift_with_ramp,
                        "Linex" = generate_linex)
    }
    
    

    
    result <- analyze_trend_with_noise(input_data, 
                                       lambda = lambda, 
                                       is_seasonal = TRUE, 
                                       trend_params = trend_params,
                                       mod_fun = mod_fun)
    
    message("finished analyze_trend_with_noise analysis")
    
    output$timeSeriesPlot <- renderPlot({
      plot_data <- result$stl_data$stl[[1]]$fit$decomposition
      
      ggplot(plot_data, aes(x = yearmon, y = final_series)) +
        geom_line() +
        labs(title = paste("Time Series with", input$series_type, "Trend"),
             x = "Time",
             y = "Value")
    })
    
    output$simcomponentsPlot <- renderPlot({
      autoplot(components(result$stl_data))
    })
    
    output$slope_est_Plot <- renderPlot({
      result$estimate_results %>% 
        ggplot(aes(x = as.numeric(lambda), y = est_slope)) +
        geom_point() +
        geom_linerange(aes(ymin = lci, ymax = uci))
    })
    
    output$summary <- renderPrint({
      result$estimate_results
    })
  })
  
  #GAM
  
  observeEvent(input$estimate_GAM_btn, {
    trend_data <- simulate_trend()
    
    if (is.null(trend_data)) {
      input_data <- data()
      mod_fun <- as.null()
    } else {
      input_data <- trend_data$input_data
      
      #trend_params <- list(modify_trend = input$series_type)
      
      if (input$series_type == "Cosine") {
        trend_params <- list()
        trend_params$initial_amplitude <- input$initial_amplitude
        trend_params$decay_rate <- input$decay_rate
        trend_params$num_peaks <- input$num_peaks
        trend_params$phase_shift <- input$phase_shift
      } else if (input$series_type == "Linear Trend") {
        trend_params <- list()
        trend_params$slope <- input$slope
        trend_params$intercept <- input$intercept
      } else if (input$series_type == "Level Shift With Ramp") {
        trend_params <- list()
        trend_params$baseline_amp <- input$baseline_amp
        trend_params$amp_change <- input$amp_change
        trend_params$ramp_start <- input$ramp_start
        trend_params$ramp_length <- input$ramp_length
        trend_params$steepness <- input$steepness
      } else if (input$series_type == "Linex") {
        trend_params <- list()
        trend_params$amplitude <- input$amplitude
        trend_params$x_min <- input$x_min
        trend_params$scale <- input$scale
      }
      
      mod_fun <- switch(input$series_type,
                        "Cosine" = generate_cosine_series,
                        "Linear Trend" = generate_linear_trend,
                        "Level Shift With Ramp" = generate_level_shift_with_ramp,
                        "Linex" = generate_linex)
    }
    
    
    
    
    GAMresult <- analyze_GAM_wrapper(input_data, 
                                    trend_params = trend_params,
                                    mod_fun = mod_fun)
    

    
    output$GAMPlot <- renderPlot({
      GAMresult[[2]]
      
    })
    


  })
  
  # Reactive values to store results
  #rolling_results <- reactiveVal(NULL)
  #noise_results <- reactiveVal(NULL)
  
  observeEvent(input$estimate_rolling_btn, {
    
    shiny::showNotification("Started estimate rolling button event", type = "message")
    
    trend_data <- simulate_trend()
    
    # Parse the periods input by ";", e.g, 'full_length; 10,2; 5,0'
    periods_input <- strsplit(input$rolling_periods, ";")[[1]]
    
    periods <- lapply(periods_input, function(period) {
      if (period == "full_length") {
        return("full_length")
      } else {
        as.numeric(unlist(strsplit(period, ",")))
      }
    })
    
    # 
     if (is.null(trend_data)) {
       input_data <- data()
       trend_params <- list()
       mod_fun <- as.null()
     } else {
       input_data <- trend_data$input_data


      
      
      # input_data <- if (is.null(trend_data)) {
      #   as.data.frame(data())
      # } else {
      #   trend_data$input_data
      # }
      # 
      # if (!is.data.frame(input_data)) {
      #   stop("input_data must be a data frame")
      # }
      
      # if (is.null(trend_data)) {
      #   trend_params <- list()
      #   mod_fun <- as.null()
      # } else {
      
      
      if (input$series_type == "Cosine") {
        trend_params <- list()
        trend_params$initial_amplitude <- input$initial_amplitude
        trend_params$decay_rate <- input$decay_rate
        trend_params$num_peaks <- input$num_peaks
        trend_params$phase_shift <- input$phase_shift
      } else if (input$series_type == "Linear Trend") {
        trend_params <- list()
        trend_params$slope <- input$slope
        trend_params$intercept <- input$intercept
      } else if (input$series_type == "Level Shift With Ramp") {
        trend_params <- list()
        trend_params$baseline_amp <- input$baseline_amp
        trend_params$amp_change <- input$amp_change
        trend_params$ramp_start <- input$ramp_start
        trend_params$ramp_length <- input$ramp_length
        trend_params$steepness <- input$steepness
      } else if (input$series_type == "Linex") {
        trend_params <- list()
        trend_params$amplitude <- input$amplitude
        trend_params$x_min <- input$x_min
        trend_params$scale <- input$scale
      }
      
      mod_fun <- switch(input$series_type,
                        "Cosine" = generate_cosine_series,
                        "Linear Trend" = generate_linear_trend,
                        "Level Shift With Ramp" = generate_level_shift_with_ramp,
                        "Linex" = generate_linex)
    
    
    }
    
    #browser()

    
    result <- analyze_trend_rolling(input_data, 
                                    periods = periods, 
                                    is_seasonal = TRUE,
                                    trend_params = trend_params,
                                    mod_fun = mod_fun)
    
    shiny::showNotification("Finished analyze_trend_rolling", type = "message")
    
    #OUTPUT
   output$timeSeriesPlot <- renderPlot({
     result$estimate_results[[2]]
   })

   output$simcomponentsPlot <- renderPlot({
     autoplot(components(result$stl_data$stl[[1]]$fit))
   })

   output$slope_est_Plot <- renderPlot({
     result$estimate_results[[1]] %>%
       ggplot(aes(x = period, y = est_slope)) +
       geom_point() +
       geom_linerange(aes(ymin = lci, ymax = uci)) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1))
   })

   output$rolling_period_plot <- renderPlot({
     result$estimate_results[[2]]
   })

   output$summary <- renderPrint({
     result$estimate_results[[1]]
   })
   
   # browser()
   # Download handlers
   # output$downloadRollingData <- downloadHandler(
   #   filename = function() {
   #     paste("rolling_results-", Sys.Date(), ".csv", sep="")
   #   },
   #   content = function(file) {
   # 
   #     rolling_results <-  result$estimate_results[[1]] %>% select(-data)
   # 
   #     write.csv(rolling_results, file, row.names = FALSE)
   #   }
   # )
  
   
   output$downloadRollingData <- downloadHandler(
     filename = function() {
       paste("rolling_results-", Sys.Date(), ".zip", sep="")
     },
     content = function(file) {
       # Create a temporary directory
       temp_dir <- tempdir()
       
       # File paths for the CSV files
       rolling_file <- file.path(temp_dir, "rolling_results.csv")
       input_data_file <- file.path(temp_dir, "input_data.csv")
       
       # Generate the rolling results CSV
       rolling_results <- result$estimate_results[[1]] %>% select(-data)
       write.csv(rolling_results, rolling_file, row.names = FALSE)
       
       # Generate the input data CSV
       write.csv(input_data, input_data_file, row.names = FALSE)
       
       # Create a zip file
       zip::zip(file, files = c(rolling_file, input_data_file))
     }
   )
   
   
  
  })
}
