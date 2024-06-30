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
  
  #functions to clear outputs
  # clear_gam_outputs <- function() {
  #   output$gamPlot <- renderPlot({ NULL })
  #   output$gamSummary <- renderPrint({ NULL })
  # 
  
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
  

  
  #GAM
  
  observeEvent(input$estimate_GAM_btn, {
    
    shiny::showNotification("Started GAM analysis", type = "message")
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
    
    
    
    
    # Use tryCatch to handle any errors during the GAM analysis
    GAMresult <- tryCatch({
      analyze_GAM_wrapper(input_data, trend_params = trend_params, mod_fun = mod_fun)
    }, error = function(e) {
      # Display a user-friendly message
      showNotification("An error occurred during the GAM analysis.", type = "error")
      NULL  # Return NULL if an error occurs
    })
    

    
    output$GAMPlot <- renderPlot({
      if (!is.null(GAMresult)) {
        GAMresult[[2]]
      } else {
        plot.new()
        text(0.5, 0.5, "No data to display", cex = 1.5)
      }
      
    })
    


  })
  
 
    
  # SENSLOPE/MK WRAPPER
  
  # Reactive value to store result
  MK_result <- reactiveVal(NULL)
  
   observeEvent(input$estimate_wrapper_btn, {
     
    # browser()
     
     shiny::showNotification("Started estimate trend wrapper event", type = "message")
     
     trend_data <- simulate_trend()
     
     
     scaling_input <- strsplit(input$Scaling_Factors, ";")[[1]]
     
     scaling_factors <-  lapply(scaling_input, function(x){
                               sf <- as.numeric(unlist(strsplit(x, ",")))
                               return(sf)
                               })
     
     
       
       
       
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
       
       
     #params for modifying trend using simulate functions
       
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
     
    
     
     
     MK_result <- tryCatch({analyze_trend_wrapper(input_data, 
                                               periods = periods, 
                                               scaling_factor = scaling_factors,
                                               is_seasonal = TRUE,
                                               trend_params = trend_params,
                                               mod_fun = mod_fun)
     },error = function(e) {
       # Display a user-friendly message
       showNotification("An error occurred.Check period input format", type = "error")
       NULL  # Return NULL if an error occurs
     })
     
     
     # Update the reactive value
     MK_result(MK_result)
     
     
     output$WrapperPlot <- renderUI({
       plot_list <- map(MK_result , ~ {seasonal_sf <- unique(.x[[1]]$seasonal_sf)
       noise_sf <- unique(.x[[1]]$noise_sf)

       plot <- .x[[2]]+
         labs(title = paste("Seasonal SF:", seasonal_sf, ",", " Noise SF: ",noise_sf ))

       return(plot)
       })

       plot_outputs <- lapply(seq_along(plot_list), function(i) {
         plotname <- paste0("plot", i)
         output[[plotname]] <- renderPlot({ plot_list[[i]] })
         plotOutput(plotname, width = "100%", height = "400px") # Adjust height as needed
       })

       do.call(tagList, plot_outputs)

     })
     # 
     # 
     output$WrapperMKPlots <- renderUI({

       MK_data_list <-  map(MK_result , ~ {seasonal_sf <- unique(.x[[1]]$seasonal_sf)
                                         noise_sf <- unique(.x[[1]]$noise_sf)
                                         #MK_data <- imap_dfr(MK_data, ~tibble(period = .y,.x))
                                         MK_data <- .x[[1]]$MK
                                         names(MK_data) <- .x[[1]]$period
                                         MK_data <- imap_dfr(MK_data, ~tibble(period = .y,.x))


                                         MK_data$seasonal_sf <- seasonal_sf
                                         MK_data$noise_sf <- noise_sf

                                         return(MK_data)
                                         })

       MK_plotlist <-  map( MK_data_list, ~ {seasonal_sf <- unique(.x$seasonal_sf)
                                        noise_sf <- unique(.x$noise_sf)

                                  get_ConfCat_from_MK_data(.) %>%
                                        get_MK_plot(.) +
                                        labs(title = paste("Seasonal SF:",
                                                            seasonal_sf, ",",
                                                            " Noise SF: ",noise_sf ) )
                                              } )


       MKplot_outputs <- lapply(seq_along(MK_plotlist), function(i) {
         MKplotname <- paste0("MKplot", i)
         output[[MKplotname]] <- renderPlot({ MK_plotlist[[i]] })
         plotOutput(MKplotname, width = "100%", height = "400px") # Adjust height as needed
       })

       do.call(tagList, MKplot_outputs)
     })
     # 

     # output$summary1 <- renderPrint({
     #   result[[1]][[1]]
     # })
     # 
     # output$summary2 <- renderPrint({
     #   result[[2]][[1]]
     # })
    
     
     
   }) 
 
   
   observeEvent(input$estimate_wrapper_QR_btn, {
     
     # browser()
     
     shiny::showNotification("Started QR estimate", type = "message")
     
     trend_data <- simulate_trend()
     
     
     scaling_input <- strsplit(input$Scaling_Factors, ";")[[1]]
     
     scaling_factors <-  lapply(scaling_input, function(x){
       sf <- as.numeric(unlist(strsplit(x, ",")))
       return(sf)
     })
     
     
     
     
     
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
       
       
       #params for modifying trend using simulate functions
       
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
     
     
     
     
     QR_result <- tryCatch({analyze_trend_wrapper_QR(input_data, 
                                               periods = periods, 
                                               scaling_factor = scaling_factors,
                                               is_seasonal = TRUE,
                                               trend_params = trend_params,
                                               mod_fun = mod_fun)
     },error = function(e) {
       # Display a user-friendly message
       showNotification("An error occurred.Check period input format", type = "error")
       NULL  # Return NULL if an error occurs
     })
     
     
     
     output$QRPlots <- renderPlot({
       
       QR_df <-  QR_result %>%
         map_dfr(~ {
           bind_rows(.x, .id = "period")
         }, .id = "scaling_factor")

       QR_df %>%
         ggplot(aes(x = period, y = estimate, colour = ConfCat))+
         geom_point(aes(size = 1),show.legend = F)+
         geom_linerange(aes(ymin = conf.low, ymax = conf.high, x = period), size = 1.5)+
         geom_hline(yintercept = 0, linetype = "dashed", size = .3)+
         theme_bw()+
         scale_color_manual(values = color_mapping) +  # Use the custom color mapping
         labs(title = "Estimates with Confidence Intervals",
              x = "Period",
              y = "Estimate",
              color = "Confidence Category")+
         facet_wrap(~scaling_factor)
       
       
     })
     
   })
     
   # DOWNLOAD HANDLER
   output$downloadAll <- downloadHandler(
     
     
     filename = function() {
       paste("MK_results-", Sys.Date(), ".zip", sep = "")
     },
     content = function(file) {
       # Create a temporary directory
       temp_dir <- tempdir()
       
       # Shorten file paths for the CSV files
       MK_file <- file.path(temp_dir, "MK_results.csv")
       MK_data_rds<- file.path(temp_dir, "MK_data.rds")
       
       QR_file <- file.path(temp_dir, "QR_results.csv")
       
       # orig_data_file <- file.path(temp_dir, "orig_data.csv")
       # mod_data_file <- file.path(temp_dir, "sim_data.csv")
       # plot_file <- file.path(temp_dir, "slope_plot.png")
       
       
       # MK results
       
       req(MK_result())
       
       MK_list <- map(MK_result(), ~ select(.[[1]], 1:4, MK) %>% unnest(MK) )
       MK_df <-  map_df(MK_list, ~ tibble(.))
       
       write.csv(MK_df, MK_file, row.names = FALSE)
       
       #data estimates are baesd on
       MK_data_list <- map(MK_result(), ~ select(.[[1]], 1:4, data) %>% unnest(data) )
       
       saveRDS(MK_data_list, MK_data_rds)
       
       
       # QR results
       # QR_df <-  QR_result %>%
       #   map_dfr(~ {
       #     bind_rows(.x, .id = "period")
       #   }, .id = "scaling_factor")
       # 
       # QR_df %>% write_csv(.,  QR_file, row.names = FALSE)
       
       
       # GAM Result
       
       # # Generate the data CSVs
       # orig_data <- result$stl_data$orig_data[[1]]
       # write.csv(orig_data, orig_data_file, row.names = FALSE)
       # 
       # mod_data <- components(result$stl_data)
       # write.csv(mod_data, mod_data_file, row.names = FALSE)
       # 
       # # Save the plot
       # ggsave(plot_file, plot = result$estimate_results[[2]], device = "png")
       
       # Create a zip file
       zip::zipr(file, files = c(MK_file, MK_data_rds))
       
       
     }
   )
   

}
