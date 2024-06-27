# ui.R
#library(shiny)

ui <- fluidPage(
  titlePanel("Time Series Simulation and Trend Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("site", "Select Site Data:", choices = NULL),
      selectInput("series_type", "Select Trend Type:",
                  choices = c("NULL", "Cosine", "Linear Trend", "Level Shift With Ramp", "Linex")),
      checkboxInput("is_seasonal", "Is Seasonal", TRUE),

      conditionalPanel(
        condition = "input.series_type == 'Cosine'",
        sliderInput("initial_amplitude", "Initial Amplitude:", min = -5, max = 5, value = .5, step = .005),
        sliderInput("decay_rate", "Decay Rate:", min = -.1, max = .1, value = 0.01, step = 0.005),
        sliderInput("num_peaks", "Number of Peaks:", min = -5, max = 10, value = 3, step = 0.25),
        sliderInput("phase_shift", "Phase Shift (radians):", min = 0, max = 2 * pi, value = pi / 4)
      ),
      
      conditionalPanel(
        condition = "input.series_type == 'Linear Trend'",
        sliderInput("slope", "Slope:", min = -1, max = 1, value = 0.1, step = 0.001),
        sliderInput("intercept", "Intercept:", min = -10, max = 10, value = 0, step = .01)
      ),
      
      conditionalPanel(
        condition = "input.series_type == 'Level Shift With Ramp'",
        sliderInput("baseline_amp", "Baseline Amplitude:", min = -5, max = 5, value = 0.3, step = 0.05),
        sliderInput("amp_change", "Amplitude Change:", min = -5, max = 5, value = 0.6, step = 0.05),
        sliderInput("ramp_start", "Ramp Start:", min = 0, max = 1, value = 0.3, step = 0.05),
        sliderInput("ramp_length", "Ramp Length:", min = 0, max = 1, value = 0.1),
        sliderInput("steepness", "Steepness:", min = 0.05, max = 2, value = 0.7, step = 0.05)
      ),
      
      conditionalPanel(
        condition = "input.series_type == 'Linex'",
        sliderInput("amplitude", "Amplitude:", min = -5, max = 5, value = 0.3, step = 0.05),
        sliderInput("x_min", "X Min:", min = -5, max = 100, value = 0.6, step = 0.05),
        sliderInput("scale", "Scale:", min = 0, max = 1, value = 0.3, step = 0.05)
      ),
      
      wellPanel(
        titlePanel("Rolling Period Parameters"),
        textInput("rolling_periods", "Enter Rolling Periods (e.g., full_length; 10,2; 5,0):", value = "full_length")
      ),
      wellPanel(
        titlePanel("Component Scaling Factors"),
        textInput("Scaling_Factors", "Enter Scaling Factors (Seasonal,Remainder): e.g., 1,1; 2,1", value = "1,1")
      ),
      
     # actionButton("estimate_noise_btn", "Estimate with Noise Scaling"),
      #actionButton("estimate_rolling_btn", "Estimate with Rolling Trend"),
      actionButton("estimate_GAM_btn", "Estimate with GAM"),
      actionButton("estimate_wrapper_btn", "Estimate with parameters"),
      

      #downloadButton("downloadData", "Download Results")
      #downloadButton("downloadRollingData", "Download Rolling Results"),
     # downloadButton("downloadNoiseData", "Download Noise Results"),
      #downloadButton("downloadPlot", "Download Plot"),
      downloadButton("downloadAll", "Download All")

    ),
    
    mainPanel(
      plotOutput("originalDataPlot"),
      plotOutput("GAMPlot"),
      plotOutput("simulatedTrendPlot"),
     # plotOutput("simcomponentsPlot"),
      #plotOutput("rolling_period_plot"),
      #plotOutput("slope_est_Plot"),
      uiOutput("WrapperPlot"),
      uiOutput("WrapperMKPlots")
      #plotOutput("WrapperPlot"),
      #plotOutput("MKplot"),
          #verbatimTextOutput("summary")
      #verbatimTextOutput("summary1"),
      #verbatimTextOutput("summary2")

    )
  )
)
