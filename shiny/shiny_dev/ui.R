# ui.R
library(shiny)

ui <- fluidPage(
  titlePanel("Time Series Simulation and Trend Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("site", "Select Site Data:", choices = NULL),
      selectInput("series_type", "Select Trend Type:",
                  choices = c("NULL", "Cosine", "Linear Trend", "Level Shift With Ramp", "Linex")),
      
      wellPanel(
        titlePanel("Lambda Parameters"),
        sliderInput("lambda_min", "Lambda Min:", min = -5, max = 5, value = -1, step = 0.1),
        sliderInput("lambda_max", "Lambda Max:", min = -5, max = 5, value = 3, step = 0.1),
        sliderInput("lambda_step", "Lambda Step:", min = 0.01, max = 1, value = 0.5, step = 0.01)
      ),
      
      conditionalPanel(
        condition = "input.series_type == 'Cosine'",
        sliderInput("initial_amplitude", "Initial Amplitude:", min = 0, max = 10, value = 5),
        sliderInput("decay_rate", "Decay Rate:", min = -1, max = 1, value = 0.01, step = 0.005),
        sliderInput("num_peaks", "Number of Peaks:", min = -5, max = 10, value = 3, step = 0.25),
        sliderInput("phase_shift", "Phase Shift (radians):", min = 0, max = 2 * pi, value = pi / 4)
      ),
      
      conditionalPanel(
        condition = "input.series_type == 'Linear Trend'",
        sliderInput("slope", "Slope:", min = -1, max = 1, value = 0.1, step = 0.01),
        sliderInput("intercept", "Intercept:", min = -10, max = 10, value = 0)
      ),
      
      conditionalPanel(
        condition = "input.series_type == 'Level Shift With Ramp'",
        sliderInput("baseline_amp", "Baseline Amplitude:", min = -5, max = 5, value = 0.3, step = 0.05),
        sliderInput("amp_change", "Amplitude Change:", min = -5, max = 5, value = 0.6, step = 0.05),
        sliderInput("ramp_start", "Ramp Start:", min = 0, max = 1, value = 0.3, step = 0.05),
        sliderInput("ramp_length", "Ramp Length:", min = 0, max = 1, value = 0.1),
        sliderInput("steepness", "Steepness:", min = -2, max = 2, value = 0.8, step = 0.1)
      ),
      
      conditionalPanel(
        condition = "input.series_type == 'Linex'",
        sliderInput("amplitude", "Amplitude:", min = -5, max = 5, value = 0.3, step = 0.05),
        sliderInput("x_min", "X Min:", min = -5, max = 100, value = 0.6, step = 0.05),
        sliderInput("scale", "Scale:", min = 0, max = 1, value = 0.3, step = 0.05)
      ),
      

      
      actionButton("estimate_btn", "Estimate")
    ),
    
    mainPanel(
      plotOutput("originalDataPlot"),
      plotOutput("simulatedTrendPlot"),
      plotOutput("timeSeriesPlot"),
      plotOutput("simcomponentsPlot"),
      plotOutput("slope_est_Plot"),
      verbatimTextOutput("summary")
    )
  )
)
