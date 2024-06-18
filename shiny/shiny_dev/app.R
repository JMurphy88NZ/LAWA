library(shiny)


#site_data <- readRDS("~/My-Local-Workspace/Github_Repos/lawa-trends/For_Andrew/N03N_filtered.rds")


site_data <- readRDS( here::here( "N03N_filtered.rds"))

# Source the UI and server files

source(here::here("shiny/shiny_dev/ui.R"))
#source(here::here("shiny/shiny_dev/server.R"))

source(here::here("shiny/shiny_dev/server_Alt.R"))

#source("For_Andrew/shiny/shiny_dev/ui.R")
#source("For_Andrew/shiny/shiny_dev/server.R")


options(shiny.reactlog = FALSE)

# Run the application
shinyApp(ui = ui, server = server)
