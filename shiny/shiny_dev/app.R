

#Build categorisation into it?? old vs new??? based on CI. Add that to interactive...
#colour based on the categorisation????

# Check if pacman is installed; if not, install it
# if (!requireNamespace("pacman", quietly = TRUE)) {
#    install.packages("pacman")
# }
# 
# # Load required packages using pacman
# pacman::p_load(
#    plyr,
#    NADA,
#    shiny,
#    ggplot2,
#    withr,
#    tsibble,
#    imputeTS,
#    tibble,
#    purrr,
#    dplyr,
#    fable,
#    feasts
# )

#TODO GAM, QR, Exports



#site_data <- readRDS("~/My-Local-Workspace/Github_Repos/lawa-trends/For_Andrew/N03N_filtered.rds")


site_data <- readRDS( here::here( "N03N_filtered.rds"))

# Source the UI and server files


#source(here::here("shiny/shiny_dev/server.R"))

source(here::here("shiny/shiny_dev/server_Alt.R"))


source(here::here("Gam_functions.R"))

source(here::here("shiny/shiny_dev/ui.R"))
#source("For_Andrew/shiny/shiny_dev/ui.R")
#source("For_Andrew/shiny/shiny_dev/server.R")


#options(shiny.reactlog = FALSE)

# Run the application
shinyApp(ui = ui, server = server)
