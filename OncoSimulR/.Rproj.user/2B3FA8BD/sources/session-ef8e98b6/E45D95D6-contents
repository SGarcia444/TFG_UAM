library(shiny)
library(OncoSimulR)
library(readxl)
library(DT)
library(periscope)
library(BiocManager)
library(rsconnect)
source("OncoSimulIndivServer.R")
source("OncoSimulIndivUI.R")
source("OncoSimulPopUI.R")
source("OncoSimulPopServer.R")

# Define UI ----
ui <- navbarPage(
  "ONCOSIMULR",
  includeCSS("css/oncoSimulRPretty.css"),
  tabPanel("OncoSimulIndiv", tabPanelOncoSimulIndiv),
  tabPanel("OncoSimulPop", tabPanelOncoSimulPop),
  
)


# Define server logic ----
server <- function(input, output) {
  # Logic functions that manage how inputs and outputs interact
  oncoSimulIndivLogic(input, output)
  
  # Add new functions Logic ++++++++
  oncoSimulPopLogic(input,output)
}


# Run the app ----
shinyApp(ui = ui, server = server)

