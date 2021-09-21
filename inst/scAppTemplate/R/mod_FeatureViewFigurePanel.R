#' FeatureViewFigurePanel UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_FeatureViewFigurePanel_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      textOutput("dev_text")
    ),
    fluidRow(
      column(8,
             uiOutput("multi_plot_ui")
      )
    )
  )
}
    
#' FeatureViewFigurePanel Server Functions
#'
#' @noRd 
mod_FeatureViewFigurePanel_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_FeatureViewFigurePanel_ui("FeatureViewFigurePanel_1")
    
## To be copied in the server
# mod_FeatureViewFigurePanel_server("FeatureViewFigurePanel_1")
