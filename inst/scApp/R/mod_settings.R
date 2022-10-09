#' scFeatureView UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_settings_ui <- function( id , title ){
    ns <- NS( id )
    
    tabPanel(
      title = title,
      sidebarPanel(
        ## Item sidepanel
        helpText(
          paste0(
            "Select Colors"
          )
        ),
      
      ),
      mainPanel(
        theme <- bslib::bs_theme(
          bg = "#d3d3d3",
          fg = "#EEE8D5",
          primary = "#2AA198" #,
          # bslib also makes it easy to import CSS fonts
          #base_font = bslib::font_google("Pacifico")
        )
      )
    )
}
    
#' Settings Server Functions
#'
#' @noRd 
mod_settings_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
  })
}
    
## To be copied in the UI
# mod_scFeatureView_ui("scFeatureView_1")
    
## To be copied in the server
# mod_scFeatureView_server("scFeatureView_1")
