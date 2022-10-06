# Module UI
  
#' @title   mod_about_ui and mod_about_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_about
#'
#' @keywords internal
#' @export 
#' @importFrom shiny NS tagList includeMarkdown
mod_about_ui <- function(id, title = "About"){
  ns <- NS(id)
  
  tabPanel(
        title = title,
        col_6(
          shiny::includeMarkdown(
            system.file("app/www/about.md", package = "biologicViewerSC")
          )
        )
      
  )
}
    
# Module Server
    
#' @rdname mod_about
#' @export
#' @keywords internal
  
mod_about_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
  })
}
    
## To be copied in the UI
# 
    
## To be copied in the server
# callModule(mod_about_server, "about_ui_1")
 
