#' The application server-side
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import DBI
#' @import RMySQL
#' @import ggplot2
#' @import colourpicker
#' @import scales
#' @import ggridges




###############################################################################
## Main App Server                                                           ##

app_server <- function(input, output, session) {
  
    ## FeatureView Module ##
    mod_scFeatureView_server("scFeatureView_1")
    ## Violin Plots > FeatureView Module ##
    mod_scFeatureView_server("scFeatureView_2")
    ## Ridge Plots > FeatureView Module ##
    mod_scFeatureView_server("scFeatureView_3")
    ## About Module ##
    mod_about_server("about_ui_1")
  
}
## Done                                                                      ##
###############################################################################