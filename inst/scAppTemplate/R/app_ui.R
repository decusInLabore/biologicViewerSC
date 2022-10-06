#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import ggplot2
#' @import DT
#' @import colourpicker
#' @import DBI
#' @import RMySQL
#' @import RSQLite
#' @noRd

###############################################################################
## Main App - server                                                         ##       

app_ui <- function(request) {
  
    # theme <- bslib::bs_theme(
    #   bg = "#d3d3d3", 
    #   fg = "#EEE8D5", 
    #   primary = "#2AA198" #,
    #   # bslib also makes it easy to import CSS fonts
    #   #base_font = bslib::font_google("Pacifico")
    # )
    
    #theme = bs_theme(version = 4, bootswatch = "cosmo")

    tagList(
        # Leave this function for adding external resources
        golem_add_external_resources(),
                                                                         ##       
        ## Start new
        navbarPage(
            title = "biologic SC",   
            mod_scFeatureView_ui("scFeatureView_1", title = "FeatureView"),
            
            mod_scFeatureView_ui(
                "scFeatureView_2", 
                title = "Violin Plots", 
                xSel = "clusterName", 
                ySel = "lg10Expr",
                colorSel = "clusterName"
            ),
            
            mod_scFeatureView_ui(
              "scFeatureView_3", 
              title = "Ridge Plots", 
              xSel = "lg10Expr", 
              ySel = "Ridgeplot",
              colorSel = "clusterName"
            ),
        
            mod_about_ui("about_ui_1", title = "About this Dataviewer"),
        
            navbarMenu("Settings", 
                tabPanel("Color Setting", "Color Settings"),
                tabPanel("Other Settings", "Other Settings")
            )
        ) ## End Navbar Parge
    
    )
}

##                                                                           ##
###############################################################################

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
    add_resource_path(
        'www', app_sys('app/www')
    )

    tags$head(
        tags$link(rel="shortcut icon", href="www/favicon.ico"),
        bundle_resources(
            app_title = "biologicSC",
            path = app_sys('app/www')
        ),
        tags$link(
            rel="stylesheet", 
            type="text/css", 
            href="www/custom.css"
        ) 
    )
}

