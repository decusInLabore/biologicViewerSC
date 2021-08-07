#' FeatureViewSidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
#' 

#violinURL <- a("Violin Plot", href="localhost/?_inputs_&gene=%22MYL5%22&y_axis=%22lg10Expr%22&x_axis=%22UMAP_2%22")


mod_FeatureViewSidebar_ui <- function(id){
  ns <- NS(id)
  tagList(
    helpText(paste0("To create a Violin Plot plot, select, for example, as x-axis: seurat clusters, as y-axis: log10Expr and as colorBy: seurat clusters. \n\n To view averaged expression values for signature gene categories, start typing cat_ in the search box to see category suggestions. ")),
    
    conditionalPanel(
      condition = "input.colorBy == 'lg10Expr'|| input.x_axis == 'lg10Expr' || input.y_axis == 'lg10Expr'",
      selectizeInput("gene", 
                   label = "Gene or Category Selection",
                   choices = NULL, #c(as.vector(sort(unique(allGenes)))),
                   selected = geneDefault,
                   options = list(maxOptions = 50, create=TRUE))
    ),
    
    
    
    selectInput("x_axis", 
                label = "Choose an X-axis",
                choices =unique(c("Log10 Expression" = "lg10Expr", allColorOptions, XYsel)),
                selected = "UMAP_1"),
    selectInput("y_axis", 
                label = "Choose an Y-axis",
                choices =unique(c("Log10 Expression" = "lg10Expr", XYsel)),
                selected = "UMAP_2"),
    
    
    selectInput("splitByColumn", 
                label = "Split Plots By",
                choices = splitOptions,
                selected = splitOptions[1]),
    
    selectInput("colorBy", 
                label = "Color Plots By",
                choices = allColorOptions,
                selected = names(allColorOptions)[1]),
    
    #####################################################################
    ## Query color input based on 'colorBy' selection                  ##
    conditionalPanel(
      condition = paste0("input.colorBy == '",numOptions,"'", collapse = "||"),
      colourInput("dotcolor", "Select Low Color", "darkblue"),
      colourInput("lowColor", "Select High color", "#D3D3D3")
    ),
    
    
    conditionalPanel(
      condition = paste0("input.colorBy == '",splitOptions,"'", collapse = "||"),
      
      uiOutput("clusterColorPanel")
      
    ),
    
    
    ## Done                                                            ##
    #####################################################################
    
    selectInput("background",
                label = "Select Background",
                choices =c("Grey" = "grey", "White" = "white","Minimal" = "minimal", "Plain" =  "plain"),
                selected = "white"),
    
    
    sliderInput("dotsize", "Choose a Dotsize",
                min = 0.01, max = 2, value = 1
    ),
    checkboxInput("showPlotLegend", "Show Plot Legends", value = TRUE, width = NULL),
    
    
    bookmarkButton(),
    br(),br(),
    downloadButton('plotDLall', "Download Plot Images"),
    br(),br(),
    conditionalPanel(
      condition = "input.colorBy != 'lg10Expr'",
      downloadButton("downloadData", "Download Color Selection")
    )
    
    
  )
}
    
#' FeatureViewSidebar Server Function
#'
#' @noRd 
mod_FeatureViewSidebar_server <- function(input, output, session){
  ns <- session$ns
 
}
    
## To be copied in the UI
# mod_FeatureViewSidebar_ui("FeatureViewSidebar_ui_1")
    
## To be copied in the server
# callModule(mod_FeatureViewSidebar_server, "FeatureViewSidebar_ui_1")
 
