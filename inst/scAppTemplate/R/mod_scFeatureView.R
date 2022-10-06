#' scFeatureView UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_scFeatureView_ui <- function( 
    id , title, 
    xSel = NULL, 
    ySel = NULL, 
    colorSel = NULL 
){
    ns <- NS( id )
    ## Load parameters ##
    startUpList <- golem::get_golem_options( which = "startUpList" )
    geneDefault <- startUpList$keyList$geneDefault
    dropDownList <- startUpList$utilityList$dropDownList
    numCols <- startUpList$utilityList$numCols
    nonNumCols <- startUpList$utilityList$nonNumCols
    splitOptions <- startUpList$utilityList$dropDownList$splitByColumn$selOptions
  
    if (is.null(xSel)){
        xSel <- as.vector(dropDownList[["x_axis"]][["default"]])
    }
    
    if (is.null(ySel)){
      ySel <- as.vector(dropDownList[["y_axis"]][["default"]])
    }
    
    if (is.null(colorSel)){
      colorSel <- as.vector(dropDownList[["colorBy"]][["default"]])
    }
  
    tabPanel(
        title = title, 
        # Sidebar panels
        sidebarPanel(
            ## Item sidepanel
            helpText(
                paste0(
                    "To create a Violin Plot plot, select, for example, as x-axis: seurat clusters, as y-axis: log10Expr and as colorBy: seurat clusters. \n\n To view averaged expression values for signature gene categories, start typing cat_ in the search box to see category suggestions. "
                )
            ),
          
            ## Item sidepanel
            conditionalPanel(
              condition = paste0(
                'input[\'', ns('colorBy'), "\'] == \'lg10Expr\' || input[\'", ns('x_axis'), "\'] == \'lg10Expr\' || input[\'", ns('y_axis'), "\'] == \'lg10Expr\'  "
              ),
              
              selectizeInput(
                ns("gene"), 
                label = as.vector(dropDownList[["gene"]][["displayName"]]),
                choices = NULL, 
                options = list(maxOptions = 50))
            ),
            ## Item sidepanel
            selectInput(
              ns("x_axis"),
              label = as.vector(dropDownList[["x_axis"]][["displayName"]]),
              choices =dropDownList[["x_axis"]][["selOptions"]],
              selected = xSel
            ),
            ## Item sidepanel
            selectInput(
              ns("y_axis"),
              label = as.vector(dropDownList[["y_axis"]][["displayName"]]),
              choices =dropDownList[["y_axis"]][["selOptions"]],
              selected = ySel
            ),
            ## Item sidepanel
            selectInput(
              ns("splitByColumn"),
              label = as.vector(dropDownList[["splitByColumn"]][["displayName"]]),
              choices =dropDownList[["splitByColumn"]][["selOptions"]],
              selected = as.vector(dropDownList[["splitByColumn"]][["default"]])
            ),
            ## Item sidepanel
            selectInput(
              ns("colorBy"),
              label = as.vector(dropDownList[["colorBy"]][["displayName"]]),
              choices =dropDownList[["colorBy"]][["selOptions"]],
              selected = colorSel
            ),
            
            ## Item sidepanel
            conditionalPanel(
              condition = paste0(paste0('input[\'', ns('colorBy'), "\'] == \'",numCols,"\'"), collapse = " || "),
              colourpicker::colourInput(ns("dotcolor"), "Select High Color", "darkblue"),
              colourpicker::colourInput(ns("lowColor"), "Select Low color", "#D3D3D3")
            ),
            ## Item sidepanel
            conditionalPanel(
              condition = paste0(paste0('input[\'', ns('colorBy'), "\'] == \'",nonNumCols,"\'"), collapse = " || "),
              #condition = paste0("input.colorBy == '",nonNumCols,"'", collapse = "||"),
              uiOutput(ns("clusterColorPanelNew"))
            ),
            ## Item sidepanel
            selectInput(
              ns("background"),
              label = "Select Background",
              choices =c("Grey" = "grey", "White" = "white","Minimal" = "minimal", "Plain" =  "plain"),
              selected = "white"
            ),
            
            ## Item sidepanel
            sliderInput(
              ns("dotsize"), 
              "Choose a Dotsize",
              min = 0.01, 
              max = 2, 
              value = 1
            ),
            ## Item sidepanel
            checkboxInput(
              ns("showPlotLegend"), 
              "Show Plot Legends", 
              value = TRUE, 
              width = NULL
            ),
            
            ## Item sidepanel
            bookmarkButton(),
            br(),
            br(),
            ## Item sidepanel
            downloadButton(ns('plotDLall'), "Download Plot Images"),
            br(),
            br(),
            ## Item sidepanel
            conditionalPanel(
              condition = paste0('input[\'', ns('lg10Expr'), "\'] != \'lg10Expr\'"),
              downloadButton(ns("downloadData"), "Download Color Selection")
            )
          ), # End sidepanel
          mainPanel(
              uiOutput(ns("plotsFV"))
          )
  )
}
    
#' scFeatureView Server Functions
#'
#' @noRd 
mod_scFeatureView_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    ################################
    ## Set up data list           ##
    rv <- reactiveValues()
    ## Done                       ##
    ################################
    
    ########################################
    ## Create/update color table                 ##
    
    observeEvent(input$colorBy, {
      rv$dfColorTable <-  createColorTable(
        startUpList = startUpList, 
        colorBy = input$colorBy
      )
    })
    
    ## Done color table                   ##
    ########################################
    
    
    ########################################
    ## Create/update main data table             ##
    toListen <- reactive({
      list(
        input$gene,
        input$splitByColumn,
        input$colorBy,
        input$x_axis,
        input$y_axis
      )
    })
    
    observeEvent(toListen, {
      rv$plotList <-  createDfTemp(
        startUpList = startUpList,
        gene = input$gene,
        splitByColumn = input$splitByColumn,
        colorBy = input$colorBy,
        x_axis = input$x_axis, 
        y_axis = input$y_axis
      )
    })
    
    ## Done main table                    ##
    ########################################
    
    ################################
    ## Load gene list server side ##
    startUpList <- golem::get_golem_options(which = "startUpList")
    allGenes <- startUpList$utilityList$allGenes
    geneDefault <- startUpList$keyList$geneDefault
    
    updateSelectizeInput(
      session, 
      'gene', 
      choices = allGenes, 
      selected = geneDefault, 
      server = TRUE
    )
    
    ## Done loading gene list server side ##
    ########################################
    
    
    
    ########################################
    ## Create Dynamic Color Selectors     ##
    observe({
      startUpList <- golem::get_golem_options(which = "startUpList")
      numCols <- startUpList$utilityList$numCols
      
      if (!(input$colorBy %in% numCols)){
        #dfColorTable <-  createColorTable(startUpList = startUpList, colorBy = input$colorBy)
        dfColorTable <- rv$dfColorTable
        #print(dfColorTable)
        
        nameCol <- names(dfColorTable)[1]
        nameColCol <- "dotColor"
        labelID <- names(dfColorTable)[1]
        
        
        dfColSel <- dfColorTable
        colVec <- as.vector(dfColSel[,nameColCol])
        names(colVec) <- as.vector(dfColSel[,nameCol])
        colVec <- colVec[colVec != ""]
        
        
        dfColSel <- dfColSel[order(dfColSel[,nameCol]),]
        # 
        
        
        output$clusterColorPanelNew = renderUI({
          dfColSel[["label"]] <- paste0(dfColSel[,nameCol], " ", labelID," Color" )
          input_list <- lapply(1:nrow(dfColSel), function(i) {
            # for each dynamically generated input, give a different name
            clusterName <- as.vector(dfColSel[i,nameCol])
            clusterColor <- as.vector(dfColSel[i,nameColCol])
            label <- paste0(as.vector(dfColSel[i,nameCol]), " ",labelID," Color")
            colourpicker::colourInput(inputId = ns(clusterName), label = label, value = clusterColor)
          })
          
          do.call(tagList, input_list)
        })
      }
      
    }) # end observe
    
    ## Done with dynamic color display ##
    #####################################
    
    #####################################
    ## Create dynamic plots            ##
    
    
    
    # Call renderPlot for each one. Plots are only actually generated when they
    # are visible on the web page.
   
    
    observeEvent(reactiveValuesToList(input), {
        startUpList <- golem::get_golem_options( which = "startUpList" )
        numCols <- startUpList$utilityList$numCols
        
        plotList <- createDfTemp(
            startUpList = startUpList,
            gene = input$gene,
            splitByColumn = input$splitByColumn,
            colorBy = input$colorBy,
            x_axis = input$x_axis, 
            y_axis = input$y_axis
        )
        
        dfTemp <- plotList[["dfTemp"]]
        plot_data <- plotList[["plot_data"]]
        plot_data_names <- plotList[["plot_data_names"]]
        maxExpr <- plotList[["maxExpr"]]
        
        req(plot_data)
        
        dimVec <- plotList[["dimVec"]]
        maxX = dimVec[2]
        minX = dimVec[1]
        maxY = dimVec[4]
        minY = dimVec[3]
        
        
        
        #dynamically create the right number of htmlOutput
        output$plotsFV <- renderUI({
          
          
          plot_output_list <- lapply(1:length(plot_data_names), function(i) {
            plotname <- paste("plot", i, sep="")
            plotOutput(ns(plotname))
          })
          
          # Convert the list to a tagList - this is necessary for the list of items
          # to display properly.
          do.call(tagList, plot_output_list)
        })
        
        
        
        
        startUpList <- golem::get_golem_options( which = "startUpList" )
        numCols <- startUpList$utilityList$numCols
        
        plotList <- rv$plotList
        
        dfTemp <- plotList[["dfTemp"]]
        
        if (!(input$colorBy %in% numCols)){
            dfColourSel <- createColorTable(startUpList = startUpList, colorBy = input$colorBy)
            dfColourSel <- dfColourSel[order(dfColourSel[[input$colorBy]]), ]
            cols <- as.vector(dfColourSel[,"dotColor"])
            names(cols) <- as.vector(dfColourSel[,input$colorBy])
            
            inInput <- names(input)[names(input) %in% as.vector(dfColourSel[,input$colorBy])]
            
            
            if (length(inInput) > 0){
              for (k in 1:length(inInput)){
                dfColourSel[dfColourSel[,input$colorBy] == inInput[k], "dotColor"] <- input[[inInput[[k]]]]
              }
              
            }
          
            cols <- as.vector(dfColourSel$dotColor)
            names(cols) <- as.vector(dfColourSel[,input$colorBy])
            
        } else {
            cols <- NULL
        }
        
        
    
        plotResList <- lapply(1:length(plot_data_names), function(i){ 
          featureViewPlot(
              df = plot_data[[i]],
              plot_name = paste0(plot_data_names[i]), 
              colorBy = input$colorBy,
              dotsize = input$dotsize,
              lowColor = input$lowColor, 
              dotcolor = input$dotcolor,
              background = input$background,
              x_axis = input$x_axis,
              y_axis = input$y_axis,
              maxX = maxX,
              minX = minX,
              maxY = maxY,
              minY = minY,
              geneSel = input$gene,
              maxExpr = maxExpr,
              showPlotLegend = input$showPlotLegend,
              colVec = cols
          ) 
          }
        )
        
    
    
    lapply(1:length(plot_data_names), function(i) {
      ## Create data frame for this plot
      
      # Need local so that each item gets its own number. Without it, the value
      # of i in the renderPlot() will be the same across all instances, because
      # of when the expression is evaluated.
      #local({
          my_i <- i
          plotname <- paste("plot", my_i, sep="")
          names(plot_data[[i]])
          
          output[[plotname]] <- renderPlot({
            featureViewPlot(
              df = plot_data[[i]],
              plot_name = paste0(plot_data_names[i]), 
              colorBy = input$colorBy,
              dotsize = input$dotsize,
              lowColor = input$lowColor, 
              dotcolor = input$dotcolor,
              background = input$background,
              x_axis = input$x_axis,
              y_axis = input$y_axis,
              maxX = maxX,
              minX = minX,
              maxY = maxY,
              minY = minY,
              geneSel = input$gene,
              maxExpr = maxExpr,
              showPlotLegend = input$showPlotLegend,
              colVec = cols
            ) 
          })
          
        #})
    }) # end lapply
        
        ## Create downloads
        output$plotDLall <- downloadHandler(
          filename = function() {
            randomString <- function(n = 5000) {
              a <- do.call(paste0, replicate(1, sample(LETTERS, n, TRUE), FALSE))
              paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
            }
            paste0("FeatureViewPlot.",randomString(1),".pdf")
          },
          content = function(file) {
            pdf(file)
            lapply(plotResList, print)
            dev.off()
          }
        )
        
    })
    
    ## Done create dynamic plots       ##
    #####################################
    
   
 
  })
}
    
## To be copied in the UI
# mod_scFeatureView_ui("scFeatureView_1")
    
## To be copied in the server
# mod_scFeatureView_server("scFeatureView_1")
