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
 



###############################################################################
## Main App Server                                                           ##



app_server <- function(input, output, session) {
  
###############################################################################
## Load dfCoord from db                                                      ##
    
    createDfCoord <- reactive({
        oldw <- getOption("warn")
        options(warn = -1)
        startUpList <- golem::get_golem_options(which = "startUpList")
        keyList <- startUpList$keyList
        if (keyList[["dataMode"]] == "SQLite"){
            
            dbDB <- DBI::dbConnect(
                drv = RSQLite::SQLite(),
                dbname=keyList[["dbname"]]
            )
            
        } else {
            
            dbDB <- DBI::dbConnect(
                drv = RMySQL::MySQL(),
                user = keyList[["user"]], 
                password = keyList[["DBpd"]], 
                host = keyList[["host"]], 
                dbname=keyList[["dbname"]]
                
            )
            
        }
        
        #dbDB <- RMySQL::dbConnect(RMySQL::MySQL(), user = user, password = DBpd, host = host, dbname=dbname)
        query <- paste0("SELECT DISTINCT * FROM ", keyList$coordinateTbName)
        dfCoordSel <- DBI::dbGetQuery(dbDB, query)
        
        RMySQL::dbDisconnect(dbDB)
        
        dfCoordSel$row_names <- NULL
        dfCoordSel[["all"]] <- "all"  
        
        
        
        dfCoordSel$cellID <- gsub("[.]", "-", dfCoordSel$cellID)
        dfCoordSel$cellID <- gsub("-", "_", dfCoordSel$cellID)
        
        dfCoordSel
        
    })
    #end_time <- Sys.time()
    #print(paste0("Q S1 DBQ Coordinates: ",end_time - start_time))
    ##                                                                           ##
    ###############################################################################
    
    ###############################################################################
    ## Create Color Table                                                        ##
    
    dfColorTable <- reactive({
        dfDL <- createDfCoord()
        startUpList <- golem::get_golem_options(which = "startUpList")
        dfColOptions <- startUpList$dfColOptions
        #######################################################################
        ## Check if colors are available                                     ##
        colorAnnoFound <- FALSE
        
        if (!is.null(dfColOptions)){
            dfPlotCols <- dfColOptions[dfColOptions$menuName == input$colorBy, c("colOption", "colSel")]
            
            if (nrow(dfPlotCols) > 0){
                checkDLvec <- sort(na.omit(as.vector(unique(dfDL[,input$colorBy]))))
                checkColvec <- sort(na.omit(unique(dfPlotCols$colOption)))
                if (identical(checkDLvec, checkColvec)){
                    dfAddCol <- unique(dfPlotCols)
                    names(dfAddCol) <- gsub("colOption", input$colorBy, names(dfAddCol))
                    names(dfAddCol) <- gsub("colSel", "dotColor", names(dfAddCol))
                    
                    dfDL[,input$colorBy] <- as.character(dfDL[,input$colorBy])
                    dfAddCol[,input$colorBy] <- as.character(dfAddCol[,input$colorBy])
                    
                    
                    dfDL <- dplyr::full_join(
                        dfDL,
                        dfAddCol,
                        by=input$colorBy
                    )
                    
                    
                    dfDL[is.na(dfDL)] <- ""
                    selVec <- c(input$colorBy, "dotColor")
                    dfDL <- unique(dfDL[,selVec])
                    colorAnnoFound <- TRUE
                }
                
            }
            
        } 
        
        
        ## Done                                                              ##
        #######################################################################
        
        
       
        dfDL[["lg10Expr"]] <- "A1"
        
        
        if(!colorAnnoFound) {
            dfDL[["dotColor"]] <- "#000000"
            selVec <- c(input$colorBy, "dotColor")
            dfDL <- unique(dfDL[,selVec])
            dfDL <- dfDL[order(dfDL[,1], decreasing = F),]
            dfDL[,input$colorBy] <- factor(dfDL[,input$colorBy])
                
            dfDL[["dotColor"]] <- scales::hue_pal()(nrow(dfDL))
        }
            
        #}
        
        dfDL <- dfDL[dfDL[,input$colorBy] != "", ]
        dfDL <- dfDL[!is.na(dfDL[,input$colorBy]), ]
        
        dfDL
    })
    
    ##                                                                       ##
    ###########################################################################
    
    #########################################################################
    ## Retrieve Coordinates for this query
    
    ## Done retrieving Coordinates
    #########################################################################
    
    ###########################################################################
    ## Database query for dfExpr                                             ##
    ## create agl315_gene_expr_tb
    #start_time <- Sys.time()
    createDfExprSel <- reactive({
        oldw <- getOption("warn")
        options(warn = -1)
        startUpList <- golem::get_golem_options(which = "startUpList")
        keyList <- startUpList$keyList
        
        if (keyList[["dataMode"]] == "SQLite"){
            
            dbDB <- DBI::dbConnect(
                drv = RSQLite::SQLite(),
                dbname=keyList[["dbname"]]
            )
            
        } else {
            
            dbDB <- DBI::dbConnect(
                drv = RMySQL::MySQL(),
                user = keyList[["user"]], 
                password = keyList[["DBpd"]], 
                host = keyList[["host"]], 
                dbname=keyList[["dbname"]]
                
            )
            
        }
        
        if ( is.null(input$gene) | input$gene == "" ){
            query <- paste0("SELECT * FROM ",keyList$exprTbName," WHERE gene = '",keyList$geneDefault,"'" )  
        } else {
            query <- paste0("SELECT * FROM ",keyList$exprTbName," WHERE gene = '",input$gene,"'" )
        }
        
        #query <- paste0("SELECT DISTINCT * FROM agl315_gene_expr_tb WHERE gene = 'GFAP'" )
        dfExprSel <- DBI::dbGetQuery(dbDB, query)
        DBI::dbDisconnect(dbDB)
        
        names(dfExprSel) <- gsub("condition", "cellID", names(dfExprSel))
        names(dfExprSel) <- gsub("^expr$", "lg10Expr", names(dfExprSel))
        dfExprSel$cellID <- gsub("[.]", "-", dfExprSel$cellID)
        dfExprSel$cellID <- gsub("-", "_", dfExprSel$cellID)
        dfExprSel$cellID <- gsub("^X", "", dfExprSel$cellID)
        dfExprSel
    })
    
    #end_time <- Sys.time()
    #print(paste0("Q S2 agl315_gene_expr_tb: ",end_time - start_time))
    #paste0("SELECT DISTINCT gene, condition, expr FROM agl315_gene_expr_tb WHERE gene = '",input$gene,"'" )
    ## Done db query                                                         ##
    ###########################################################################
    
    ###############################################################################
    ## Create dfTemp                                                             ##       
    createDfTemp <- reactive({
        
        dfTemp <- dplyr::full_join(
            createDfCoord(), 
            createDfExprSel(), 
            by="cellID"
        )
        
        
        #dfTemp2 <- merge(createDfCoord(), createDfExprSel(), by.x = "cellID", by.y="cellID", all=TRUE)
        dfTemp[is.na(dfTemp)] <- 0
        dfTemp <- data.frame(dfTemp, stringsAsFactors = FALSE)
        dfTemp$gene <- as.character(dfTemp$gene)
        
        conditionVec <- sort(unique(dfTemp[,input$splitByColumn]))  
        
        #######################################################################
        ## Check if custom colors are to be used                             ##
        
        dfTemp[["Dcolor"]] <- dfTemp[,input$colorBy]
        
        dfTemp[["dotColor"]] <- "#000000"
        
        #if (input$colorBy == "clusterName"){
        #    dfTemp[["dotColor"]] <- dfTemp[["clusterColor"]]
        #} else if (input$colorBy == "sampleName"){
        #    dfTemp[["dotColor"]] <- dfTemp[["sampleColor"]]
        #}
        
        inInput <- names(input)[names(input) %in% unique(dfTemp[[input$colorBy]])]
        if (length(inInput) > 0){
            for (k in 1:length(inInput)){
                dfTemp[dfTemp[,input$colorBy] == inInput[k], "dotColor"] <- input[[inInput[[k]]]]
            }

        }
        
       
        ## Done                                                              ##
        #######################################################################
        
        pos <- grep(paste0("^", input$x_axis, "$"), names(dfTemp))
        if (length(pos) > 0){
            dfTemp[["x_axis"]] <- dfTemp[,input$x_axis]
        } else {
            dfTemp[["x_axis"]] <- input$x_axis
        }
        
        
        if (!is.numeric(dfTemp$x_axis)){
            dfTemp$x_axis <- factor(dfTemp$x_axis, levels = sort(unique(dfTemp$x_axis)))
        }
        
        
        ## We need to consider cases like Densityplot and Histogram, where input$y_axis is not a column
        if (input$y_axis %in% names(dfTemp)){
            dfTemp[["y_axis"]] <- dfTemp[,input$y_axis]
        } else {
            dfTemp[["y_axis"]] <- input$y_axis
        }
        
        
        #dfTemp <- dfTemp[,selVec]  
        dfTemp <- dfTemp[(dfTemp$x_axis != 0 | dfTemp$y_axis != 0),] 
        #dfTemp
        
        #################
        ## Create plot select
        #plot_select <- reactive({
        df <- dfTemp
        df[["all"]] <- "all"
        plot_select <-  as.vector(unique(df[, input$splitByColumn]))
        #})
        ## Done Creating plot select
        ####################
        
        ####################
        ## Create plot data names
        plot_data_names <- sort(as.vector(unique(dfTemp[, input$splitByColumn])))
        ##
        ####################
        
        ###################
        ## get max expr
        #maxExpr <- reactive({
        #  dfTemp <- createDfTemp()
          
          if (is.numeric(dfTemp$Dcolor)){
            maxExpr <- max(as.vector(dfTemp$Dcolor))
          } else{
            maxExpr <- NULL
          }
      #    return(maxExpr)
      #  })
        
        ##
        ####################
        
        ###################
        ## plot data
        #plot_data <- reactive({
        #  dfTemp <- createDfTemp()
          
          #plot_select <- plot_data_names()
          
          plot_data <- lapply(plot_data_names, function(x) dfTemp[dfTemp[,input$splitByColumn] == x,])
        #})
        ##
        ###################
        
        ######################
        ## min/max
        #determinePlotDims <- reactive({
        #  dfTemp <- createDfTemp()
          
          if (!is.numeric(dfTemp$x_axis)){
            minX <- 0
            maxX <- length(unique(dfTemp$x_axis)) + 1
          } else {
            maxX <- 1.1*max(dfTemp$x_axis, na.rm = T)
            minX <- 1.1*min(dfTemp$x_axis, na.rm = T)
          }
          
          if (!is.numeric(dfTemp$y_axis)){
            minY <- 0
            maxY <- length(unique(dfTemp$y_axis)) + 1
          } else {
            minY <- 1.1*min(dfTemp$y_axis, na.rm = T)
            maxY <- 1.1*max(dfTemp$y_axis, na.rm = T)
          }
          
          
          dimVec <- c(minX, maxX, minY, maxY)
          dimVec
        #})
        
        ##
        #######################
        
        returnList <- list(
            "dfTemp" = dfTemp,
            "plot_select" = plot_select,
            "plot_data_names" = plot_data_names,
            "plot_data" = plot_data,
            "maxExpr" = maxExpr,
            "dimVec" = dimVec
        )
        
        return(returnList)
        
    })
    ##                                                                           ##
    ###############################################################################
    
    
    
    
###############################################################################
##                                                                           ##
    
    
    startUpList <- golem::get_golem_options(which = "startUpList")
    allGenes <- startUpList$utilityList$allGenes
    geneDefault <- startUpList$keyList$geneDefault
    
    observe({
        updateSelectizeInput(session, 'gene', choices = allGenes, server = TRUE, selected=geneDefault) 
    })
    
    onRestored(function(state) {
        updateSelectizeInput(session, "gene", selected=state$input$gene, allGenes, server=TRUE)
    })
    
    ###################################################################
    ## Add dropdownselection                                         ##
    
    observe({
        startUpList <- golem::get_golem_options(which = "startUpList")
        dropDownList <- startUpList$utilityList$dropDownList
        
        if (length(dropDownList) > 0){
            
            output$dropDownPanel = renderUI({
                #dfColSel[["label"]] <- paste0(dfColSel[,nameCol], " ", labelID," Color" )
                dropdown_list <- lapply(1:length(dropDownList), function(i) {
                    # for each dynamically generated input, give a different name
                    displayID <- names(dropDownList)[i]
                    displayName <- as.vector(dropDownList[[displayID]][["displayName"]])
                    selOptions <- dropDownList[[displayID]][["selOptions"]]
                    selDisplayOptions <- as.vector(dropDownList[[displayID]][["selDisplayOptions"]])
                    
                    names(selOptions) <- selDisplayOptions
                    
                    default <- as.vector(dropDownList[[displayID]][["default"]])
                    
                    selectInput(
                        displayID, 
                        label = displayName,
                        choices = selOptions,
                        selected = default
                    )
                    
                })
                
                do.call(tagList, dropdown_list)
            })
        }
        
    })
    
    ## Done dropdown                                                 ##
    ###################################################################
    
    observe({
        startUpList <- golem::get_golem_options(which = "startUpList")
        numCols <- startUpList$utilityList$numCols
        
        if (!(input$colorBy %in% numCols)){
            dfColorTable <-  dfColorTable()
            
            nameCol <- names(dfColorTable)[1]
            nameColCol <- "dotColor"
            labelID <- names(dfColorTable)[1]
            
           
            dfColSel <- dfColorTable
            colVec <- as.vector(dfColSel[,nameColCol])
            names(colVec) <- as.vector(dfColSel[,nameCol])
            colVec <- colVec[colVec != ""]
            dfColSel <- dfColSel[order(dfColSel[,nameCol]),]
                
            output$clusterColorPanel = renderUI({
                dfColSel[["label"]] <- paste0(dfColSel[,nameCol], " ", labelID," Color" )
                    input_list <- lapply(1:nrow(dfColSel), function(i) {
                        # for each dynamically generated input, give a different name
                        clusterName <- as.vector(dfColSel[i,nameCol])
                        clusterColor <- as.vector(dfColSel[i,nameColCol])
                        label <- paste0(as.vector(dfColSel[i,nameCol]), " ",labelID," Color")
                        colourInput(inputId = clusterName, label = label, value = clusterColor)
                    })
                
                do.call(tagList, input_list)
            })
        }
        
    })
    
    
###############################################################################
##                                                                           ##
    
    observeEvent(reactiveValuesToList(input), {
        plotList <- createDfTemp()
        
        plot_data <- plotList[["plot_data"]]
        plot_data_names <- plotList[["plot_data_names"]]
        maxExpr <- plotList[["maxExpr"]]
        
        req(plot_data)
        
        dimVec <- plotList[["dimVec"]]
        maxX = dimVec[2]
        minX = dimVec[1]
        maxY = dimVec[4]
        minY = dimVec[3]
        
###############################################################################
##        
        output$multi_plot_ui <- renderUI({
            
            lapply(seq_along(plot_data),
                   function(n) {
                       return(plot_prep_ui(paste0("n", n)))
                   })
        })
        
##
###############################################################################
      
        
        
        lapply(seq_along(plot_data),
               function(i){
                   callModule(plot_prep_server,
                              paste0("n", i),
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
                              showPlotLegend = input$showPlotLegend
                   )
               }
        )
        
        
        ## Make plot list for download
        res <- lapply(seq_along(plot_data),
                      function(i){
                          callModule(plot_prep_server_dl,
                                     paste0("n", i),
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
                                     showPlotLegend = input$showPlotLegend
                          )
                      }
        )
         
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
                lapply(res, print)
                dev.off()
            }
        )
        
        output$downloadData <- downloadHandler(
            
            
            filename = function() {
                paste(input$colorBy, ".color.selection.csv", sep = "")
            },
            content = function(file) {
                
                write.csv(dfColorTable(), file, row.names = FALSE)
            }
        )
        
        output$selected_var <- renderText({ 
            paste0("You have selected this: ", length(res), ". Class:", class(res[[1]]))
        })
        
    })
##    
###############################################################################    
    
    
    
    
    
}
