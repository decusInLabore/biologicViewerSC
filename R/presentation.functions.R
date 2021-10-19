###############################################################################
## Function createDfTempLocal                                                ##

#' @title createInputListFromLink
#'
#'
#' @param figureLink A figure link from the app
#' @import dplyr
#' @return input
#' @export

createInputListFromLink <- function(
    figureLink = NULL
){
    ## Create input list from link ##
    library(dplyr)
    figureLink <- gsub("%22", "", figureLink)
    figureLink <- figureLink %>% URLdecode() 
    inputList <- unlist(strsplit(figureLink, "\\?_inputs_&"))[2]
    inputList <- strsplit(inputList, "&")
    listNames <- sapply(
        strsplit(unlist(inputList), "="),
        function(x) unlist(x, 1)[1]
    )
    listValues <-  sapply(
        strsplit(unlist(inputList), "="),
        function(x) unlist(x, 1)[2]
    )
    listValues <- gsub("%22", "", listValues)
    
    input <- split(listValues, listNames)
    
    pos <- grep("showPlotLegend", names(input))
    
    if (length(pos) > 0){
        input$showPlotLegend <- as.logical(input$showPlotLegend)
    }
    
    return(input)
}

##                                                                           ##
###############################################################################

###############################################################################
## Function createDfTempLocal                                                ##

#' @title createDfExprSelReport
#'
#'
#' @param dfID An ID data frame
#' @param input An input object
#' @param dfColOptions A color options table 

#' @import DBI
#' @import RSQLite
#' @import RMySQL
#' @return dfCoordSel
#' @export

createDfCoordLocal <- function(dfID, dfColOptions = NULL, input = NULL){
    oldw <- getOption("warn")
    options(warn = -1)
    
    if (dfID[["dataMode"]] == "SQLite"){
        
        dbDB <- DBI::dbConnect(
            drv = RSQLite::SQLite(),
            dbname=dfID[["db"]]
        )
        
    } else {
        
        dbDB <- DBI::dbConnect(
            drv = RMySQL::MySQL(),
            user = dfID[["id"]], 
            password = dfID[["id2"]], 
            host = dfID[["url"]], 
            dbname=dfID[["db"]]
            
        )
        
    }
    
    #dbDB <- RMySQL::dbConnect(RMySQL::MySQL(), user = user, password = DBpd, host = host, dbname=dbname)
    query <- paste0("SELECT DISTINCT * FROM ", dfID$coordTb)
    dfCoordSel <- DBI::dbGetQuery(dbDB, query)
    
    RMySQL::dbDisconnect(dbDB)
    
    dfCoordSel$row_names <- NULL
    dfCoordSel[["all"]] <- "all"  
    
    
    
    dfCoordSel$cellID <- gsub("[.]", "-", dfCoordSel$cellID)
    dfCoordSel$cellID <- gsub("-", "_", dfCoordSel$cellID)
    
    ## Add dotcolor for non-numeric ##
    colorAnnoFound <- FALSE
    if (!is.null(dfColOptions)){
        dfPlotCols <- dfColOptions[dfColOptions$menuName == input$colorBy, c("colOption", "colSel")]
        
        if (nrow(dfPlotCols) > 0){
            checkDLvec <- sort(na.omit(as.vector(unique(dfCoordSel[,input$colorBy]))))
            checkColvec <- sort(na.omit(unique(dfPlotCols$colOption)))
            if (identical(checkDLvec, checkColvec)){
                dfAddCol <- unique(dfPlotCols)
                names(dfAddCol) <- gsub("colOption", input$colorBy, names(dfAddCol))
                names(dfAddCol) <- gsub("colSel", "dotColor", names(dfAddCol))
                
                dfCoordSel[,input$colorBy] <- as.character(dfCoordSel[,input$colorBy])
                dfAddCol[,input$colorBy] <- as.character(dfAddCol[,input$colorBy])
                
                
                dfCoordSel <- dplyr::full_join(
                    dfCoordSel,
                    dfAddCol,
                    by=input$colorBy
                )
                
                
                dfCoordSel[is.na(dfCoordSel)] <- ""
                colorAnnoFound <- TRUE
            }
            
        }
    }
    
    
    
    ## Done                                                              ##
    #######################################################################
    
    
    
    dfCoordSel[["lg10Expr"]] <- "A1"
    
    
    if(!colorAnnoFound) {
        dfCoordSel[["dotColor"]] <- "#000000"
        #selVec <- c(input$colorBy, "dotColor")
        #dfCoordSel <- unique(dfCoordSel[,selVec])
        dfCoordSel <- dfCoordSel[order(dfCoordSel[,1], decreasing = F),]
        dfCoordSel[,input$colorBy] <- factor(dfCoordSel[,input$colorBy])
        
        dfCoordSel[["dotColor"]] <- scales::hue_pal()(nrow(dfCoordSel))
    }
    
    dfCoordSel$lg10Expr <- NULL
    
    return(dfCoordSel)
    
}
#end_time <- Sys.time()
#print(paste0("Q S1 DBQ Coordinates: ",end_time - start_time))
##                                                                           ##
###############################################################################


###############################################################################
## Function createDfTempLocal                                                ##

#' @title dfColorTableReport
#'
#'
#' @param dfMeta A meta data table
#' @param input An input object
#' @param dfColOptions A color options table

#' @import dplyr
#' @return dfDL
#' @export


dfColorTableReport <- function(dfMeta, input, dfColOptions = NULL){
    dfDL <- dfMeta
    dfDL$dotColor <- NULL
    #startUpList <- golem::get_golem_options(which = "startUpList")
    #dfColOptions <- startUpList$dfColOptions
    
    #dfColOptions 
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
    
    return(dfDL)
}

##                                                                       ##
###########################################################################

###############################################################################
## Function createDfTempLocal                                                ##

#' @title createDfExprSelReport
#'
#'
#' @param dfID An ID data frame
#' @param input An input object

#' @import DBI
#' @import RSQLite
#' @import RMySQL
#' @return dfExprSel
#' @export


createDfExprSelReport <- function(dfID, input){
    oldw <- getOption("warn")
    options(warn = -1)
    
    if (dfID[["dataMode"]] == "SQLite"){
        
        dbDB <- DBI::dbConnect(
            drv = RSQLite::SQLite(),
            dbname=dfID[["exprTb"]]
        )
        
    } else {
        
        dbDB <- DBI::dbConnect(
            drv = RMySQL::MySQL(),
            user = dfID[["id"]], 
            password = dfID[["id2"]], 
            host = dfID[["url"]], 
            dbname=dfID[["db"]]
            
        )
        
    }
    
    if ( is.null(input$gene) | input$gene == "" ){
        query <- paste0("SELECT * FROM ",dfID$exprTb," WHERE gene = '",keyList$geneDefault,"'" )  
    } else {
        query <- paste0("SELECT * FROM ",dfID$exprTb," WHERE gene = '",input$gene,"'" )
    }
    
    #query <- paste0("SELECT DISTINCT * FROM agl315_gene_expr_tb WHERE gene = 'GFAP'" )
    dfExprSel <- DBI::dbGetQuery(dbDB, query)
    DBI::dbDisconnect(dbDB)
    
    names(dfExprSel) <- gsub("condition", "cellID", names(dfExprSel))
    names(dfExprSel) <- gsub("^expr$", "lg10Expr", names(dfExprSel))
    dfExprSel$cellID <- gsub("[.]", "-", dfExprSel$cellID)
    dfExprSel$cellID <- gsub("-", "_", dfExprSel$cellID)
    dfExprSel$cellID <- gsub("^X", "", dfExprSel$cellID)
    return(dfExprSel)
}

#end_time <- Sys.time()
#print(paste0("Q S2 agl315_gene_expr_tb: ",end_time - start_time))
#paste0("SELECT DISTINCT gene, condition, expr FROM agl315_gene_expr_tb WHERE gene = '",input$gene,"'" )
## Done db query                                                         ##
###########################################################################

###############################################################################
## Function createDfTempLocal                                                ##

#' @title createDfTempLocal
#'
#'
#' @param dfMeta A Metadata table
#' @param dfExpr An expression value table
#' @param dfExpr A color table
#' @import dplyr
#' @export

createDfTempLocal <- function(dfMeta, dfExpr, input, dfColorTable = NULL){
    
    dfTemp <- dplyr::full_join(
        dfMeta, 
        dfExpr, 
        by="cellID"
    )
    
    
    #dfTemp2 <- merge(createDfCoord(), createDfExprSel(), by.x = "cellID", by.y="cellID", all=TRUE)
    dfTemp[is.na(dfTemp)] <- 0
    dfTemp <- data.frame(dfTemp, stringsAsFactors = FALSE)
    dfTemp$gene <- as.character(dfTemp$gene)
    
    #conditionVec <- sort(unique(dfTemp[,input$splitByColumn]))  
    #conditionVec <- conditionVec[conditionVec != 0]
    
    #######################################################################
    ## Check if custom colors are to be used                             ##
    
    dfTemp[["Dcolor"]] <- dfTemp[,input$colorBy]
    
    #dfTemp[["dotColor"]] <- "#000000"
    
    #if (input$colorBy == "clusterName"){
    #    dfTemp[["dotColor"]] <- dfTemp[["clusterColor"]]
    #} else if (input$colorBy == "sampleName"){
    #    dfTemp[["dotColor"]] <- dfTemp[["sampleColor"]]
    #}
    
    # inInput <- names(input)[names(input) %in% unique(dfTemp[[input$colorBy]])]
    # if (length(inInput) > 0){
    #     for (k in 1:length(inInput)){
    #         dfTemp[dfTemp[,input$colorBy] == inInput[k], "dotColor"] <- input[[inInput[[k]]]]
    #     }
    # 
    # }
    
    baseCols <- names(dfTemp)[!(names(dfTemp) %in% c("row_names", "cellID"))]
    numCols <- names(dfMeta)[unlist(lapply(dfMeta, is.numeric))]
    numCols <- numCols[!(numCols %in% c("row_names", "cellID"))]
    
    
    numCols <- c("lg10Expr", numCols)
    nonNumCols <- baseCols[!(baseCols %in% numCols)]
    
    
    if (!(input$colorBy %in% numCols)){
        dfTemp$dotColor <- NULL
        dfColorTable$lg10Expr <- NULL
        
        dfTemp <- dplyr::full_join(
            dfTemp,
            dfColorTable,
            by=input$colorBy
        )
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
    #  dfTemp <- createDfTempLocal()
    
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
    #  dfTemp <- createDfTempLocal()
    
    #plot_select <- plot_data_names()
    
    plot_data <- lapply(plot_data_names, function(x) dfTemp[dfTemp[,input$splitByColumn] == x,])
    #})
    ##
    ###################
    
    ######################
    ## min/max
    #determinePlotDims <- reactive({
    #  dfTemp <- createDfTempLocal()
    
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
    #dimVec
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
    
}

## Enf of function                                                           ##
###############################################################################


###############################################################################
## Helper function

#' @title plot_prep_report
#'
#'
#' @import dplyr
#' @import ggplot2
#' @export
#' 
plot_prep_report <- function(
    #input,
    #output,
    #session, 
    df,
    plot_name,
    colorBy = "lg10Expr",
    dotsize = 1,
    lowColor = "grey", 
    dotcolor = "darkblue",
    x_axis = "UMAP_1",
    y_axis = "UMAP_2",
    background = "grey",
    maxX = NULL,
    minX = NULL,
    maxY = NULL,
    minY = NULL,
    geneSel = NULL,
    maxExpr = NULL,
    showPlotLegend = FALSE
) {
    if (is.null(maxX)){
        maxX <- 1.1*max(df$x_axis, na.rm = T)  
    } 
    
    if (is.null(maxY)){
        maxY <- 1.1*max(df$y_axis, na.rm = T)  
    }
    
    if (is.null(minX)){
        minX <- 1.1*min(df$x_axis, na.rm = T)  
    } 
    
    if (is.null(minY)){
        minY <- 1.1*min(df$y_axis, na.rm = T)  
    }
    
    
    
    nCellsTotal <- nrow(df)
    nExpr <- nrow(df[df$gene != 0,])
    percExpr <- 100*round(nrow(df[df$gene != 0,])/nCellsTotal, 3)
    qGene <- unique(na.omit(df$gene))
    qGene <- qGene[qGene != 0]
    
    baseCols <- names(df)[!(names(df) %in% c("row_names", "cellID"))]
    numCols <- names(df)[unlist(lapply(df, is.numeric))]
    numCols <- numCols[!(numCols %in% c("row_names", "cellID"))]
    
    
    numCols <- c("lg10Expr", numCols)
    nonNumCols <- baseCols[!(baseCols %in% numCols)]
    
    df[["Dcolor"]] <- "A1"
    df[["Dcolor"]] <- df[,colorBy]
    
    if (colorBy %in% nonNumCols ){
        df$Dcolor[df$Dcolor == ""] <- "Rest"
        df$Dcolor <- factor(df$Dcolor)
    } else if( is.numeric( df$Dcolor ) ) {
        minExpr <- floor(min(df$Dcolor, na.rm = T))
        
        if (is.null(maxExpr)){
            maxExpr <- ceiling(max(df$Dcolor, na.rm = T))   
            if (maxExpr == 1){
                ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)
                maxExpr <- ceiling_dec(max(df$Dcolor, na.rm = T),2)   
            }
        }
        
    } else {
        df$Dcolor[df$Dcolor == ""] <- "Rest"
        df$Dcolor <- factor(df$Dcolor)
    }     
    
    
    
    ###########################################################################
    ## Determine plot logic                                                  ##
    if (is.numeric(df$x_axis)){
        #######################################################################
        ## Decide on factorial display logic                                 ##
        if (df$y_axis[1] == "Densityplot"){
            plotLogic <- "density"
            p <- ggplot2::ggplot(
                data = df, ggplot2::aes(x=x_axis, y=..density.., color=Dcolor,fill=Dcolor)
            ) + ggplot2::geom_density(alpha=0.3, position="stack") 
        } else if (df$y_axis[1] == "Histogram"){
            plotLogic <- "histogram"
            Nbin <- ceiling(length(df$x_axis)/5)
            p <- ggplot2::ggplot(
                data = df, ggplot2::aes(x=x_axis, color=Dcolor,fill=Dcolor)
            ) + ggplot2::geom_histogram(alpha=0.3, position="stack", bins = Nbin)
        } else if (df$y_axis[1] == "Ridgeplot"){
            plotLogic <- "ridgeplot"
            p <- ggplot2::ggplot(data = df, ggplot2::aes_string(x = "x_axis", y = colorBy, fill=colorBy, color=colorBy)
            ) + ggridges::geom_density_ridges()
            
        } else {
            plotLogic <- "point"
            p <- ggplot2::ggplot(
                data = df, ggplot2::aes(x_axis, y_axis, color=Dcolor)
            ) + ggplot2::geom_point(
                shape = 16,
                size = as.numeric(dotsize)
            ) 
        }
        ## Done deciding factorial display logic
        #########################################################################  
    } else {
        
        plotLogic <- "violin"
        p <- ggplot2::ggplot(
            data = df, ggplot2::aes(x_axis, y_axis, color=Dcolor)
        ) + ggplot2::geom_violin(trim=FALSE, fill="#E8E8E8"
        )+ ggplot2::geom_jitter(height = 0, size = as.numeric(dotsize))
    }
    ## Done plot logic                                                       ##
    ###########################################################################
    
    
    p <- p + ggplot2::xlab(x_axis) + ggplot2::ylab(y_axis)
    
    if (colorBy %in% nonNumCols ){
        dfCol <- unique(df[,c(colorBy, "dotColor")])
        colVec <- dfCol$dotColor
        names(colVec) <- as.character(dfCol[,colorBy])
        colVec <- colVec[colVec != ""]
        
        
        p <- p + ggplot2::scale_colour_manual(colorBy ,values = colVec
        ) + ggplot2::guides(col = ggplot2::guide_legend(override.aes = list(shape = 16, size = 5))
        )
        
        if (plotLogic %in% c("ridgeplot","density", "histogram")){
            p <- p + ggplot2::scale_fill_manual(colorBy ,values = colVec
            ) + ggplot2::guides(
                col = ggplot2::guide_legend(
                    override.aes = list(shape = 16, size = 5)
                )
            )
        }
        
    } else if (is.numeric( df$Dcolor )){
        if (minExpr < 0){
            p <- p + ggplot2::scale_color_gradient2("Expr",low= lowColor, mid = "white", high= dotcolor, midpoint = 0, limits=c(minExpr,maxExpr)
            )
            
        } else {
            p <- p + ggplot2::scale_color_gradient("Expr",low= lowColor, high= dotcolor, limits=c(minExpr,maxExpr)
            )
        }
        
    } 
    
    
    if (background == "white"){
        p <- p + ggplot2::theme_bw()
    } else if (background == "minimal"){
        p <- p + ggplot2::theme_minimal()
    } else if (background == "plain"){
        p <- p + ggplot2::theme_void()
    } else {
        p <- p + ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = "lightgrey")
        )
    }
    
    p <- p + ggplot2::theme(
        axis.text.y   = ggplot2::element_text(size=8),
        axis.text.x   = ggplot2::element_text(size=8),
        axis.title.y  = ggplot2::element_text(size=8),
        axis.title.x  = ggplot2::element_text(size=8),
        axis.line = ggplot2::element_line(colour = "black"),
        panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12)
    ) 
    
    if (is.numeric(df$x_axis)){
        p <- p + ggplot2::xlim(minX, maxX) 
    }
    
    if (is.numeric(df$y_axis)){
        p <- p + ggplot2::ylim(minY, maxY) 
    }
    
    if (colorBy == "lg10Expr" | x_axis == "lg10Expr" | y_axis == "lg10Expr") {
        titleString <- paste0("Sample: ", plot_name, " ", nExpr, "/", nCellsTotal, " cells (",percExpr,"%) express ", geneSel)
    } else {
        titleString <-paste0("Sample: ", plot_name)
    }
    
    p <- p + ggplot2::ggtitle(titleString) 
    #+ ggtitle(paste0("Gene ", input$gene, " in sample ", conditionVec[i], " (E:",cellsExpressingGene[i],"/NE:",cellsNotExpressingGene[i], ", ",percE[i],"%)")) + scale_size_continuous(limits = c(0, maxExpr)
    #) #+ xlim(minX, maxX) + ylim(minY, maxY)
    
    posX <- grep("UMAP", x_axis)
    posY <- grep("UMAP", y_axis)
    if ( (length(posX) == 1) & (length(posY) == 1)){
        p <-  p + ggplot2::coord_fixed()
    }
    
    if (!showPlotLegend){
        p <- p + ggplot2::theme(legend.position = "none")
    } 
    
    
    return(p)
    
}
## Done
###############################################################################


###############################################################################
## Function createAppPlotList                                                ##

#' @title createAppPlotList
#'
#'
#' @param input An input list
#' @param dfID A project parameter data frame
#' @param dfColOptions A color options data frame
#' @return gg plot list
#' @export

createAppPlotList <- function(
    input = "inputList",
    dfID = NULL,
    dfColOptions = NULL
    
){
    #source("presentation.functions.R")
    ## Load app parameters ##
    dfID <- dfID[1,]
    
    
    dfMeta <- createDfCoordLocal(dfID = dfID, dfColOptions = dfColOptions, input = input)
    
    dfExpr <- createDfExprSelReport(dfID = dfID, input = input)
    
    dfColorTable <- dfColorTableReport(dfMeta = dfMeta, input = input, dfColOptions = dfColOptions)
    #dfColorTable$lg10Expr <- NULL
    
    resList <- createDfTempLocal(
        dfMeta = dfMeta, 
        dfExpr = dfExpr,
        input = input,
        dfColorTable = dfColorTable
    )
    
    dfTemp <- resList$dfTemp
    plot_data <- resList$plot_data
    plot_data_names <- resList$plot_data_names
    dimVec <- resList$dimVec
    maxX = dimVec[2]
    minX = dimVec[1]
    maxY = dimVec[4]
    minY = dimVec[3]
    maxExpr <- resList$maxExpr
    
    plotList <- lapply(seq_along(plot_data),
                       function(i){
                           #callModule(
                           plot_prep_report(
                               #paste0("n", i),
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
    
    
    return(plotList)
    
}

## End function                                                              ##
###############################################################################


###############################################################################
## Function createPowerpointPresentation                                     ##

#' @title createPowerpointPresentation
#'
#'
#' @param slideList A slide list
#' @param outputfile Specifies a powerpoint output file. Default is singleCell.powerpoint.presentation.pptx
#' @import officer
#' @import cowplot
#' @export

createPowerpointPresentation <- function(
    slideList, 
    outputfile="singleCell.powerpoint.presentation.pptx"
){
    
    doc <- officer::read_pptx()
    ## Check template presentation ##
    # layout_summary(doc)
    # layout_properties (doc)
    
    for (i in 1:length(slideList)){
        appParameterFile <- slideList[[i]]$appParameterFile
        pathToAppColorParameters <- slideList[[i]]$pathToAppColorParameters
        figureLink <- slideList[[i]]$figureLink
        
        input <- createInputListFromLink(figureLink = figureLink)
        
        dfID <- read.delim(
            appParameterFile,
            stringsAsFactors = F, 
            sep = "\t"
        )
        
        dfColOptions <- read.delim(
            pathToAppColorParameters,
            header = T,
            sep = "\t"
        ) 
        
        tempPlotList <- createAppPlotList(
            input = input,
            dfID = dfID,
            dfColOptions = dfColOptions
        )
        
        plotList <- list()
        
        if (length(tempPlotList) > 1){
            #plotList[[tag]] <- gridExtra::grid.arrange(grobs = tempPlotList, ncol=2)
            plot <- cowplot::plot_grid(plotlist = tempPlotList, ncol=2)
        } else {
            plotList <- tempPlotList
            plot <- plotList[[1]]
        }
        
        
        
        
        #doc <- add_slide(doc)
        doc <- officer::add_slide(
            doc, layout = "Title and Content", 
            master = "Office Theme"
        )
        
        doc <- officer::ph_with(
            x = doc, 
            value = names(slideList)[i], 
            type = "sldNum", 
            location = officer::ph_location_type(
                type = "sldNum"
            )
        )
        
        
        doc <- officer::ph_with(
            x = doc, value = plot, 
            location = officer::ph_location_type(
            type = "body") 
        )
        
        doc <- officer::ph_with(
            x = doc, 
            slideList[[i]]$Slidetitle, 
            location = officer::ph_location_type(
                type = "title"
            ) 
        )
        doc <- officer::ph_with(
            doc, "biologicViewerSC Powerpoint Output", 
            location = officer::ph_location_type(type = "ftr")
        )
    }
    
    print(doc, target = outputfile)
    print(paste0("Presentation generated and saved as ", outputfile))
}
## Function end                                                              ##
###############################################################################
