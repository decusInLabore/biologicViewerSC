###############################################################################
## Load dfCoord from db                                                      ##
#' @title addDf2seuratMetaData
#' 
#' @description Plot function
#' 
#' @param startUpList
#' 
#' @import DBI RSQLite RMySQL
#'
#' @return The return value, if any, from executing the function.
#' 

createDfCoord <- function(
  startUpList
){
  oldw <- getOption("warn")
  options(warn = -1)
  #startUpList <- golem::get_golem_options(which = "startUpList")
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
  
  return(dfCoordSel) 
  
}
##                                                                           ##
###############################################################################

###############################################################################
## Create Color Table                                                        ##
#' @title createColorTable
#' 
#' @description Create color table
#' 
#' @param startUpList
#' @param colorBy
#' 
#' @import dplyr
#'
#' @return The return value, if any, from executing the function.

createColorTable <- function(
  startUpList,
  colorBy = "lg10Expr"
){
  dfDL <- createDfCoord(startUpList = startUpList)
  #startUpList <- golem::get_golem_options(which = "startUpList")
  dfColOptions <- startUpList$dfColOptions
  #######################################################################
  ## Check if colors are available                                     ##
  colorAnnoFound <- FALSE
  
  if (!is.null(dfColOptions)){
    dfPlotCols <- dfColOptions[dfColOptions$menuName == colorBy, c("colOption", "colSel")]
    
    if (nrow(dfPlotCols) > 0){
      checkDLvec <- sort(na.omit(as.vector(unique(dfDL[,colorBy]))))
      checkColvec <- sort(na.omit(unique(dfPlotCols$colOption)))
      if (identical(checkDLvec, checkColvec)){
        dfAddCol <- unique(dfPlotCols)
        names(dfAddCol) <- gsub("colOption", colorBy, names(dfAddCol))
        names(dfAddCol) <- gsub("colSel", "dotColor", names(dfAddCol))
        
        dfDL[,colorBy] <- as.character(dfDL[,colorBy])
        dfAddCol[,colorBy] <- as.character(dfAddCol[,colorBy])
        
        
        dfDL <- dplyr::full_join(
          dfDL,
          dfAddCol,
          by=colorBy
        )
        
        
        dfDL[is.na(dfDL)] <- ""
        selVec <- c(colorBy, "dotColor")
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
    selVec <- c(colorBy, "dotColor")
    dfDL <- unique(dfDL[,selVec])
    dfDL <- dfDL[order(dfDL[,1], decreasing = F),]
    dfDL[,colorBy] <- factor(dfDL[,colorBy])
    
    dfDL[["dotColor"]] <- scales::hue_pal()(nrow(dfDL))
  }
  
  #}
  
  dfDL <- dfDL[dfDL[,colorBy] != "", ]
  dfDL <- dfDL[!is.na(dfDL[,colorBy]), ]
  
  return(dfDL)
}

##                                                                       ##
###########################################################################

###########################################################################
## Database query for dfExpr                                             ##
#' @title createDfExprSel
#' 
#' @description Create Expression data frame
#' 
#' @param startUpList
#' @param gene
#' 
#' @import dplyr DBI RSQLite RMySQL
#'
#' @return The return value, if any, from executing the function.
createDfExprSel <- function(
  startUpList,
  gene
  
){
  oldw <- getOption("warn")
  options(warn = -1)
  # startUpList <- golem::get_golem_options(which = "startUpList")
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
  
  if ( is.null(gene) | gene == "" ){
    query <- paste0("SELECT * FROM ",keyList$exprTbName," WHERE gene = '",keyList$geneDefault,"'" )  
  } else {
    query <- paste0("SELECT * FROM ",keyList$exprTbName," WHERE gene = '",gene,"'" )
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

## Done db query                                                         ##
###########################################################################

###############################################################################
## Create dfTemp                                                             ##  
#' @title createDfTemp
#' 
#' @description Create Expression data frame
#' 
#' @param startUpList
#' @param gene
#' @param splitByColumn,
#' @param colorBy,
#' @param x_axis, 
#' @param y_axis
#' 
#' @import dplyr
#'
#' @return The return value, if any, from executing the function.      
createDfTemp <- function(
  startUpList,
  gene,
  splitByColumn,
  colorBy,
  x_axis, 
  y_axis
){
  
  dfTemp <- dplyr::full_join(
    createDfCoord(startUpList = startUpList), 
    createDfExprSel(startUpList = startUpList, gene = gene), 
    by="cellID"
  )
  
  
  #dfTemp2 <- merge(createDfCoord(), createDfExprSel(), by.x = "cellID", by.y="cellID", all=TRUE)
  dfTemp[is.na(dfTemp)] <- 0
  dfTemp <- data.frame(dfTemp, stringsAsFactors = FALSE)
  dfTemp$gene <- as.character(dfTemp$gene)
  
  conditionVec <- unique(dfTemp[,splitByColumn])
  
  #######################################################################
  ## Check if custom colors are to be used                             ##
  
  dfTemp[["Dcolor"]] <- dfTemp[,colorBy]
  
  dfTemp[["dotColor"]] <- "#000000"
  
  #if (colorBy == "clusterName"){
  #    dfTemp[["dotColor"]] <- dfTemp[["clusterColor"]]
  #} else if (colorBy == "sampleName"){
  #    dfTemp[["dotColor"]] <- dfTemp[["sampleColor"]]
  #}
  
  
  
  
  ## Done                                                              ##
  #######################################################################
  
  pos <- grep(paste0("^", x_axis, "$"), names(dfTemp))
  if (length(pos) > 0){
    dfTemp[["x_axis"]] <- dfTemp[,x_axis]
  } else {
    dfTemp[["x_axis"]] <- x_axis
  }
  
  
  if (!is.numeric(dfTemp$x_axis)){
    dfTemp$x_axis <- factor(dfTemp$x_axis, levels = sort(unique(dfTemp$x_axis)))
  }
  
  
  ## We need to consider cases like Densityplot and Histogram, where y_axis is not a column
  if (y_axis %in% names(dfTemp)){
    dfTemp[["y_axis"]] <- dfTemp[,y_axis]
  } else {
    dfTemp[["y_axis"]] <- y_axis
  }
  
  
  #dfTemp <- dfTemp[,selVec]  
  dfTemp <- dfTemp[(dfTemp$x_axis != 0 | dfTemp$y_axis != 0),] 
  #dfTemp
  
  #################
  ## Create plot select
  #plot_select <- reactive({
  df <- dfTemp
  df[["all"]] <- "all"
  plot_select <-  as.vector(unique(df[, splitByColumn]))
  #})
  ## Done Creating plot select
  ####################
  
  ####################
  ## Create plot data names
  plot_data_names <- sort(as.vector(unique(dfTemp[, splitByColumn])))
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
  
  plot_data <- lapply(plot_data_names, function(x) dfTemp[dfTemp[,splitByColumn] == x,])
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
  
}
##                                                                           ##
###############################################################################


######################################################################$#########
## FeatureViewPlot Function                                                   ## 

#' @title featureViewPlot
#' 
#' @description Create Expression data frame
#' 
#' @param df
#' @param plot_name
#' @param colorBy,
#' @param dotsize,
#' @param lowColor, 
#' @param dotcolor
#' @param x_axis
#' @param y_axis
#' @param background
#' @param maxX
#' @param minX
#' @param maxY
#' @param minY
#' @param geneSel
#' @param maxExpr
#' @param showPlotLegend
#' @param colVec
#' 
#' @import dplyr
#'
#' @return The return value, if any, from executing the function.

featureViewPlot <- function(
  df,
  plot_name,
  colorBy = "lg10Expr",
  dotsize = "dotsize",
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
  showPlotLegend = FALSE,
  colVec = NULL
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
  
  #plotInput <- reactive({
  startUpList <- golem::get_golem_options(which = "startUpList")
  nonNumCols <- startUpList$utilityList$nonNumCols
  
  if (colorBy %in% startUpList$nonNumCols ){
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
    if (df$y_axis[1] == "Barchart" | df$y_axis[1] == "Histogram"){
      plotLogic <- "barchart"
      p <- ggplot2::ggplot(
        data = df, ggplot2::aes(x= x_axis, fill=Dcolor)) + ggplot2::geom_bar(color="black")  
      if (showPlotLegend){
        p <- p + geom_text(stat='count', aes(label=..count..), position = position_stack(vjust = 0.5))
      }
    } else {
      plotLogic <- "violin"
      
      ## Add noise to y_axis to mimit Seurat displays as per
      # https://github.com/satijalab/seurat/issues/3322
      noise <- rnorm(n = length(x = df[, "y_axis"])) / 100000
      df$y_axis <- df$y_axis + noise
      
      p <- ggplot2::ggplot(
        data = df, ggplot2::aes(x_axis, y_axis, color=Dcolor, fill=Dcolor)
      ) + ggplot2::geom_violin(trim=FALSE,  alpha = 0.3)
      
      if (showPlotLegend){
        p <- p + ggplot2::geom_jitter(height = 0, size = as.numeric(dotsize))
      }
    }
  }
  ## Done plot logic                                                       ##
  ###########################################################################
  
  
  p <- p + ggplot2::xlab(x_axis) + ggplot2::ylab(y_axis)
  
  if (colorBy %in% nonNumCols ){
    # dfCol <- unique(df[,c(colorBy, "dotColor")])
    # colVec <- dfCol$dotColor
    # names(colVec) <- as.character(dfCol[,colorBy])
    # colVec <- colVec[colVec != ""]
    ## Colvec is provided if cols are non numeric
    
    
    p <- p + ggplot2::scale_colour_manual(colorBy ,values = colVec
    ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
    )
    
    if (plotLogic %in% c("ridgeplot","density", "histogram", "barchart", "violin")){
      p <- p + ggplot2::scale_fill_manual(colorBy ,values = colVec
      ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
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
    axis.text.x   = ggplot2::element_text(size=8, angle=90, vjust=0.5, hjust=1),
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
    p <- p + theme(legend.position = "none")
  } 
  
  
  p
}
##
##################




