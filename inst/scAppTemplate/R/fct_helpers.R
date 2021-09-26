###############################################################################
## Create plot namespace                                                     ##

plot_prep_ui <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("my_plot"), 
               width = "100%")
  )
}

##                                                                           ##
###############################################################################


###############################################################################
##                                                                           ##


plot_prep_server <- function(
  input,
  output,
  session, 
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
  
  plotInput <- reactive({
    startUpList <- golem::get_golem_options(which = "startUpList")
    nonNumCols <- startUpList$utilityList$nonNumCols
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
    startUpList <- golem::get_golem_options(which = "startUpList")
    nonNumCols <- startUpList$utilityList$nonNumCols
    if (colorBy %in% nonNumCols ){
        dfCol <- unique(df[,c(colorBy, "dotColor")])
        colVec <- dfCol$dotColor
        names(colVec) <- as.character(dfCol[,colorBy])
        colVec <- colVec[colVec != ""]
        
        
        p <- p + ggplot2::scale_colour_manual(colorBy ,values = colVec
        ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
        )
        
        if (plotLogic %in% c("ridgeplot","density", "histogram")){
          p <- p + ggplot2::scale_fill_manual(colorBy ,values = colVec
          ) + ggplot2::guides(
              col = guide_legend(
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
      p <- p + theme(legend.position = "none")
    } 
    
    
    p
  })
  
  output$my_plot <- renderPlot({
    
    print(plotInput())
    
  })
  
  
  
  
}


plot_prep_server_dl <- function(
  input,
  output,
  session, 
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
      ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
      )
      
      if (plotLogic %in% c("ridgeplot","density", "histogram")){
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
      p <- p + theme(legend.position = "none")
    } 
    
    
    p
  #})
  # })
  ## End of download functions   
}
##                                                                           ##
###############################################################################