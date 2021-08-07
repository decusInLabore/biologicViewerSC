#' The application server-side
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import DBI
#' @import RMySQL
#' @import ggplot2
#' @import colourpicker
 


library(RMySQL)
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
    #library(ggplot2)
    
    
    
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
    
    ###########################################################################
    ## Determine split options                                               ##
    splitOptions <- names(df)
    
    rmVec <- c(
        grep("orig_", splitOptions),
        grep("sampleID", splitOptions),
        grep("old_ident", splitOptions),
        grep("hmIdent", splitOptions),
        grep("color", tolower(splitOptions))
        
    )
    
    if (length(rmVec) > 0){
        splitOptions <- splitOptions[-rmVec]
    }
    
    
    ## Remove all split options with more than 20 options ##
    Nopt <- apply(df[,splitOptions], 2, function(x) length(unique(x)))
    Nopt <- sort(Nopt[Nopt < 25], decreasing = F)
    
    
    splitOptions <- as.vector(names(Nopt))
    ## Done                                                                  ##
    ###########################################################################
    
    
    nCellsTotal <- nrow(df)
    nExpr <- nrow(df[df$gene != 0,])
    percExpr <- 100*round(nrow(df[df$gene != 0,])/nCellsTotal, 3)
    qGene <- unique(na.omit(df$gene))
    qGene <- qGene[qGene != 0]
    
    plotInput <- reactive({
        
        if (colorBy %in% splitOptions ){
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
        
        
        
        
        p <- ggplot2::ggplot( data = df, ggplot2::aes(x_axis, y_axis, color=Dcolor)
        )+ ggplot2::geom_point(
            shape = 16,
            size = as.numeric(dotsize)
        ) + xlab(x_axis) + ylab(y_axis)
        
        if (colorBy %in% splitOptions ){
            dfCol <- unique(df[,c(colorBy, "dotColor")])
            colVec <- dfCol$dotColor
            names(colVec) <- as.character(dfCol[,colorBy])
            colVec <- colVec[colVec != ""]
            
            
            p <- p + ggplot2::scale_colour_manual(colorBy ,values = colVec
            ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
            )
        } else if (is.numeric( df$Dcolor )){
            if (minExpr < 0){
                p <- p + ggplot2::scale_color_gradient2("Expr",low= lowColor, mid = "white", high= dotcolor, midpoint = 0, limits=c(minExpr,maxExpr)
                )
                
            } else {
                p <- p + ggplot2::scale_color_gradient("Expr",low= lowColor, high= dotcolor, limits=c(minExpr,maxExpr)
                )
            }
            
        } else if (colorBy == "DF_Classification" & length(unique(df$Dcolor)) == 2) {
            p <- p + ggplot2::scale_colour_manual("Doublet Class",values = c("red","black")
            ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
            )
        } else if (colorBy == "all") {
            p <- p + ggplot2::scale_colour_manual("All",values = c("black")
            ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
            )
        }  else if (colorBy == "clusterName"){
            dfCol <- unique(df[,c("clusterName", "dotColor")])
            colVec <- dfCol$dotColor
            names(colVec) <- dfCol$clusterName
            colVec <- colVec[colVec != ""]
            
            
            p <- p + ggplot2::scale_colour_manual("Cluster Names" ,values = colVec
            ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
            )
        } else if (colorBy == "subClusterName"){  
            df$subClusterName <- gsub("^$", "Rest",df$subClusterName)
            dfCol <- unique(df[,c("subClusterName", "subClusterColor")])
            dfCol[dfCol$subClusterName == "Rest", "subClusterColor"] <- "#d3d3d3"
            colVec <- dfCol$subClusterColor
            names(colVec) <- dfCol$subClusterName
            
            colVec <- colVec[colVec != ""]
            p <- p + ggplot2::scale_colour_manual("Sub-cluster Names" ,values = colVec
            ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
            )
            
        } else if (colorBy == "sampleName"){  
            dfCol <- unique(df[,c("sampleName", "dotColor")])
            colVec <- dfCol$dotColor
            names(colVec) <- dfCol$sampleName
            colVec <- colVec[colVec != ""]
            p <- p + ggplot2::scale_colour_manual("Sample Names" ,values = colVec
            ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
            )
            
        }
        
        if (!is.numeric(df$x_axis)){
            p <- p + ggplot2::geom_jitter(height = 0) 
            
            
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
    #library(ggplot2)
    
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
    
    ###########################################################################
    ## Determine split options                                               ##
    splitOptions <- names(df)
    
    rmVec <- c(
        grep("orig_", splitOptions),
        grep("sampleID", splitOptions),
        grep("old_ident", splitOptions),
        grep("hmIdent", splitOptions),
        grep("color", tolower(splitOptions))
        
    )
    
    if (length(rmVec) > 0){
        splitOptions <- splitOptions[-rmVec]
    }
    
    
    ## Remove all split options with more than 20 options ##
    Nopt <- apply(df[,splitOptions], 2, function(x) length(unique(x)))
    Nopt <- sort(Nopt[Nopt < 25], decreasing = F)
    
    
    splitOptions <- as.vector(names(Nopt))
    ## Done                                                                  ##
    ###########################################################################
    
    
    nCellsTotal <- nrow(df)
    nExpr <- nrow(df[df$gene != 0,])
    percExpr <- 100*round(nrow(df[df$gene != 0,])/nCellsTotal, 3)
    qGene <- unique(na.omit(df$gene))
    qGene <- qGene[qGene != 0]
    
    #plotInput <- reactive({
        
        if (colorBy %in% splitOptions ){
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
        
        
        
        
        p <- ggplot2::ggplot( data = df, ggplot2::aes(x_axis, y_axis, color=Dcolor)
        )+ ggplot2::geom_point(
            shape = 16,
            size = as.numeric(dotsize)
        ) + xlab(x_axis) + ylab(y_axis)
        
        if (colorBy %in% splitOptions ){
            dfCol <- unique(df[,c(colorBy, "dotColor")])
            colVec <- dfCol$dotColor
            names(colVec) <- as.character(dfCol[,colorBy])
            colVec <- colVec[colVec != ""]
            
            
            p <- p + ggplot2::scale_colour_manual(colorBy ,values = colVec
            ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
            )
        } else if (is.numeric( df$Dcolor )){
            if (minExpr < 0){
                p <- p + ggplot2::scale_color_gradient2("Expr",low= lowColor, mid = "white", high= dotcolor, midpoint = 0, limits=c(minExpr,maxExpr)
                )
                
            } else {
                p <- p + ggplot2::scale_color_gradient("Expr",low= lowColor, high= dotcolor, limits=c(minExpr,maxExpr)
                )
            }
            
        } else if (colorBy == "DF_Classification" & length(unique(df$Dcolor)) == 2) {
            p <- p + ggplot2::scale_colour_manual("Doublet Class",values = c("red","black")
            ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
            )
        } else if (colorBy == "all") {
            p <- p + ggplot2::scale_colour_manual("All",values = c("black")
            ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
            )
        }  else if (colorBy == "clusterName"){
            dfCol <- unique(df[,c("clusterName", "dotColor")])
            colVec <- dfCol$dotColor
            names(colVec) <- dfCol$clusterName
            colVec <- colVec[colVec != ""]
            
            
            p <- p + ggplot2::scale_colour_manual("Cluster Names" ,values = colVec
            ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
            )
        } else if (colorBy == "subClusterName"){  
            df$subClusterName <- gsub("^$", "Rest",df$subClusterName)
            dfCol <- unique(df[,c("subClusterName", "subClusterColor")])
            dfCol[dfCol$subClusterName == "Rest", "subClusterColor"] <- "#d3d3d3"
            colVec <- dfCol$subClusterColor
            names(colVec) <- dfCol$subClusterName
            
            colVec <- colVec[colVec != ""]
            p <- p + ggplot2::scale_colour_manual("Sub-cluster Names" ,values = colVec
            ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
            )
            
        } else if (colorBy == "sampleName"){  
            dfCol <- unique(df[,c("sampleName", "dotColor")])
            colVec <- dfCol$dotColor
            names(colVec) <- dfCol$sampleName
            colVec <- colVec[colVec != ""]
            p <- p + ggplot2::scale_colour_manual("Sample Names" ,values = colVec
            ) + ggplot2::guides(col = guide_legend(override.aes = list(shape = 16, size = 5))
            )
            
        }
        
        if (!is.numeric(df$x_axis)){
            p <- p + ggplot2::geom_jitter(height = 0) 
            
            
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
    # })
## End of download functions   
}
# -------------------------------------------------------------------------

#dfkey <- read.delim("inst/exdata/connect/db.txt", header = T, sep="\t", stringsAsFactors = F)
#load("data/dfkey.rda")

FNkey <- "data/connect/db.txt"
FNrda <- "data/dfkey.rda"
if (file.exists(FNkey)){
    dfkey <- read.delim(FNkey, stringsAsFactors = F, sep="\t")
} else if (file.exists(FNrda)){
    load(FNrda)
} else {
    data("dfkey")
}



geneDefault = as.character(dfkey$default)

host <- as.character(dfkey$url)
user <- as.character(dfkey$id)
DBpd <- as.character(dfkey$id2)
dbname <- as.character(dfkey$db)
coordinateTbName <- as.character(dfkey$coordTb)
exprTbName <- as.character(dfkey$exprTb)
geneID_TbName <- as.character(dfkey$geneTb)

oldw <- getOption("warn")
options(warn = -1)

dbDB <- RMySQL::dbConnect(RMySQL::MySQL(), user = user, password = DBpd, host = host, dbname=dbname)
query <- paste0("SELECT DISTINCT gene FROM ", geneID_TbName)
dfGene <- DBI::dbGetQuery(dbDB, query)

allGenes <- as.vector(dfGene[,"gene"])
allGenes <- c(geneDefault, allGenes)
RMySQL::dbDisconnect(dbDB) 


###############################################################################
## Determine numeric and factorial display                                   ##
       
oldw <- getOption("warn")
options(warn = -1)
dbDB <- DBI::dbConnect(RMySQL::MySQL(), user = user, password = DBpd, host = host, dbname=dbname)
query <- paste0("SELECT DISTINCT * FROM ", coordinateTbName)
dfCoordSel <- DBI::dbGetQuery(dbDB, query)
DBI::dbDisconnect(dbDB)

dfCoordSel[["all"]] <- "all"

pos <- grep("sampleOrder", names(dfCoordSel))

if (length(pos) > 0){
    dfOrder <- unique(dfCoordSel[,c("sampleName", "sampleOrder")])
    dfOrder <- dfOrder[order(dfOrder$sampleOrder, decreasing = F),]
    conditionVec <- as.vector(dfOrder$sampleName)
} else {
    conditionVec <- unique(sort(dfCoordSel$sampleName))  
}

Nsamples <- length(conditionVec)


allOptions <- names(dfCoordSel)

rmNameVec <-c(
    "^DC",
    "uniquecellID",
    "hmIdent",
    "old_ident",
    "cellID", 
    "sample_group",
    "DF_pANN",
    "clusterColor",
    "sampleColor",
    "clustIdent",
    "G2M_Score",
    #"DM_Pseudotime",
    "^Sub_clusters_ExNeurons$",
    "sample_group_colors",
    "row_names",
    "sampleID"
)

rmVec <- as.vector(NULL, mode = "numeric")
for (i in 1:length(rmNameVec)){
    rmVec <- c(
        rmVec,
        grep(rmNameVec[i], names(dfCoordSel))
    )
}

XYsel <- allOptions
if (length(rmVec) > 0){
    XYsel <- XYsel[-rmVec]
}

## Reorder
XYsel <- c(
    XYsel[grep("UMAP_", XYsel)],
    XYsel[grep("tSNE_", XYsel)],
    XYsel[grep("sampleName", XYsel)],
    XYsel[grep("clusterName", XYsel)],
    XYsel[grep("ClusterTame", XYsel)],
    XYsel[grep("ClusterTest", XYsel)],
    XYsel[grep("PC_", XYsel)],
    XYsel[grep("DM_Pseudotime", XYsel)],
    XYsel[grep("meta_", XYsel)],
    #XYsel[grep("DF_Classification", XYsel)],
    XYsel[grep("nCount", XYsel)],
    XYsel[grep("nFeatures", XYsel)],
    XYsel[grep("nCount", XYsel)],
    XYsel[grep("percent", XYsel)],
    XYsel[grep("nCount", XYsel)],
    XYsel[grep("nCount", XYsel)]
)


#############################################
## Make Color Selections 
## Get color selection ##
allColorOptions <- c(
    #"Log10 Expresson" = "lg10Expr",
    #"DM Pseudotime"  = "DM_Pseudotime",
    "Sample" = "sampleName",
    "Cluster" = "clusterName",
    "Subcluster" = "subClusterName",
    # "WT vs. IDH" = "WT_IDH",
    "Gender" = "Gender",
    #  "Norm vs Hyp" = "Norm_Hyp",
    #  "Con Prad AZ" = "Con_Prad_AZ",
    "Cells From Tumor" = "CellFromTumor",
    "Patient" = "Patient",
    "Region" = "Region",
    "Article Cell Type" = "Article_Cell_Type",
    "Doublet Classification" = "DF_Classification" ,
    "nCount_RNA" = "nCount_RNA",
    "nFeature_RNA" = "nFeature_RNA",
    "percent_mt" = "percent_mt",
    "S Phase Score" = "S_Score",
    "G2M Score" = "G2M_Score",
    "Cell Cycle Phase" = "Phase",
    "Uniform" = "all"
)

colAddvec <- c(
    XYsel[grep("meta_", XYsel)],
    XYsel[grep("ClusterTestRes", XYsel)]
)

names(colAddvec) <- colAddvec

allColorOptions <- c(
    allColorOptions, 
    colAddvec
)


allColorOptions <- allColorOptions[allColorOptions %in% names(dfCoordSel)]
allColorOptions <- 
    c(
        "Log10 Expression" = "lg10Expr",
        allColorOptions
    )


splitOptions <- names(dfCoordSel)

rmVec <- c(
    grep("orig_", splitOptions),
    grep("sampleID", splitOptions),
    grep("old_ident", splitOptions),
    grep("hmIdent", splitOptions),
    grep("color", tolower(splitOptions))
    
)

if (length(rmVec) > 0){
    splitOptions <- splitOptions[-rmVec]
}


## Remove all split options with more than 20 options ##
Nopt <- apply(dfCoordSel[,splitOptions], 2, function(x) length(unique(x)))
Nopt <- sort(Nopt[Nopt < 51], decreasing = F)


splitOptions <- as.vector(names(Nopt))

factorOptions <- splitOptions
numericOptions <- allColorOptions[!(allColorOptions %in% factorOptions)]
numericOptions <- numericOptions[numericOptions != "lg10Expr"]
numericRes <- as.vector(NULL, mode = "character")
for (i in 1:length(numericOptions)){
    if (is.numeric(dfCoordSel[,numericOptions[i]])){
        numericRes <- c(
            numericRes,
            numericOptions[i]
        )
    }
}
numericOptions <- numericRes

## Done                                                                      ##
###############################################################################

app_server <- function(input, output, session) {
    
    
    
    
    
    
    ###############################################################################
    ## Load dfCoord from db                                                      ##
    
    createDfCoord <- reactive({
        oldw <- getOption("warn")
        options(warn = -1)
        dbDB <- RMySQL::dbConnect(RMySQL::MySQL(), user = user, password = DBpd, host = host, dbname=dbname)
        query <- paste0("SELECT DISTINCT * FROM ", coordinateTbName)
        dfCoordSel <- DBI::dbGetQuery(dbDB, query)
        
        RMySQL::dbDisconnect(dbDB)
        
        dfCoordSel$row_names <- NULL
        dfCoordSel[["all"]] <- "all"  
        
        # clusterCols <- unique(
        #   c(
        #     names(dfCoordSel)[grep("cluster", names(dfCoordSel))],
        #     names(dfCoordSel)[grep("Cluster", names(dfCoordSel))]
        #   )
        # )
        # 
        # 
        # 
        # if (length(clusterCols) > 0){
        #     for (m in 1:length(clusterCols)){
        #         clusters <- sort(unique(dfCoordSel[, clusterCols[m]]))
        #         tag <- paste0(clusterCols[m], "_number")
        #         dfCoordSel[[tag]] <- -1
        #         for (k in 1:length(clusters)){
        #             dfCoordSel[dfCoordSel[,clusterCols[m]] == clusters[k],  tag ] <- k
        #         }
        #     }
        # }
        # 
        # output$dev_text <- renderPrint({
        #   cat(getwd())
        # })
        
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
        
        dfDL[["dotColor"]] <- ""
        dfDL[["lg10Expr"]] <- "A1"
        
        if (input$colorBy == "sampleName"){
            selVec <- unique(c("sampleName", "sampleColor", "dotColor"))
            dfDL <- unique(dfDL[,selVec])
            dfDL$dotColor <- dfDL$sampleColor
            
        } else if (input$colorBy == "clusterName"){
            selVec <- unique(c("clusterName", "seurat_clusters", "clusterColor", "dotColor"))
            dfDL <- unique(dfDL[,selVec])
            dfDL$dotColor <- dfDL$clusterColor
        } else {
            selVec <- c(input$colorBy, "dotColor")
            dfDL <- unique(dfDL[,selVec])
            dfDL <- dfDL[order(dfDL[,1], decreasing = F),]
            dfDL[,input$colorBy] <- factor(dfDL[,input$colorBy])
            
            dfDL[["dotColor"]] <- scales::hue_pal()(nrow(dfDL))
            
        }
        
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
        dbDB <- RMySQL::dbConnect(RMySQL::MySQL(), user = user, password = DBpd, host = host, dbname=dbname)
        
        if ( is.null(input$gene) | input$gene == "" ){
            query <- paste0("SELECT * FROM ",exprTbName," WHERE gene = '",geneDefault,"'" )  
        } else {
            query <- paste0("SELECT * FROM ",exprTbName," WHERE gene = '",input$gene,"'" )
        }
        
        #query <- paste0("SELECT DISTINCT * FROM agl315_gene_expr_tb WHERE gene = 'GFAP'" )
        dfExprSel <- DBI::dbGetQuery(dbDB, query)
        dbDisconnect(dbDB)
        
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
        
        
        dfTemp <- merge(createDfCoord(), createDfExprSel(), by.x = "cellID", by.y="cellID", all=TRUE)
        dfTemp[is.na(dfTemp)] <- 0
        dfTemp <- data.frame(dfTemp, stringsAsFactors = FALSE)
        dfTemp$gene <- as.character(dfTemp$gene)
        
        conditionVec <- sort(unique(dfTemp[,input$splitByColumn]))  
        
        #######################################################################
        ## Check if custom colors are to be used                             ##
        
        dfTemp[["Dcolor"]] <- dfTemp[,input$colorBy]
        
        dfTemp[["dotColor"]] <- "#000000"
        
        if (input$colorBy == "clusterName"){
            dfTemp[["dotColor"]] <- dfTemp[["clusterColor"]]
        } else if (input$colorBy == "sampleName"){
            dfTemp[["dotColor"]] <- dfTemp[["sampleColor"]]
        }
        
        inInput <- names(input)[names(input) %in% unique(dfTemp[[input$colorBy]])]
        if (length(inInput) > 0){
            for (k in 1:length(inInput)){
                dfTemp[dfTemp[,input$colorBy] == inInput[k], "dotColor"] <- input[[inInput[[k]]]]
            }

        }
        
       
        ## Done                                                              ##
        #######################################################################
        
        #dfTemp[["Dcolor"]] <- dfTemp[,input$colorBy]
        
        # if (!(input$colorBy %in% c("lg10Expr"))){
        #     dfTemp$Dcolor <- factor(dfTemp$Dcolor) 
        # } else {
        #     dfTemp$Dcolor <- as.numeric(dfTemp$Dcolor)        
        # } 
        
        
        
        # if (input$x_axis == "clusterName"){
        #     clusters <- sort(unique(dfTemp[,input$x_axis]))
        #     
        #     dfTemp[["x_axis"]] <- dfTemp[,paste0( input$x_axis, "_number")]
        #    
        #     
        # } else {
        #     dfTemp[["x_axis"]] <- dfTemp[,input$x_axis]    
        # }
        
        dfTemp[["x_axis"]] <- dfTemp[,input$x_axis]   
        
        if (!is.numeric(dfTemp$x_axis)){
            dfTemp$x_axis <- factor(dfTemp$x_axis, levels = sort(unique(dfTemp$x_axis)))
        }
        
        dfTemp[["y_axis"]] <- dfTemp[,input$y_axis]
        # clusterColorColName <- "clusterName"
        # dfCol <- unique(dfTemp[,c("clusterName", "clusterColor")])
        # colVec <- dfCol$clusterColor
        # names(colVec) <- dfCol$clusterName
        # 
        # subClusterColorColName <- "subClusterName"
        # dfCol <- unique(dfTemp[,c("subClusterName", "subClusterColor")])
        # sColVec <- dfCol$subClusterColor
        # names(sColVec) <- dfCol$subClusterName
        # sColVec <- sColVec[sColVec != ""]
        # 
        # levels <- 
        # dfTemp[["Cluster"]] <- factor(dfTemp[,clusterColorColName], levels = sort(unique(dfTemp[,clusterColorColName])))   
        # #dfTemp$clusterName <- as.numeric(dfTemp$clusterName)
        
        
        # if (input$colorBy == "lg10Expr"){
        #     selVec <- unique(c( "gene", "lg10Expr", "x_axis", "y_axis", "Dcolor", "cellID", "sampleID", input$splitByColumn))
        # } else {
        #     selVec <- unique(c( "gene", "lg10Expr", "x_axis", "y_axis", "Dcolor", "cellID", "sampleID", input$colorBy, input$splitByColumn))
        # }
        # 
        
        
        
        #dfTemp <- dfTemp[,selVec]  
        dfTemp <- dfTemp[(dfTemp$x_axis != 0 | dfTemp$y_axis != 0),] 
        dfTemp
    })
    ##                                                                           ##
    ###############################################################################
    
    
    plot_select <- reactive({
        df <- createDfTemp()
        df[["all"]] <- "all"
        as.vector(unique(df[, input$splitByColumn]))
    })
    
    
    
    
    
    # library(DT)
    # 
    # output$table5 <- DT::renderDataTable({
    #     plot_data()[[1]]
    # }) 
    observe({
        updateSelectizeInput(session, 'gene', choices = allGenes, server = TRUE,selected=geneDefault) 
    })
    
    
    # onBookmarked(function(url) {
    #     updateQueryString(url)
    # })
    
    onRestored(function(state) {
        updateSelectizeInput(session, "gene", selected=state$input$gene, choices=allGenes, server=TRUE)
    })
    #https://github.com/rstudio/shiny/issues/1375
    # onRestore(function(state) {
    #     if (!is.null(state$input$addButton) && state$input$addButton > 0) {
    #         updateSelectizeInput(session, "gene", choices = allGenes,
    #                              selected = state$input$gene, server = TRUE)
    #     }
    # })
    # 
    
    
    observe({
        if (!(input$colorBy %in% numericOptions)){
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
    
    
    
    
    
    
    # toListen <- reactive({
    #     list(
    #         input$gene,
    #         input$x_axis,
    #         input$y_axis,
    #         input$splitByColumn,
    #         input$dotsize,
    #         input$colorBy,
    #         input$lowColor, 
    #         input$dotcolor,
    #         input$background,
    #         input$Epip6
    #     )
    # })
    
    
    plot_data_names <- reactive({
        dfTemp <- createDfTemp()
        
        plot_select <- sort(as.vector(unique(dfTemp[, input$splitByColumn])))
        
        wtPos <- unique(c(
            grep("wt", plot_select),
            grep("WT", plot_select),
            grep("Wt", plot_select),
            grep("Ctrl", plot_select),
            grep("CTRL", plot_select)
        ))
        
        if (length(wtPos) > 0){
            plot_select <- c(
                plot_select[wtPos],
                plot_select[-wtPos]
            )
        }
        
        plot_select
    })
    
    maxExpr <- reactive({
        dfTemp <- createDfTemp()
        
        if (is.numeric(dfTemp$Dcolor)){
            maxExpr <- max(as.vector(dfTemp$Dcolor))
        } else{
            maxExpr <- NULL
        }
        return(maxExpr)
    })
    
    plot_data <- reactive({
        dfTemp <- createDfTemp()
        
        plot_select <- plot_data_names()
        
        lapply(plot_select, function(x) dfTemp[dfTemp[,input$splitByColumn] == x,])
    })
    
    
    
    
    determinePlotDims <- reactive({
        dfTemp <- createDfTemp()
        
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
    })
    
    
    
    
    observeEvent(reactiveValuesToList(input), {
    #observeEvent(toListen(), {
        #req(!is.null(input$splitByColumn))
        req(plot_data())
        
        
        dimVec <- determinePlotDims()
        maxX = dimVec[2]
        minX = dimVec[1]
        maxY = dimVec[4]
        minY = dimVec[3]
        
        
        output$multi_plot_ui <- renderUI({
            
            lapply(seq_along(plot_data() ),
                   function(n) {
                       return(plot_prep_ui(paste0("n", n)))
                   })
        })
        
        
        
        
        lapply(seq_along(plot_data()),
               function(i){
                   callModule(plot_prep_server,
                              paste0("n", i),
                              df = plot_data()[[i]],
                              plot_name = paste0(plot_data_names()[i]), 
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
                              maxExpr = maxExpr(),
                              showPlotLegend = input$showPlotLegend
                   )
               }
        )
        
        ## Make plot list for download
        res <- lapply(seq_along(plot_data()),
                      function(i){
                          callModule(plot_prep_server_dl,
                                     paste0("n", i),
                                     df = plot_data()[[i]],
                                     plot_name = paste0(plot_data_names()[i]), 
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
                                     maxExpr = maxExpr(),
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
        
        # for (i in seq_along(input$selected_sample)) {
        #   callModule(plot_prep_server,
        #              paste0("n", i),
        #              spec = plot_data()[[i]],
        #              plot_name = i)
        # }
        
    }
    
    
    
    
    )
    
    
}
