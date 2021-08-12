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
## Load category color file if available                                     ##
FNcolParameters <- "parameters/colorParameters.txt"

if (file.exists(FNcolParameters)){
    dfColOptions <- read.delim(
        FNcolParameters, 
        header = T, 
        sep = "\t",
        stringsAsFactors = F
    )
    
    colorFileLoaded <- TRUE
    ## Check file integrity ##
    if (!(sum(names(dfColOptions) %in% c("menuName", "displayName", "colSel", "displayOrder")))){
        rm(dfColOptions)
        colorFileLoaded <- FALSE
    }
    
} else {
    colorFileLoaded <- FALSE
}


##                                                                           ##
###############################################################################

###############################################################################
## Data Access Module                                                        ##

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

pos <- grep("dataMode", names(dfkey))
if (length(pos) == 1){
    if (dfkey$dataMode == "SQLite"){
        dataMode <- dfkey$dataMode
    } else {
        dataMode <- "MySQL"
    }
} else {
    dataMode <- "MySQL"
}

## Done Data Access Module                                                   ##
###############################################################################

oldw <- getOption("warn")
options(warn = -1)

if (dataMode == "SQLite"){
    
    dbDB <- DBI::dbConnect(
        drv = RSQLite::SQLite(),
        dbname=dbname
    )
    
} else {
    
    dbDB <- DBI::dbConnect(
        drv = RMySQL::MySQL(),
        user = user, 
        password = DBpd, 
        host = host, 
        dbname=dbname
        
    )
    
}

#dbDB <- RMySQL::dbConnect(RMySQL::MySQL(), user = user, password = DBpd, host = host, dbname=dbname)



query <- paste0("SELECT DISTINCT gene FROM ", geneID_TbName)
dfGene <- DBI::dbGetQuery(dbDB, query)

allGenes <- as.vector(dfGene[,"gene"])
allGenes <- c(geneDefault, allGenes)
RMySQL::dbDisconnect(dbDB) 


###############################################################################
## Determine numeric and factorial display                                   ##
       
oldw <- getOption("warn")
options(warn = -1)
# dbDB <- DBI::dbConnect(RMySQL::MySQL(), user = user, password = DBpd, host = host, dbname=dbname)

if (dataMode == "SQLite"){
    
    dbDB <- DBI::dbConnect(
        drv = RSQLite::SQLite(),
        dbname=dbname
    )
    
} else {
    
    dbDB <- DBI::dbConnect(
        drv = RMySQL::MySQL(),
        user = user, 
        password = DBpd, 
        host = host, 
        dbname=dbname
        
    )
    
}

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


## Remove all split options with more than 50 options ##
Nopt <- apply(dfCoordSel[,splitOptions], 2, function(x) length(unique(x)))
Nopt <- sort(Nopt[Nopt < 51], decreasing = F)


splitOptions <- as.vector(names(Nopt))

factorOptions <- splitOptions
numericOptions <- allColorOptions[!(allColorOptions %in% factorOptions)]

    
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
        
        
        if (dataMode == "SQLite"){
            
            dbDB <- DBI::dbConnect(
                drv = RSQLite::SQLite(),
                dbname=dbname
            )
            
        } else {
            
            dbDB <- DBI::dbConnect(
                drv = RMySQL::MySQL(),
                user = user, 
                password = DBpd, 
                host = host, 
                dbname=dbname
                
            )
            
        }
        
        #dbDB <- RMySQL::dbConnect(RMySQL::MySQL(), user = user, password = DBpd, host = host, dbname=dbname)
        query <- paste0("SELECT DISTINCT * FROM ", coordinateTbName)
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
        
        #######################################################################
        ## Check if colors are available                                     ##
        colorAnnoFound <- FALSE
        
        if (colorFileLoaded){
            dfPlotCols <- dfColOptions[dfColOptions$menuName == input$colorBy, c("colOption", "colSel")]
            
            if (nrow(dfPlotCols) > 0){
                checkDLvec <- sort(na.omit(as.vector(unique(dfDL[,input$colorBy]))))
                checkColvec <- sort(na.omit(unique(dfPlotCols$colOption)))
                if (identical(checkDLvec, checkColvec)){
                    dfAddCol <- unique(dfPlotCols)
                    names(dfAddCol) <- gsub("colOption", input$colorBy, names(dfAddCol))
                    names(dfAddCol) <- gsub("colSel", "dotColor", names(dfAddCol))
                    
                    dfDL <- merge(
                        dfDL, 
                        dfAddCol, 
                        by.x = input$colorBy,
                        by.y = input$colorBy, 
                        all = TRUE
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
        
        # if (input$colorBy == "sampleName"){
        #     selVec <- unique(c("sampleName", "sampleColor", "dotColor"))
        #     dfDL <- unique(dfDL[,selVec])
        #     dfDL$dotColor <- dfDL$sampleColor
        #     
        # } else if (input$colorBy == "clusterName"){
        #     selVec <- unique(c("clusterName", "seurat_clusters", "clusterColor", "dotColor"))
        #     dfDL <- unique(dfDL[,selVec])
        #     dfDL$dotColor <- dfDL$clusterColor
        # } else {
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
        #dbDB <- RMySQL::dbConnect(RMySQL::MySQL(), user = user, password = DBpd, host = host, dbname=dbname)
        
        if (dataMode == "SQLite"){
            
            dbDB <- DBI::dbConnect(
                drv = RSQLite::SQLite(),
                dbname=dbname
            )
            
        } else {
            
            dbDB <- DBI::dbConnect(
                drv = RMySQL::MySQL(),
                user = user, 
                password = DBpd, 
                host = host, 
                dbname=dbname
                
            )
            
        }
        
        if ( is.null(input$gene) | input$gene == "" ){
            query <- paste0("SELECT * FROM ",exprTbName," WHERE gene = '",geneDefault,"'" )  
        } else {
            query <- paste0("SELECT * FROM ",exprTbName," WHERE gene = '",input$gene,"'" )
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
        
        
        dfTemp <- merge(createDfCoord(), createDfExprSel(), by.x = "cellID", by.y="cellID", all=TRUE)
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
        
        
        
        dfTemp[["x_axis"]] <- dfTemp[,input$x_axis]   
        
        if (!is.numeric(dfTemp$x_axis)){
            dfTemp$x_axis <- factor(dfTemp$x_axis, levels = sort(unique(dfTemp$x_axis)))
        }
        
        dfTemp[["y_axis"]] <- dfTemp[,input$y_axis]
        
        
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
        updateSelectizeInput(session, 'gene', choices = allGenes, server = TRUE, selected=geneDefault) 
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
