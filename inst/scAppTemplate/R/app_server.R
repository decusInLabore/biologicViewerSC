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
###############################################################################
## Chapter I - Loading Parameters                                            ##
###############################################################################
###############################################################################

###############################################################################
## Load parameter file if available                                          ##
FNparameters <- "parameters/menuParameters.txt"

if (file.exists(FNparameters)){
    dfParam <- read.delim(
        FNparameters, 
        header = T, 
        sep = "\t",
        stringsAsFactors = F
    )
    
    parameterFileLoaded <- TRUE
    ## Check file integrity ##
    if (!(sum(names(dfParam) %in% c("menuName", "displayName", "colSel", "displayOrder")))){
        rm(dfParam)
        parameterFileLoaded <- FALSE
    }
    
} else {
    parameterFileLoaded <- FALSE
}

## Done                                                                      ##
###############################################################################

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
    ## Option to be used if installed as local app ##
    dfkey <- read.delim(FNkey, stringsAsFactors = F, sep="\t")
} else if (file.exists(FNrda)){
    load(FNrda)
} else {
    ## Option to be used if this is installed as R-package ##
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

## Decide on database mode ##
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

###############################################################################
## Query all genes to be listed                                              ##
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

## Done gene query                                                           ##
###############################################################################

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

numOptions <- c("lg10Expr", names(dfCoordSel)[unlist(lapply(dfCoordSel, is.numeric))])
numOptions <- c(
  "lg10Expr",
  numOptions
)




###############################################################################
## Create order in which samples are displayed                               ##

if (colorFileLoaded){
    dfSO <- dfColOptions[grep("^SAMPLENAME$", toupper(dfColOptions$menuName)), ]  
    if (nrow(dfSO) > 0){
        dfSO <- unique(dfSO[,c("colOption", "displayOrder")])
        dfSO <- dfSO[order(dfSO$displayOrder, decreasing = F),]
        conditionVec <- as.vector(dfSO$colOption)
    } else {
        conditionVec <- unique(sort(dfCoordSel$sampleName))  
    }
} else {
    pos <- grep("sampleOrder", names(dfCoordSel))
    
    if (length(pos) > 0){
        dfCoordSel <- unique(dfCoordSel[,c("sampleName", "sampleOrder")])
        dfCoordSel <- dfCoordSel[order(dfCoordSel$sampleOrder, decreasing = F),]
        conditionVec <- as.vector(dfCoordSel$sampleName)
    } else {
        conditionVec <- unique(sort(dfCoordSel$sampleName))  
    }
}


Nsamples <- length(conditionVec)

##                                                                           ##
###############################################################################


###############################################################################
## Select Display options                                                    ##       

## Create x and y axis selections if no parameterfile is loaded ##
if (parameterFileLoaded){
  Xsel <- as.vector(dfParam[dfParam$menuName == "x_axis","colSel"])
  Ysel <- as.vector(dfParam[dfParam$menuName == "y_axis","colSel"])
  
  xDisplayName <- gsub("_", " ", unique(dfParam[dfParam$menuName == "x_axis", "displayName"]))
  yDisplayName <- gsub("_", " ", unique(dfParam[dfParam$menuName == "y_axis", "displayName"]))
  
  
} else {
  allOptions <- names(dfCoordSel)
  
  ## Remove common undesirable colums ##
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
    "row_names"
  )
  
  rmVec <- as.vector(NULL, mode = "numeric")
  for (i in 1:length(rmNameVec)){
    rmVec <- c(
      rmVec,
      grep(rmNameVec[i], allOptions)
    )
  }
  
  XYsel <- allOptions
  if (length(rmVec) > 0){
    XYsel <- XYsel[-rmVec]
  }
  
  ## Reorder
  headSel <- c(
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
  
  restSel <- XYsel[!(XYsel %in% headSel)]
  
  XYsel <- c(
    "lg10Expr",
    headSel, 
    restSel
  )
  
  Xsel <- XYsel
  Ysel <- XYsel
  xDisplayName <- "Choose a X-axis"
  yDisplayName <- "Choose a Y-axis"
}

## check if all column names are valid ##
check <- c("lg10Expr", names(dfCoordSel))
Xsel <- Xsel[Xsel %in% check]
Ysel <- Ysel[Ysel %in% check]

defaultX <- "UMAP_1"
if (length(defaultX %in% Xsel) != 1){
  defaultX <- XYsel[2]
}

defaultY <- "UMAP_2"
if (length(defaultX %in% Ysel) != 1){
  defaultX <- XYsel[3]
}


## Add to dropdownlist 
dropDownList <- list()

dropDownList[["x_axis"]] <- list(
  "displayName" = xDisplayName,
  "selOptions" = Xsel,
  "selDisplayOptions" = gsub("_", " ", Xsel),
  "default" = defaultX
)

## Add to dropdownlist 
Ysel <- c(Ysel, "Densityplot", "Histogram")
dropDownList[["y_axis"]] <- list(
  "displayName" = yDisplayName,
  "selOptions" = Ysel,
  "selDisplayOptions" = gsub("_", " ", Ysel),
  "default" = defaultY
)


##                                                                           ##
###############################################################################

###############################################################################
## Set color options                                                         ##

if(parameterFileLoaded){
  ## If paramsfile is loaded 
  allColorOptions <- unique(dfParam[dfParam$displayName == "Color Plots By", "colSel"])
  names(allColorOptions) <- gsub("_", " ", unique(dfParam[dfParam$displayName == "Color Plots By", "colSel"]))
  
  
  
} else {
  allColorOptions <- c(
    "log10 Expression" = "lg10Expr",
    names(dfCoordSel)
  )
  
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
    "row_names"
  )
  
  rmVec <- as.vector(NULL, mode = "numeric")
  for (i in 1:length(rmNameVec)){
    rmVec <- c(
      rmVec,
      grep(rmNameVec[i], allColorOptions)
    )
  }
  
  
  if (length(rmVec) > 0){
    allColorOptions <- allColorOptions[-rmVec]
  }
  
  
  allColorOptions <- allColorOptions[allColorOptions %in% names(dfCoordSel)]
  
  allColorOptions <- c(
    "Log10 Expression" = "lg10Expr",
    allColorOptions
  )  
}


## Organise order ##
headVec <- unique(
  c(
    grep("LG10EXPR", toupper(allColorOptions)),
    grep("CLUSTERNAME", toupper(allColorOptions)),
    grep("CLUSTER", toupper(allColorOptions)),
    grep("SAMPLENAME", toupper(allColorOptions)),
    grep("META_", toupper(allColorOptions)),
    grep("CLUSTERNAME", toupper(allColorOptions)),
    grep("subClusterName", toupper(allColorOptions))
  )
)

if (length(headVec) > 0){
  headOptions <- allColorOptions[headVec]
  restVec <- allColorOptions[-headVec]
  allColorOptions <- c(
    headOptions,
    restVec
  )
}


pos <- grep("^all$", allColorOptions)

if (length(pos) > 0){
    names(allColorOptions)[pos] <- "Unicolor"
}



## check if all column names are valid ##
check <- c("lg10Expr", names(dfCoordSel))
allColorOptions <- allColorOptions[allColorOptions %in% check]


defaultCol <- "lg10Expr"
if (length(defaultCol %in% allColorOptions) != 1){
  defaultCol <- allColorOptions[1]
}

#cDisplayName <- gsub("_", " ", unique(dfParam[dfParam$menuName == "colorPlotsBy", "displayName"]))
cDisplayName <- "Color Plots By"

## Add to dropdownlist 
dropDownList[["colorBy"]] <- list(
  "displayName" = cDisplayName,
  "selOptions" = allColorOptions,
  "selDisplayOptions" = gsub("_", " ", names(allColorOptions)),
  "default" = defaultCol
)
## Done with color options                                                   ##
###############################################################################



###############################################################################
## Set split options                                                         ##
if (parameterFileLoaded){
  ## If paramsfile is loaded 
  splitOptions <- unique(dfParam[dfParam$menuName == "splitPlotsBy", "colSel"])
  names( splitOptions) <- gsub("_", " ", unique(dfParam[dfParam$menuName == "splitPlotsBy", "colSel"]))
  
} else {
  splitOptions <- names(dfCoordSel)
  
  
  
}

splitOptions <- splitOptions[splitOptions %in% names(dfCoordSel)]

## Remove all split options with more than 20 options ##
Nopt <- apply(dfCoordSel[,splitOptions], 2, function(x) length(unique(x)))
Nopt <- sort(Nopt[Nopt < 42], decreasing = F)


splitOptions <- as.vector(names(Nopt))

pos <- c(
  grep("sampleName", names(dfCoordSel)),
  grep("orig_ident", names(dfCoordSel))
)

if (length(pos) > 0){
  Nsamples <- length(unique(dfCoordSel[,pos[1]]))
} else {
  Nsamples <- 1000
}


if (Nsamples > 3 | nrow(dfCoordSel) < 5000){
  headVec <- c(
    grep("all", splitOptions),
    grep("sampleName", splitOptions),
    grep("meta_", tolower(splitOptions)),
    grep("clusterName", splitOptions),
    grep("subClusterName", splitOptions)
  )  
} else {
  headVec <- c(
    grep(names(dfCoordSel)[pos[1]], splitOptions),
    grep("meta_", tolower(splitOptions)),
    grep("all", splitOptions),
    grep("clusterName", splitOptions),
    grep("subClusterName", splitOptions)
  )  
}



if (length(headVec) > 0){
  headOptions <- splitOptions[headVec]
  restVec <- splitOptions[-headVec]
  splitOptions <- c(
    headOptions,
    restVec
  )
}

names(splitOptions) <- splitOptions
names(splitOptions) <- gsub("all", "None", names(splitOptions) )


numOptions <- c("lg10Expr", names(dfCoordSel)[unlist(lapply(dfCoordSel, is.numeric))])
numOptions <- c(
  "lg10Expr",
  numOptions
)



## check if all column names are valid ##
check <- c(names(dfCoordSel))
splitOptions <- splitOptions[splitOptions %in% check]


defaultS <- splitOptions[1]

sDisplayName <- "Split Plots By"

## Add to dropdownlist 
dropDownList[["splitByColumn"]] <- list(
  "displayName" = sDisplayName,
  "selOptions" = splitOptions,
  "selDisplayOptions" = gsub("_", " ", names(splitOptions)),
  "default" = defaultS
)
## Done setting split options                                                ##
###############################################################################


numOptions <- names(dfCoordSel)[!(names(dfCoordSel)) %in% splitOptions]
numOptions <- c(
  "lg10Expr",
  numOptions
)

##                                                                           ##
###############################################################################




## Remove all split options with more than 50 options ##
Nopt <- apply(dfCoordSel[,splitOptions], 2, function(x) length(unique(x)))
Nopt <- sort(Nopt[Nopt < 51], decreasing = F)


factorOptions <- splitOptions
numericOptions <- numOptions
# numericOptions <- numericOptions[numericOptions != "lg10Expr"]
#     
# numericRes <- as.vector(NULL, mode = "character")
# for (i in 1:length(numericOptions)){
#     if (is.numeric(dfCoordSel[,numericOptions[i]])){
#         numericRes <- c(
#             numericRes,
#             numericOptions[i]
#         )
#     }
# }
# numericOptions <- numericRes


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
        ## Create dummy log10 Expr
        dfDL[["lg10Expr"]] <- "A1"
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
                    
                    dfDL <- dplyr::full_join(
                        dfDL,
                        dfAddCol,
                        by= input$colorBy
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
        
        
        if(!colorAnnoFound) {
            dfDL[["dotColor"]] <- "#000000"
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
    
###############################################################################
## Database query for dfExpr                                                 ##
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
## Done db query                                                             ##
###############################################################################
    
###############################################################################
## Create dfTemp                                                             ##       
    createDfTemp <- reactive({
        
        dfTemp <- dplyr::full_join(
            createDfCoord(), 
            createDfExprSel(), 
            by="cellID"
        )
        
        dfTemp[["Dcolor"]] <- dfTemp[,input$colorBy]
        
        if (!(input$colorBy %in% numOptions)){
            dfTemp <- dplyr::full_join(
              dfTemp, 
              dfColorTable(), 
              by= input$colorBy
            )
        } else {
            dfTemp[["dotColor"]] <- "#000000"
        }
        
        
        #dfTemp2 <- merge(createDfCoord(), createDfExprSel(), by.x = "cellID", by.y="cellID", all=TRUE)
        dfTemp[is.na(dfTemp)] <- 0
        dfTemp <- data.frame(dfTemp, stringsAsFactors = FALSE)
        dfTemp$gene <- as.character(dfTemp$gene)
        
        conditionVec <- sort(unique(dfTemp[,input$splitByColumn]))  
        
        #######################################################################
        ## Check if custom colors are to be used                             ##
        
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
    
    
    
    observe({
        updateSelectizeInput(session, 'gene', choices = allGenes, server = TRUE, selected=geneDefault) 
    })
    
    onRestored(function(state) {
        updateSelectizeInput(session, "gene", selected=state$input$gene, choices=allGenes, server=TRUE)
    })
    
    
    ###################################################################
    ## Add dropdownselection                                         ##
    
    observe({
        if (length(dropDownList) > 0){
            
            output$dropDownPanel = renderUI({
                #dfColSel[["label"]] <- paste0(dfColSel[,nameCol], " ", labelID," Color" )
                dropdown_list <- lapply(1:length(dropDownList), function(i) {
                    # for each dynamically generated input, give a different name
                    displayID <- names(dropDownList)[i]
                    displayName <- as.vector(dropDownList[[displayID]][["displayName"]])
                    selOptions <- as.vector(dropDownList[[displayID]][["selOptions"]])
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
        if (!(input$colorBy %in% numOptions)){
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
        
        ## Call UI module
        output$multi_plot_ui <- renderUI({
            
        lapply(seq_along(plot_data),
                   function(n) {
                       return(plot_prep_ui(paste0("n", n)))
                   })
        })
        
        ## Call server side module
        lapply(
            seq_along(plot_data),
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
        
        #######################################################################
        ## Prepare plot for file output                                      ##
        
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
        
        # for (i in seq_along(input$selected_sample)) {
        #   callModule(plot_prep_server,
        #              paste0("n", i),
        #              spec = plot_data()[[i]],
        #              plot_name = i)
        # }
        
    }
    
    
    
    
    )
    
    
    } # end server 


