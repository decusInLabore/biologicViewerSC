###############################################################################
## Load data required on startup                                             ##


###############################################################################
## Load key file                                                             ##
assembleKeyList <- function(){
    
    keyList <- list()
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
    
    
    keyList[["geneDefault"]] = as.character(dfkey$default)
    keyList[["host"]] <- as.character(dfkey$url)
    keyList[["user"]] <- as.character(dfkey$id)
    keyList[["DBpd"]] <- as.character(dfkey$id2)
    keyList[["dbname"]] <- as.character(dfkey$db)
    keyList[["coordinateTbName"]] <- as.character(dfkey$coordTb)
    keyList[["exprTbName"]] <- as.character(dfkey$exprTb)
    keyList[["geneID_TbName"]] <- as.character(dfkey$geneTb)
    
    pos <- grep("experiment", names(dfkey))
    if (length(pos) == 0){
        dfkey[["experiment"]] <- paste0("Experiment_", 1:nrow(dfkey))
    }
    
    ## Decide on database mode ##
    pos <- grep("dataMode", names(dfkey))
    if (length(pos) == 1){
        if (dfkey$dataMode == "SQLite"){
            keyList[["dataMode"]] <- dfkey$dataMode
        } else {
            keyList[["dataMode"]] <- "MySQL"
        }
    } else {
        keyList[["dataMode"]] <- "MySQL"
    }
    
    keyList[["dfkey"]] <- dfkey
    return(keyList)
}

## Done keyList                                                              ##
###############################################################################

###############################################################################
## Load parameter file                                                       ##

loadParameterFile <- function(){
    
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
            dfParam <- NULL
        }
        
    } else {
        parameterFileLoaded <- FALSE
    }
    
    
    return(dfParam)
    
}

##                                                                           ##
###############################################################################

###############################################################################
## Load color file                                                           ##

loadColorFile <- function(){
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
            dfColOptions <- NULL
        }
        
    } else {
        colorFileLoaded <- FALSE
    }
    
    return(dfColOptions)
}

##                                                                           ##
###############################################################################




###############################################################################
## Create dropdown menues                                                    ##
createDropdownMenuList <- function(){
    
    collectionList <- list()
    
    oldw <- getOption("warn")
    options(warn = -1)
    
    keyList <- assembleKeyList()
    
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
    
    
    query <- paste0("SELECT DISTINCT * FROM ", keyList[["coordinateTbName"]])
    dfCoordSel <- DBI::dbGetQuery(dbDB, query)
    DBI::dbDisconnect(dbDB)
    
    ## Add column for all values
    dfCoordSel[["all"]] <- "all"
    
    baseCols <- names(dfCoordSel)[names(dfCoordSel) != c("row_names", "cellID")]
    numCols <- names(dfCoordSel)[unlist(lapply(dfCoordSel, is.numeric))]
    numCols <- numCols[!(numCols %in% c("row_names", "cellID"))]
    
    
    collectionList[["numCols"]] <- c("lg10Expr", numCols)
    collectionList[["nonNumCols"]] <- baseCols[!(baseCols %in% numCols)]
    ##                                                                           ##
    ###############################################################################
    
    ###############################################################################
    ## Create order in which samples are displayed                               ##
    dfColOptions <- loadColorFile()
    
    if (!is.null(dfColOptions)){
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
            colSel <- c(
                grep("orig_ident", names(dfCoordSel)),
                grep("sampleName", names(dfCoordSel))
            )
            
            if (length(colSel) > 0){
                conditionVec <- unique(sort(dfCoordSel[,colSel[1]])) 
            } else {
                conditionVec <- unique(sort(dfCoordSel[,1]))  
            }
        }
    }
    
    
    Nsamples <- length(conditionVec)
    
    collectionList[["conditionVec"]] <- conditionVec
    collectionList[["Nsamples"]] <- Nsamples
    
    ##                                                                           ##
    ###############################################################################
    
    ###############################################################################
    ## Select Display options                                                    ##       
    
    ## Create x and y axis selections if no parameterfile is loaded ##
    
    dfParam <- loadParameterFile()
    
    if (!is.null(dfParam)){
        Xsel <- as.vector(dfParam[dfParam$menuName == "x_axis","colSel"])
        Ysel <- as.vector(dfParam[dfParam$menuName == "y_axis","colSel"])
        
        xDisplayName <- gsub("_", " ", unique(dfParam[dfParam$menuName == "x_axis", "displayName"]))
        yDisplayName <- gsub("_", " ", unique(dfParam[dfParam$menuName == "y_axis", "displayName"]))
        
        
    } else {
        XYsel <- c(
            collectionList[["numCols"]], 
            collectionList[["nonNumCols"]]
        )
        
        
        ## Reorder
        headSel <- c(
            XYsel[grep("UMAP_", toupper(XYsel))],
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
    if (length(defaultX %in% toupper(Xsel)) != 1){
        defaultX <- XYsel[2]
    }
    
    defaultY <- "UMAP_2"
    if (length(defaultX %in% toupper(Ysel)) != 1){
        defaultX <- XYsel[3]
    }
    
    firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
    }
    
    
    ## Add to dropdownlist 
    dropDownList <- list()
    
    Xdisplay <- gsub("_", " ", Xsel)
    Xdisplay <- firstup(Xdisplay)
    names(Xsel) <- Xdisplay
    
    dropDownList[["x_axis"]] <- list(
        "displayName" = xDisplayName,
        "selOptions" = Xsel,
        "selDisplayOptions" = Xdisplay,
        "default" = defaultX
    )
    
    ## Add to dropdownlist 
    Ysel <- c("Densityplot", "Histogram", Ysel)
    
    Ydisplay <- gsub("_", " ", Ysel)
    Ydisplay <- firstup(Ydisplay)
    names(Ysel) <- Ydisplay
    
    dropDownList[["y_axis"]] <- list(
        "displayName" = yDisplayName,
        "selOptions" = Ysel,
        "selDisplayOptions" = Ydisplay,
        "default" = defaultY
    )
    
    
    ##                                                                           ##
    ###############################################################################
    
    ###############################################################################
    ## Set color options                                                         ##
    
    if(!is.null(dfParam)){
        ## If paramsfile is loaded 
        allColorOptions <- unique(dfParam[dfParam$displayName == "Color Plots By", "colSel"])
        names(allColorOptions) <- gsub("_", " ", unique(dfParam[dfParam$displayName == "Color Plots By", "colSel"]))
    } else {
        allColorOptions <- c(
            collectionList$nonNumCols,
            collectionList$numCols
        )
        
        names(allColorOptions) <- firstup(
            gsub(
                "_", " ", allColorOptions
            )
        )
        
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
    
    firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
    }
    
    colorDisplayOptions = names(allColorOptions)
    
    ## Add to dropdownlist 
    dropDownList[["colorBy"]] <- list(
        "displayName" = cDisplayName,
        "selOptions" = allColorOptions,
        "selDisplayOptions" = colorDisplayOptions,
        "default" = defaultCol
    )
    
    
    ###############################################################################
    ## Set split options                                                         ##
    if (!is.null(dfParam)){
        ## If paramsfile is loaded 
        splitOptions <- unique(dfParam[dfParam$menuName == "splitPlotsBy", "colSel"])
        names( splitOptions) <- gsub("_", " ", unique(dfParam[dfParam$menuName == "splitPlotsBy", "colSel"]))
        
    } else {
        splitOptions <- c(
            collectionList$nonNumCols,
            collectionList$numCols
        )
        
    }
    
    splitOptions <- splitOptions[splitOptions %in% names(dfCoordSel)]
    
    ## Remove all split options with more than 20 options ##
    Nopt <- apply(dfCoordSel[,splitOptions], 2, function(x) length(unique(x)))
    Nopt <- sort(Nopt[Nopt < 42], decreasing = F)
    
    
    splitOptions2 <- as.vector(names(Nopt))
    
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
            grep("all", splitOptions2),
            grep("sampleName", splitOptions2),
            grep("meta_", tolower(splitOptions2)),
            grep("clusterName", splitOptions2),
            grep("subClusterName", splitOptions2)
        )  
    } else {
        headVec <- c(
            grep(names(dfCoordSel)[pos[1]], splitOptions2),
            grep("meta_", tolower(splitOptions2)),
            grep("all", splitOptions2),
            grep("clusterName", splitOptions2),
            grep("subClusterName", splitOptions2)
        )  
    }
    
    
    
    if (length(headVec) > 0){
        headOptions <- splitOptions2[headVec]
        restVec <- splitOptions2[-headVec]
        splitOptions2 <- c(
            headOptions,
            restVec
        )
    }
    
    names(splitOptions2) <- splitOptions2
    names(splitOptions2) <- gsub("all", "None", names(splitOptions2) )
    
    
    ## check if all column names are valid ##
    check <- names(dfCoordSel)
    splitOptions2 <- splitOptions2[splitOptions2 %in% check]
    
    
    defaultS <- splitOptions2[1]
    
    sDisplayName <- "Split Plots By"
    
    firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
    }
    
    selDisplayOptions = gsub("_", " ", names(splitOptions2))
    selDisplayOptions <- gsub("all", "None", selDisplayOptions)
    names(splitOptions2) <- firstup(selDisplayOptions)
    
    
    ## Add to dropdownlist 
    dropDownList[["splitByColumn"]] <- list(
        "displayName" = sDisplayName,
        "selOptions" = splitOptions2,
        "selDisplayOptions" = selDisplayOptions,
        "default" = defaultS
    )
    ## Done setting split options                                                ##
    ###############################################################################
    
    collectionList[["dropDownList"]] <- dropDownList
     
    return(collectionList)  
}

## Done dropdownList                                                         ##
###############################################################################

###############################################################################
## Create utility list                                                       ##


###############################################################################
## Get all genes to list                                                     ##

getAllEntriesToList <- function(){
    ###############################################################################
    ## Query all genes to be listed                                              ##
    oldw <- getOption("warn")
    options(warn = -1)
    
    keyList <- assembleKeyList()
    
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
    
    
    
    query <- paste0("SELECT DISTINCT gene FROM ", keyList[["geneID_TbName"]])
    dfGene <- DBI::dbGetQuery(dbDB, query)
    
    allGenes <- as.vector(dfGene[,"gene"])
    allGenes <- c(keyList[["geneDefault"]], allGenes)
    RMySQL::dbDisconnect(dbDB) 
    
    return(allGenes)
}

## Done                                                                      ##
###############################################################################





###############################################################################
## Collection function                                                       ##

loadStartUpData <- function(){
    
    ###############################################################################
    ## Assemble utility list                                                     ##
    utilityList <- list()
    utilityList[["allGenes"]] <- getAllEntriesToList()
    
    collectionList <- createDropdownMenuList()
    utilityList[["dropDownList"]] <- collectionList$dropDownList
    utilityList[["numCols"]] <- collectionList$numCols
    utilityList[["nonNumCols"]] <- collectionList$nonNumCols
    ## Done                                                                      ##
    ###############################################################################
    
    startUpList <- list()
    startUpList[["utilityList"]] <- utilityList
    
    startUpList[["keyList"]] <- assembleKeyList()
    startUpList[["dfParam"]] <- loadParameterFile()
    startUpList[["dfColOptions"]] <- loadColorFile()
    
    return(startUpList)
}

## Done                                                                      ##
###############################################################################

