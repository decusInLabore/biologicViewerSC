###############################################################################
## Load key file                                                             ##
getDataAccess <- function(){
    
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
    ## Done Data Access Module                                                   ##
    ###############################################################################
    return(keyList)
}


##                                                                           ##
###############################################################################


###############################################################################
## Load parameter file                                                       ##

loadParameterFile <- function(){
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
            dfParam <- NULL
        }
        
    } else {
        parameterFileLoaded <- FALSE
    }
    
    
    return(dfParam)
    ## Done                                                                      ##
    ###############################################################################
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
## Get all genes to list                                                     ##

getAllEntriesToList <- function(){
    ###############################################################################
    ## Query all genes to be listed                                              ##
    oldw <- getOption("warn")
    options(warn = -1)
    
    keyList <- getDataAccess()
    
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
    
    ## Done gene query                                                           ##
    ###############################################################################
    return(allGenes)
}

###############################################################################
## Create dropdown menues                                                    ##
createDropdownMenuList <- function(){
    oldw <- getOption("warn")
    options(warn = -1)
    
    keyList <- getDataAccess()
    
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
    
    numCols <- c("lg10Expr", names(dfCoordSel)[unlist(lapply(dfCoordSel, is.numeric))])
    nonNumCols <- names(dfCoordSel)[!(names(dfCoordSel) %in% numCols)]
    ##                                                                           ##
    ###############################################################################
    
    ###############################################################################
    ## Create order in which samples are displayed                               ##
    
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
            conditionVec <- unique(sort(dfCoordSel$sampleName))  
        }
    }
    
    
    Nsamples <- length(conditionVec)
    
    ##                                                                           ##
    ###############################################################################
    
    ###############################################################################
    ## Select Display options                                                    ##       
    
    ## Create x and y axis selections if no parameterfile is loaded ##
    if (!is.null(dfParam)){
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
    
    if(!is.null(dfParam)){
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
    
    firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
    }
    
    colorDisplayOptions = gsub("_", " ", names(allColorOptions))
    colorDisplayOptions = gsub("all", "Uniform", colorDisplayOptions)
    names(allColorOptions) <- firstup(colorDisplayOptions)
    
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
    
    firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
    }
    
    selDisplayOptions = gsub("_", " ", names(splitOptions))
    selDisplayOptions <- gsub("all", "None", selDisplayOptions)
    names(splitOptions) <- firstup(selDisplayOptions)
    
    
    ## Add to dropdownlist 
    dropDownList[["splitByColumn"]] <- list(
        "displayName" = sDisplayName,
        "selOptions" = splitOptions,
        "selDisplayOptions" = selDisplayOptions,
        "default" = defaultS
    )
    ## Done setting split options                                                ##
    ###############################################################################
    
    
    dropDownList[["numCols"]] <- numCols
    
    ## Remove non-numerical entries with more than 100 options ##
    Nopt <- apply(dfCoordSel[,nonNumCols], 2, function(x) length(unique(x)))
    Nopt <- sort(Nopt[Nopt < 100], decreasing = F)
    
    dropDownList[["nonNumCols"]] <- names(Nopt)
    
    ## Done with color options                                                   ##
    ###############################################################################
    return(dropDownList)  
}

## Done                                                                      ##
###############################################################################
