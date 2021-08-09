

###############################################################################
## seuratObjectToLocalViewer                                                 ##

#' @title seuratObjectToLocalViewer
#'
#'
#' @param project_id Project id
#' @param OsC Seurat object
#' 
#' @import Seurat RMySQL RSQLite biologicSeqTools
#' @export

seuratObjectToLocalViewer <- function(
  project_id = "testApp",
  projectPath = "./",
  OsC = NULL,
  dataMode = "SQLite"
  #host = host,
  #user = db.user,
  #password = db.pwd
  
    
){
    ###############################################################################
    ## Create Single-cell Application                                            ##
    
    shinyBasePath <- projectPath
    shinyProjectPath <-paste0(
      shinyBasePath, 
      project_id, "_app"
    )
    if (!dir.exists(shinyProjectPath)){
      dir.create(shinyProjectPath)
    }
    
    shinyProjectPath <-paste0(
      shinyBasePath, 
      project_id, "_app"
    )
    if (!dir.exists(shinyProjectPath)){
      dir.create(shinyProjectPath)
    }
    
    shinyDataPath <-paste0(shinyProjectPath, "/data/connect/")
    if (!dir.exists(shinyDataPath)){
      dir.create(shinyDataPath, recursive = T)
    }
    
    shinyLocalDbPath <-paste0(shinyProjectPath, "/data/")
    if (!dir.exists(shinyLocalDbPath)){
      dir.create(shinyDataPath, recursive = T)
    }
    
    
    shinyParamPath <-paste0(shinyProjectPath, "/data/parameters/")
    if (!dir.exists(shinyParamPath)){
      dir.create(shinyParamPath)
    }
    
    
    primDataDB = paste0(shinyLocalDbPath, "appDataDb.sqlite")
    
    
    PCAdbTableName <- paste0(
      project_id, 
      "_PCA"
    )
    
    expDbTable <- paste0(project_id, "_gene_expr_tb")
    
    geneTb = paste0(project_id, "_geneID_tb")
    
    ## Done                                                                      ##
    ###############################################################################
    
    ###############################################################################
    ## Create Expr table                                                         ##
    
    print("Rendering expression data...")
    ## Running this function may take a minute or two, depending on the number of cells in your dataset
    
    dfExpr <-  biologicSeqTools::createDfExpr(
      obj = OsC,
      assay = "RNA",
      #slot = "data",
      geneSel = NULL
    ) 
    
    ## In Sqlite version load via function ##
    
    
    
    ## Adding gene_ids 
    
    library("dplyr")
    
    dfIDTable <- dfExpr %>% 
      select(gene) %>% 
      distinct() %>% 
      mutate(gene_id = row_number())
    
    
    ## Upload expression table to database 
    
    print(paste0("Database to be used: ", primDataDB))
    print(paste0("Database table name to be used: ", paste0(project_id, "_geneID_tb")))
    
    biologicSeqTools::upload.datatable.to.database(
      #host = host,
      #user = db.user,
      #password = db.pwd,
      prim.data.db = primDataDB,
      dbTableName = geneTb,
      df.data = dfIDTable,
      db.col.parameter.list = list(
        "BIGINT(8) NULL DEFAULT NULL" = c("row_names"),
        "VARCHAR(100) CHARACTER SET utf8 COLLATE utf8_general_ci" = c("gene"),
        "INT(8) NULL DEFAULT NULL" = c("gene_id")
      ),
      new.table = T,
      cols2Index = c("gene"),
      mode = dataMode  # Options: "MySQL" and "SQLite"
    )
    
    
    ## Rearrange expression talbe
    ## This step may take a couple of minutes in a large dataset
    dfExpr <- dfExpr %>% 
      rename(condition = cellID)  %>%  
      mutate(lg10Expr = round(lg10Expr, 3)) %>% 
      arrange(gene) 
    
    
    
    
    
    
    
    
    ## Upload expression table to database 
    
    print(paste0("Database to be used: ", primDataDB))
    print(paste0("Database table name to be used: ", expDbTable))
    
    colCatList <- biologicSeqTools::inferDBcategories(dfExpr)
    
    biologicSeqTools::upload.datatable.to.database(
      #host = host,
      #user = db.user,
      #password = db.pwd,
      prim.data.db = primDataDB,
      dbTableName = expDbTable,
      df.data = dfExpr,
      db.col.parameter.list = colCatList,
      new.table = T,
      cols2Index = c("gene"),
      #indexName = c("idx_gene_exp"),
      mode = dataMode  # Options: "MySQL" and "SQLite"
    )
    
    
    ###############################################################################
    ## Create Metadata table                                                     ##
    print("Rendering Metadata...")
    
    
    dfCoord <- biologicSeqTools::createDfCoord(OsC)
    
    
    dupTest <- duplicated(toupper(names(dfCoord)))
    
    if (sum(dupTest) > 0 ){
      names(dfCoord)[duplicated(toupper(names(dfCoord)))] <- paste0( names(dfCoord)[duplicated(toupper(names(dfCoord)))], "_B")
    }
    
    dfCoord <- dfCoord [,!(duplicated(names(dfCoord)))]
    
    dfdbTable <- dfCoord
    
    pos <- grep("integrated_", names(dfdbTable))
    if (length(pos) > 0){
      dfdbTable <- dfdbTable[,-pos]
    }
    names(dfdbTable) <- gsub("\\.", "_", names(dfdbTable))
    
    write.table(
      dfdbTable,
      "temp.txt",
      row.names = F,
      sep = "\t"
    )
    
    dfdbTable <- read.delim(
      "temp.txt",
      sep = "\t",
      stringsAsFactors = F
    )
    
    if (file.exists("temp.txt")){
      unlink("temp.txt")
    }
    
    dfdbTable[is.na(dfdbTable)] <- ""
    
    
    
    columnDBcategoryList <- biologicSeqTools::inferDBcategories(dfData=dfdbTable)
    
    
    biologicSeqTools::upload.datatable.to.database(
      host = host,
      user = db.user,
      password = db.pwd,
      prim.data.db = primDataDB,
      dbTableName = PCAdbTableName,
      df.data = dfdbTable,
      db.col.parameter.list = columnDBcategoryList,
      new.table = TRUE,
      mode = dataMode
    )
    biologicSeqTools::killDbConnections()
    ##                                                                           ##
    ###############################################################################
    
    ###############################################################################
    ## Menu Options                                                              ##
    paramList <- biologicSeqTools::scanObjParams(OsC)
    
    ## Done                                                                      ##
    ###############################################################################
    
    ###############################################################################
    ## Select default splitBy options                                            ##
    
    ## Default all non-numeric columns
    numCols <- unlist(lapply(dfCoord, is.numeric))
    chrCols <- unlist(lapply(dfCoord, is.character))
    
    splitOptions <- names(dfCoord)[chrCols]
    
    splitOptions <- c(splitOptions[splitOptions == "all"], sort(splitOptions[splitOptions != "all"]))
    
    names(splitOptions) <- splitOptions
    names(splitOptions) <- gsub("all", "None", names(splitOptions))
    
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    
    names(splitOptions) <- sapply(names(splitOptions), function(x) firstup(x))
    
    
    paramList$splitPlotsBy <- NULL
    
    paramList$splitPlotsBy <- splitOptions
    
    ## Make None the first option
    
    ##                                                                           ##
    ###############################################################################
    
    ###############################################################################
    ## Select colorBy options                                                    ##
    
    
    ## Default all non-numeric columns
    numCols <- unlist(lapply(dfCoord, is.numeric))
    chrCols <- unlist(lapply(dfCoord, is.character))
    
    colorByOptions <- sort(names(dfCoord))
    
    ## Columns not to show in the color selection
    rmVec <- c(
      grep("^UMAP", toupper(colorByOptions)),
      grep("^TSNE", toupper(colorByOptions)),
      grep("^PC", toupper(colorByOptions)),
      grep("orig.ident", toupper(colorByOptions)),
      grep("cellID", toupper(colorByOptions))
      
    )
    
    if (length(rmVec) > 0){
      colorByOptions <- colorByOptions[-rmVec]
    }
    
    colorByOptionsPart1 <- c(
      grep("clusterName", colorByOptions),
      grep("seurat_clusters", colorByOptions),
      grep("sampleName", colorByOptions)
    )
    
    if (length(colorByOptionsPart1) > 0){
      colorByOptionsPart2 <- colorByOptions[-colorByOptionsPart1]
    }
    
    colorDisplayOptions <- c(
      "lg10Expr",
      colorByOptions[colorByOptionsPart1],
      sort(colorByOptionsPart2)
    )
    
    names(colorDisplayOptions) <- colorDisplayOptions
    names(colorDisplayOptions) <- gsub("all", "Uniform", names(colorDisplayOptions))
    names(colorDisplayOptions) <- gsub("lg10Expr", "log10 Expr", names(colorDisplayOptions))
    
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    
    names(colorDisplayOptions) <- sapply(names(colorDisplayOptions), function(x) firstup(x))
    names(colorDisplayOptions) <- gsub("_", " ", names(colorDisplayOptions))
    
    
    paramList$colorPlotsBy <- NULL
    
    paramList$colorPlotsBy <- colorDisplayOptions
    
    ## Order Options ##
    
    ##                                                                           ##
    ###############################################################################
    
    ###############################################################################
    ## Edit x-axis                                                               ##
    
    ## Done                                                                      ##
    ###############################################################################
    
    ###############################################################################
    ## Edit y-axis                                                               ##
    
    ## Done                                                                      ##
    ###############################################################################
    
    
    ###############################################################################
    ## Helper function                                                           ##
    
    createListChunk <- function(
      optionName = "Choose_an_X-axis",
      pList = paramList,
      listItem = "x_axis"
    ){
      dfTemp <- data.frame(
        menuName = rep( optionName,length(pList[listItem])), 
        displayName = names(pList[listItem]), 
        colSel=pList[[listItem]],
        displayOrder = 1:length(pList[[listItem]])
      )
      return(dfTemp)
    }
    
    ##                                                                           ##
    ###############################################################################
    
    dfTemp <- createListChunk(
      optionName = "Choose_an_X-axis",
      pList = paramList,
      listItem = "x_axis"
    )
    
    df <- dfTemp
    
    dfTemp <- createListChunk(
      optionName = "Choose_an_Y-axis",
      pList = paramList,
      listItem = "y_axis"
    )
    
    df <- rbind(
      dfTemp,
      df
    )
    
    dfTemp <- createListChunk(
      optionName = "Split_Plots_By",
      pList = paramList,
      listItem = "splitPlotsBy"
    )
    
    df <- rbind(
      dfTemp,
      df
    )
    
    dfTemp <- createListChunk(
      optionName = "Color_Plots_By",
      pList = paramList,
      listItem = "colorPlotsBy"
    )
    
    df <- rbind(
      dfTemp,
      df
    )
    
    write.table(
      df,
      paste0(shinyParamPath, "menuParameters.txt"),
      sep = "\t",
      row.names = F
    )
    
    ##                                                                           ##
    ###############################################################################
    
    ###############################################################################
    ## Make color table                                                          ##
    
    ## Done                                                                      ##
    ###############################################################################
    
    ## Done 
    ###############################################################################
    
    
    
    
    
    
    
    
    
    ########################
    ## Temp section
    
    copyVec <- paste0(
      system.file("scAppTemplate/",package="biologicViewerSC"),
      "/",
      list.files(system.file("scAppTemplate/",package="biologicViewerSC"))
    )
    
    file.copy(copyVec, shinyProjectPath,  recursive=TRUE)
    
    ##
    ##########################
    
    ################################
    ## 
    
    dfID <- data.frame(
      dataMode = dataMode,
      url = "",
      id = "",
      id2 = "",
      db =  gsub(paste0("./",project_id, "_app/"),"",primDataDB),
      coordTb = PCAdbTableName,
      exprTb = expDbTable,
      geneTb = geneTb,
      default = as.vector(dfExpr[1,"gene"])
    )
    
    FN <- paste0(shinyDataPath, "db.txt")
    write.table(dfID, FN, row.names = F, sep="\t")
    
    
    ##
    ################################
    
    pFolder <- paste0(getwd(), gsub("./", "/", shinyProjectPath))
    
    print("Local App generated in folder ", pFolder)
  
  
}

## End                                                                       ##
###############################################################################