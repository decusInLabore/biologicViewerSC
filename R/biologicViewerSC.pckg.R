###############################################################################
## Add to Seurat metadata                                                    ##

#' @title addDf2seuratMetaData
#' @param obj Seurat object
#' @return paramerer list
#' @import Seurat
#' @export

setGeneric(
    name="addDf2seuratMetaData",
    def=function(obj, dfAdd) {
    #print(paste0("Dims before addition: ", dim(obj@meta.data)))
    ###########################################################################
    ## Detach packages that interfere with AddMetaData                       ##
    detach_package <- function(pkg, character.only = FALSE) {
        if(!character.only) {
          pkg <- deparse(substitute(pkg))
        }
        
        search_item <- paste("package", pkg, sep = ":")
        
        while(search_item %in% search()){
          detach(search_item, unload = TRUE, character.only = TRUE)
        }
    }
    
    interferenceVec <- c("SeuratDisk", "SeuratObject", "DESeq2")
    detachPacks <- interferenceVec[paste0("package:", interferenceVec) %in% search()]
    
    if (length(detachPacks) > 0){
        for (n in 1:length(detachPacks)){
          detach_package(detachPacks[n], character.only = TRUE)
        }
    }
    ##
    ###########################################################################
    
    for (i in 1:ncol(dfAdd)){
        addVec <- as.vector(dfAdd[,i])
        names(addVec) <- row.names(dfAdd)
        colName <- as.vector(names(dfAdd)[i])
        obj <- Seurat::AddMetaData(
            object = obj,
            metadata = addVec,
            colName
        )
    }
    
    #print(paste0("Dims after addition: ", dim(obj@meta.data)))
    #print(paste0("Meta data column names: ", paste(names(obj@meta.data), collapse = ", ")))
    return(obj)
  }
)

## Done adding to Seurat metadata                                            ##
###############################################################################

###############################################################################
#' Scan Parameters
#'
#' @title createDfCoord
#'
#' @param obj Seurat object
#' @return paramerer list
#' @import Seurat
#' @import dplyr
#' @import tibble
#' @export

setGeneric(
    name="createDfCoord",
    def=function(
        obj,
        params = NULL
    ) {
    
    if (is.null(params)){
      params <- biologicViewerSC::scanObjParams(obj)
    }
    
    
    ## Add reductions ##
    reds <- names(obj@reductions)
    
    for (i in 1:length(reds)){
        dfAdd <- data.frame(obj@reductions[[reds[i]]]@cell.embeddings)
        if (nrow(dfAdd) > 0){
            obj <- biologicViewerSC::addDf2seuratMetaData(obj = obj, dfAdd = dfAdd)
        }
    }
    
    ## Shape table
    
    
    ##
    dfdbTable <- obj@meta.data %>% 
        dplyr::select(-one_of("cellID"))  %>% 
        tibble::rownames_to_column(var = "cellID")
    
    pos <- grep("sampleID", names(dfdbTable))
    pos2 <- grep("orig.ident", names(dfdbTable))
    
    if (length(pos) == 0 | length(pos2) == 1){
      dfdbTable[["sampleID"]] <- dfdbTable[["orig.ident"]]
    } else {
      dfdbTable[["sampleID"]] <- "sampleID_TBD"
    }
    
    return(dfdbTable)
    
  }
)

###############################################################################

###############################################################################
#' Scan Parameters
#'
#' TBD.
#' @title createDfExpr
#' @param obj Seurat object
#' @return paramerer list
#' @import Seurat
#' @import tibble
#' @import tidyr
#' @import dplyr
#' @export


setGeneric(
    name="createDfExpr",
    def=function(
        obj,
        assay = "RNA",
        #slot = "data",
        geneSel = NULL
    ) {
    
    Seurat::DefaultAssay(obj) <- assay
    
    ## Breakdown into chunks to make it more memory friendly ##
    ## 2022 03 21
    cells <- row.names(OsC@meta.data)
    cellList <- split(cells, ceiling(seq_along(cells)/10000))
    
    ## Define helper function ##
    subset_fun <- function(obj, cellIDs) {
      z <- subset(obj, subset = cellID %in% cellIDs)[[assay]]@data %>%
        data.frame() %>%
        tibble::rownames_to_column(var = "gene") %>%
        tidyr::pivot_longer(
          !gene,
          names_to = "cellID",
          values_to = "lg10Expr"
        ) %>%
        filter(lg10Expr > 0)
        return(z)
    }
    ## End of helper function
    
    dfExpr <- purrr::map(cellList, function(x) subset_fun(obj=obj, cellIDs = x)) %>%
        dplyr::bind_rows() %>%
        data.frame()
    
    return(dfExpr)
  }
)


###############################################################################

###############################################################################
## Scan Seurat Parameters                                                    ##



#' @title scanObjParams
#'
#' TBD.
#'
#' @param obj Seurat object
#' @param NmaxSplit Maximum number of group members for a column/category to be listed in the splitBy column
#' @param NcatColorMax Value above which a numeric column is displayed as color gradient
#' @return paramerer list
#' @import Seurat scales
#' @export

setGeneric(
  name="scanObjParams",
  def=function(
    obj,
    NmaxSplit = 25,
    NcatColorMax = 40
  ) {
    
    ###########################################################################
    ## Helper function                                                       ##
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    ## Done                                                                  ##
    ###########################################################################
    
    obj@meta.data[["all"]] <- "all"
    
    addReductions <- function(red = "pca", obj, paramList=list()){
        reds <- names(obj@reductions)
        pos <- grep(red, reds)
        if (length(pos) > 0){
            dfTemp <- data.frame(obj@reductions[[red]]@cell.embeddings)
            if (nrow(dfTemp) > 0){
                paramList[[red]] <- names(dfTemp)
            }
        }
        return(paramList)
    }
    
    reds <- names(obj@reductions)
    t <- purrr::map(reds, function(x) addReductions(red = x, obj=obj))
    
    tempList <- list()
    
    tempList[["meta.data"]] <- c(names(obj@meta.data))
    names(tempList[["meta.data"]]) <- gsub("[.]", "_",gsub("meta.data", "", c(names(obj@meta.data))))
    reds <- names(obj@reductions)
    
    tempList <- c(
        tempList, 
        t
    )
    
    
    allOptions <- unlist(tempList, use.names = F)
    allOptions <- allOptions[allOptions != "all"]
    
    uPos <- allOptions[grep("UMAP", toupper(allOptions))]
    tPos <- allOptions[grep("TSNE", toupper(allOptions))]
    tPos <- allOptions[grep("^PC", toupper(allOptions))]
    clustPos <- allOptions[grep("CLUSTER", toupper(allOptions))]
    
    rmPos <- c(uPos, tPos, clustPos)
    
    restPos <- allOptions[!(allOptions %in% rmPos)]
    allOptions <- unique(c(uPos, tPos, clustPos, restPos))
    names(allOptions) <- gsub("[.]", "_", allOptions)
    
    
    
    
    XYsel <- allOptions
    # if (length(rmVec) > 0){
    #   XYsel <- XYsel[-rmVec]
    # }
    
    XYorder <- c(
      XYsel[grep("UMAP_", toupper(XYsel))],
      XYsel[grep("TSNE_", toupper(XYsel))],
      XYsel[grep("SAMPLENAME", toupper(XYsel))],
      XYsel[grep("CLUSTERNAME", toupper(XYsel))],
      XYsel[grep("CLUSTERTEST", toupper(XYsel))],
      XYsel[grep("PC_", toupper(XYsel))],
      XYsel[grep("DM_PSEUDOTIME", toupper(XYsel))],
      XYsel[grep("META_", toupper(XYsel))],
      #XYsel[grep("DF_Classification", XYsel)],
      XYsel[grep("NFEATURES", toupper(XYsel))],
      XYsel[grep("PERCENT", toupper(XYsel))],
      XYsel[grep("nCount", toupper(XYsel))]
    )
    
    XYrest <- sort(XYsel[!(XYsel %in% XYorder)])
    
    XYsel <- c(
        XYorder, 
        XYrest
    )
    
    names(XYsel) <- sapply(gsub("_", " ", names(XYsel)), firstup)
    
    paramList <- list()
    paramList[["x_axis"]] <- c("log10 Expr" = "lg10Expr", XYsel)
    paramList[["y_axis"]] <- c("log10 Expr" = "lg10Expr", XYsel)
    
    
    ## Now cretate split by options ##
    
    catOptions <- as.vector(NULL, mode = "character")
    fullCatOptions <- names(obj@meta.data)
    for (i in 1:ncol(obj@meta.data)){
      if (length(unique(obj@meta.data[,i])) <= NmaxSplit){
        catOptions <- c(
          catOptions,
          names(obj@meta.data)[i]
        )
      }
      
      
    }
    
    
    
    splitOptions <- catOptions
    
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
    Nopt <- apply(obj@meta.data[,splitOptions], 2, function(x) length(unique(x)))
    Nopt <- sort(Nopt[Nopt < NmaxSplit], decreasing = F)
    
    splitOptions <- as.vector(names(Nopt))
    
    Nsamples <- length(unique(obj@meta.data$sampleName))
    
    if (Nsamples > 3 | nrow(obj@meta.data) < 5000){
      headVec <- c(
        grep("all", splitOptions),
        grep("sampleName", splitOptions),
        grep("meta_", tolower(splitOptions)),
        grep("clusterName", splitOptions),
        grep("subClusterName", splitOptions)
      )  
    } else {
      headVec <- c(
        grep("sampleName", splitOptions),
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
    names(splitOptions) <- gsub("meta_", "", names(splitOptions) )
    names(splitOptions) <- gsub("META_", "", names(splitOptions) )
    names(splitOptions) <- gsub("meta_", "", names(splitOptions) )
    names(splitOptions) <- gsub("sampleName", "Sample Name", names(splitOptions) )
    names(splitOptions) <- gsub("clusterName", "Cluster Name", names(splitOptions) )
    names(splitOptions) <- gsub("all", "None", names(splitOptions) )
    
    names(splitOptions) <- sapply(names(splitOptions), firstup)
    
    paramList[["splitPlotsBy"]] <- splitOptions
    
    
    Nopt <- apply(obj@meta.data[,splitOptions], 2, function(x) length(unique(x)))
    Nopt <- sort(Nopt[Nopt < NcatColorMax], decreasing = F)
    
    numOptions <- names(obj@meta.data)[!(names(obj@meta.data)) %in% splitOptions]
    numOptions <- c(
      "lg10Expr",
      numOptions
    )
    
    ###############################################################################
    ## Select colorBy options                                                    ##
    
    
    ## Default all non-numeric columns
    numCols <- unlist(lapply(obj@meta.data, is.numeric))
    chrCols <- unlist(lapply(obj@meta.data, is.character))
    
    colorByOptions <- sort(names(obj@meta.data))
    
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
    
    
    
    colorDisplayOptions <- c(
      "lg10Expr",
      colorByOptions
    )
    
    names(colorDisplayOptions) <- colorDisplayOptions
    names(colorDisplayOptions) <- gsub("all", "Uniform", names(colorDisplayOptions))
    names(colorDisplayOptions) <- gsub("lg10Expr", "log10 Expr", names(colorDisplayOptions))
    
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    
    names(colorDisplayOptions) <- sapply(names(colorDisplayOptions), function(x) firstup(x))
    names(colorDisplayOptions) <- sapply(gsub("_", " ", names(colorDisplayOptions)), firstup)
    names(colorDisplayOptions) <- gsub("NFeature", "nFeature", names(colorDisplayOptions))
    names(colorDisplayOptions) <- gsub("NCount", "nCount", names(colorDisplayOptions))
    
    paramList[["colorPlotsBy"]] <- colorDisplayOptions
    
    ##                                                                           ##
    ###############################################################################
    
    
    ###############################################################################
    ## Create color list                                                         ## 
    sampleColorList <- list()
    
    colorVec <- paramList[["colorPlotsBy"]]
    
    Nopt <- apply(obj@meta.data, 2, function(x) length(unique(x)))
    Nopt <- sort(Nopt[Nopt < NcatColorMax], decreasing = F)
    
    
    
    numOptions <- names(obj@meta.data)[unlist(lapply(obj@meta.data, is.numeric))]
    numOptions <- c(
      "lg10Expr",
      numOptions
    )
    
    colorVec <- colorVec[!(colorVec %in% numOptions)]
    
    
    for (i in 1:length(colorVec)){
        tag <-colorVec[i]
        
        
        sampleVec <- as.vector(sort(unique(obj@meta.data[,tag])))
        sampleVec <- na.omit(sampleVec)
        sampleVec <- sampleVec[sampleVec != ""]
        
        if (length(sampleVec) == 1){
          sampleColVec <- "black"
        }  else if (length(sampleVec) == 2){
          l1 <- length(obj@meta.data[obj@meta.data[,tag] == sampleVec[1],tag])
          l2 <- length(obj@meta.data[obj@meta.data[,tag] == sampleVec[2],tag])
          if (l1 > l2){
            sampleColVec <- c("black", "red")
          } else {
            sampleColVec <- c("red", "black")
          }
          
        } else {
          library(scales)
          sampleColVec <- scales::hue_pal()(length(sampleVec))
          
        }
        
        names(sampleColVec) <- sampleVec
        sampleColorList[[colorVec[i]]] <-  sampleColVec
      }
      paramList[["catColorList"]] <- sampleColorList
      ## Done creating color list
      
      return(paramList)
    
  }
  
)

###############################################################################


###############################################################################
## Create parameter file                                                     ##


#' @title writeAppParameterFiles
#'
#'
#' @param project_id Project id
#' @param projectPath Path to project
#' @param params biologic parameterList
#' 
#' @import Seurat RMySQL RSQLite biologicSeqTools2
#' @export


writeAppParameterFiles <- function(
    project_id = "testApp",
    projectPath = "./",
    params,
    menuParametersFN = "menuParameters.txt",
    colorParametersFN = "colorParameters.txt"
){
    ###########################################################################
    ## Write menu ParameterFile                                              ##
    pos <- grep("catColorList", names(params))
    
    menuList <- params
    if (length(pos) > 0){
        menuList <- menuList[-pos]
    }
    
    mList <- list()
    for (i in 1:length(menuList)){
        mList[[i]] <- rbind(data.frame(
          menuName = rep(names(menuList)[i], length(menuList[[i]])),
          colOption = menuList[[i]],
          colOptionName = menuList[[i]],
          colSel = menuList[[i]],
          displayOrder = 1:length(menuList[[i]])
        ))
        
        
    }
    
    
    dfM <- data.frame(do.call(rbind,mList), stringsAsFactors = F)
    row.names(dfM) <- NULL
    
    dfM[["displayName"]] <- dfM$menuName
    dfM$displayName <- gsub("x_axis", "Choose a X-axis", dfM$displayName)
    dfM$displayName <- gsub("y_axis", "Choose a Y-axis", dfM$displayName)
    dfM$displayName <- gsub("splitPlotsBy", "Split Plots By", dfM$displayName)
    dfM$displayName <- gsub("colorPlotsBy", "Color Plots By", dfM$displayName)
    
    dfM <- dfM[, c("menuName", "displayName", "colOption", "colOptionName" , "displayOrder", "displayName")]
    
    outDir <- paste0(
        projectPath,
        project_id, "_app/parameters"
    )
    
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
    }
    
    FNout <- paste0(
      outDir, "/",
      menuParametersFN
    )
    
    
    write.table(
        dfM, 
        FNout, 
        row.names = F, 
        sep = "\t"
    )
    
    print(paste0("Saved Menuoption parameter file in ",FNout,"."))
    
    ## Done writing menu parameterfile                                       ##
    ###########################################################################
    
    ###########################################################################
    ## save sample color list                                                ##
    
    
    pos <- grep("catColorList", names(params))
    # 
    # menuList <- params
    # if (length(pos) > 0){
    #   menuList <- menuList[-pos]
    # }
    
    #dfCol <- 
    
    colorList <- params[[pos]]
    
    mList <- list()
    for (i in 1:length(colorList)){
      mList[[i]] <- rbind(data.frame(
        menuName = rep(names(colorList)[i], length(colorList[[i]])),
        colOption = names(colorList[[i]]),
        colOptionName = gsub("_", " ", names(colorList[[i]])),
        colSel = as.vector(colorList[[i]]),
        displayOrder = 1:length(colorList[[i]])
      ))
      
      
    }
    
    dfC <- data.frame(do.call(rbind,mList), stringsAsFactors = F)
    row.names(dfC) <- NULL
    
    dfC <- dfC[, c("menuName", "colOption", "colOptionName", "colSel", "displayOrder")]
    
    menuList <- params
    
    dfCol <- data.frame(
        columnName = names(unlist(menuList[[pos]])), 
        colOption = as.vector(unlist(menuList[[pos]]))
    )
    
    outDir <- paste0(
      projectPath,
      project_id, "_app/parameters"
    )
    
    if (!dir.exists(outDir)){
      dir.create(outDir, recursive = T)
    }
    
    FNoutC <- paste0(
      outDir, "/",
      colorParametersFN
    )
    
    
    write.table(
      dfC, 
      FNoutC, 
      row.names = F, 
      sep = "\t"
    )
    
    print(paste0("Saved color options file in ",FNoutC ,"."))
    ## done                                                                  ##
    ###########################################################################
}



## Done                                                                      ##
###############################################################################

###############################################################################
## seuratObjectToLocalViewer                                                 ##

#' @title seuratObjectToLocalViewer
#'
#'
#' @param project_id Project id
#' @param OsC Seurat object
#' @import Seurat RMySQL RSQLite dplyr
#' @export

seuratObjectToLocalViewer <- function(
    params = NULL,
    project_id = "testApp",
    projectPath = "./",
    OsC = NULL,
    dataMode = "SQLite",
    geneDefault = NULL,
    dfExpr = NULL
    #host = host,
    #user = db.user,
    #password = db.pwd
){  
  
    ###############################################################################
    ## Ensure no column is factor                                                ##
    
    for (i in 1:ncol(OsC@meta.data)){
      
      if (is.factor(OsC@meta.data[,i])){
        OsC@meta.data[,i] <- as.character(OsC@meta.data[,i])
        print(paste0("Metadata column ", names(OsC@meta.data)[i], " changed from factor to character."))
      }
      
      # print(is.factor(OsC@meta.data[,i]))
    }
    
    ##
    ###############################################################################
  
    ###############################################################################
    ## Set gene default if it isnt                                               ##
    if (is.null(geneDefault)){
        DefaultAssay(OsC) <- "RNA"
        my_genes <- rownames(x = OsC@assays$RNA)
        
        
        ## Based on https://github.com/satijalab/seurat/issues/3560 the next two lines were 
        # added/altered:
        # in large datasets fetch data produces errors. After some trial and error, 
        # breaking down the problem into chunks seems to solve the problem. 
        cells <- row.names(OsC@meta.data)
        cellList <- split(cells, ceiling(seq_along(cells)/50000))
        
        for (i in 1:length(cellList)){
          expTemp <- FetchData(OsC, vars = my_genes, cells = cellList[[i]] )
          if (i == 1){
            exp <- expTemp
          } else {
            exp <- rbind(exp, expTemp[,colnames(exp)])
          }
          #print(i)
        }
        
        #exp <- FetchData(OsC, vars = my_genes, cells = cells )
        ## End change 2022 03 21
        
        
        ExprMatrix <- round(as.matrix(colMeans(exp  > 0)) *100,1)
        colnames(ExprMatrix)[1] <- "count_cut_off"
        dfExprMatrix <- data.frame(ExprMatrix)
        dfExprMatrix[["gene"]] <- row.names(dfExprMatrix)
        
        
        
        dfPercCellsExpr <- dfExprMatrix
        #dfPercCellsExpr <- dfPercCellsExpr[dfPercCellsExpr$gene %in% Obio@dataTableList$referenceList$integrated_top30var, ]
        dfPercCellsExpr <- dfPercCellsExpr[order(dfPercCellsExpr$count_cut_off, decreasing = T),]
        
        geneDefault <- as.vector(dfPercCellsExpr[1,"gene"])
    }
    ## Done 
    ###############################################################################
    
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
    
    if (is.null(dfExpr)){
        dfExpr <-  biologicViewerSC::createDfExpr(
          obj = OsC,
          assay = "RNA",
          #slot = "data",
          geneSel = NULL
        ) 
    }
    
    ## In Sqlite version load via function ##
    
    
    
    ## Adding gene_ids 
    
    library("dplyr")
    
    dfIDTable <- dfExpr %>% 
      dplyr::select(gene) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(gene_id = dplyr::row_number())
    
    
    ## Upload expression table to database 
    
    print(paste0("Database to be used: ", primDataDB))
    print(paste0("Database table name to be used: ", paste0(project_id, "_geneID_tb")))
    
    biologicSeqTools2::upload.datatable.to.database(
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
        dplyr::rename(condition = cellID)  %>%  
        dplyr::mutate(lg10Expr = round(lg10Expr, 3)) %>% 
        dplyr::arrange(gene) 
    
    ## Upload expression table to database 
    
    print(paste0("Database to be used: ", primDataDB))
    print(paste0("Database table name to be used: ", expDbTable))
    
    colCatList <- biologicSeqTools2::inferDBcategories(dfExpr)
    
    biologicSeqTools2::upload.datatable.to.database(
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
    
    
    dfCoord <- biologicViewerSC::createDfCoord(OsC)
    
    
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
    
    
    
    columnDBcategoryList <- biologicSeqTools2::inferDBcategories(dfData=dfdbTable)
    
    
    biologicSeqTools2::upload.datatable.to.database(
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
    biologicSeqTools2::killDbConnections()
    ##                                                                           ##
    ###############################################################################
    
    ###############################################################################
    ## Create params if they don't exist                                         ##
    if (is.null(params)){
      params <- scanObjParams(OsC)
    }
    ##                                                                           ##
    ###############################################################################
    
    ###############################################################################
    ## Select default splitBy options                                            ##
    
    
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
    
    
    if (!(exists("paramList"))){
        paramList <- list()
    }
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
    
    colorDisplayOptions <- c(
      "lg10Expr",
      colorByOptions
    )
    
    # colorByOptionsPart1 <- c(
    #   grep("clusterName", colorByOptions),
    #   grep("seurat_clusters", colorByOptions),
    #   grep("sampleName", colorByOptions)
    # )
    # 
    # if (length(colorByOptionsPart1) > 0){
    #   colorByOptionsPart2 <- colorByOptions[-colorByOptionsPart1]
    #   colorDisplayOptions <- c(
    #     "lg10Expr",
    #     colorByOptions[colorByOptionsPart1],
    #     colorByOptionsPart2
    #   )
    # } else {
    #   colorDisplayOptions <- c(
    #     "lg10Expr",
    #     colorByOptions[colorByOptionsPart1]
    #   )
    # }
    
    
    
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
    ## Create sample order and color specification files                         ##
    
    writeAppParameterFiles(
        project_id = project_id,
        projectPath = projectPath,
        params = params,
        menuParametersFN = "menuParameters.txt",
        colorParametersFN = "colorParameters.txt"
    )
    
    ## Done 
    ###########################################################################
    
    
    ###########################################################################
    ## Temp section                                                          ##
    
    copyVec <- paste0(
      system.file("scAppTemplate/",package="biologicViewerSC"),
      "/",
      list.files(system.file("scAppTemplate/",package="biologicViewerSC"))
    )
    
    a <- file.copy(copyVec, shinyProjectPath,  recursive=TRUE)
    
    ##                                                                       ##
    ###########################################################################
    
    ###########################################################################
    ##                                                                       ##
    
    dfID <- data.frame(
        dataMode = dataMode,
        url = "",
        id = "",
        id2 = "",
        db =  gsub(paste0("./",project_id, "_app/"),"",primDataDB),
        coordTb = PCAdbTableName,
        exprTb = expDbTable,
        geneTb = geneTb,
        default = geneDefault
    )
    
    FN <- paste0(shinyDataPath, "db.txt")
    write.table(dfID, FN, row.names = F, sep="\t")
    
    
    ##
    ################################
    
    pFolder <- paste0(getwd(), gsub("./", "/", shinyProjectPath))
    
    print(paste0("Local App generated in folder ", pFolder))
  
  
}

## End                                                                       ##
###############################################################################

###############################################################################
## seuratObjectToLocalViewer                                                 ##

#' @title seuratObjectToLocalViewer
#'
#'
#' @param project_id Project id
#' @param OsC Seurat object
#' 
#' @import Seurat RMySQL RSQLite
#' @export

seuratObjectToViewer <- function(
    params = NULL,
    project_id = "testApp",
    projectPath = "./",
    OsC = NULL,
    dataMode = "MySQL",
    host = "dbHostURL",
    dbname = "dbname_db",
    db.pwd = "dbAdminPassword",
    db.user = "boeings",
    appDomains = c("shiny-bioinformatics.crick.ac.uk","10.%"),
    geneDefault = NULL,
    dfExpr = NULL
){  
  ###############################################################################
  ## Ensure no column is factor                                                ##
  
  for (i in 1:ncol(OsC@meta.data)){
    if (is.factor(OsC@meta.data[,i])){
      OsC@meta.data[,i] <- as.character(OsC@meta.data[,i])
      print(paste0("Metadata column ", names(OsC@meta.data)[i], " changed from factor to character."))
    }
    
    # print(is.factor(OsC@meta.data[,i]))
  }
  
  ##
  ###############################################################################
  
  ###############################################################################
  ## Set gene default if it isnt                                               ##
  if (is.null(geneDefault)){
    DefaultAssay(OsC) <- "RNA"
    my_genes <- rownames(x = OsC@assays$RNA)
    
    ## Based on https://github.com/satijalab/seurat/issues/3560 the next two lines were 
    # added/altered:
    cells <- as.factor(row.names(OsC@meta.data))
    
    exp <- FetchData(OsC, my_genes, cells = cells )
    # exp <- FetchData(OsC, my_genes)
    ## End change 2022 03 21
    
    ExprMatrix <- round(as.matrix(colMeans(exp  > 0)) *100,1)
    colnames(ExprMatrix)[1] <- "count_cut_off"
    dfExprMatrix <- data.frame(ExprMatrix)
    dfExprMatrix[["gene"]] <- row.names(dfExprMatrix)
    
    
    
    dfPercCellsExpr <- dfExprMatrix
    dfPercCellsExpr <- dfPercCellsExpr[dfPercCellsExpr$gene %in% Obio@dataTableList$referenceList$integrated_top30var, ]
    dfPercCellsExpr <- dfPercCellsExpr[order(dfPercCellsExpr$count_cut_off, decreasing = T),]
    
    geneDefault <- as.vector(dfPercCellsExpr[1,"gene"])
  }
  ## Done 
  ###############################################################################
  
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
  
  print("Formating expression data. This may take a few minutes, particulary with larger datasets...")
  ## Running this function may take a minute or two, depending on the number of cells in your dataset
  
  if (is.null(dfExpr)){
      dfExpr <-  biologicViewerSC::createDfExpr(
        obj = OsC,
        assay = "RNA",
        #slot = "data",
        geneSel = NULL
      ) 
  }
  
  ## In Sqlite version load via function ##
  
  
  
  ## Adding gene_ids 
  
  library("dplyr")
  
  dfIDTable <- dfExpr %>% 
    dplyr::select(gene) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(gene_id = dplyr::row_number())
  
  
  ## Upload expression table to database 
  
  print(paste0("Database to be used: ", dbname))
  print(paste0("Database table name to be used: ", paste0(project_id, "_geneID_tb")))
  
  biologicSeqTools2::upload.datatable.to.database(
    host = host,
    user = db.user,
    password = db.pwd,
    prim.data.db = dbname,
    dbTableName = geneTb,
    df.data = dfIDTable,
    db.col.parameter.list = biologicSeqTools2::inferDBcategories(dfIDTable),
    new.table = T,
    cols2Index = c("gene"),
    mode = dataMode  # Options: "MySQL" and "SQLite"
  )
  
  
  ## Rearrange expression talbe
  ## This step may take a couple of minutes in a large dataset
  dfExpr <- dfExpr %>% 
    dplyr::rename(condition = cellID)  %>%  
    dplyr::mutate(lg10Expr = round(lg10Expr, 3)) %>% 
    dplyr::arrange(gene) 
  
  ## Upload expression table to database 
  
  print(paste0("Database to be used: ", dbname))
  print(paste0("Database table name to be used: ", expDbTable))
  
  colCatList <- biologicSeqTools2::inferDBcategories(dfExpr)
  
  
  if (nrow(dfExpr) > 100000 & dataMode == "MySQL"){
      ## load infile  
      biologicSeqTools2::uploadDbTableInfile(
          host = host,
          user = db.user,
          password = db.pwd,
          prim.data.db = dbname,
          dbTableName = expDbTable,
          df.data = dfExpr,
          db.col.parameter.list = colCatList,
          new.table = TRUE,
          cols2Index = c("gene"),
          indexName = "idx_gene",
          mode = "MySQL",
          tempFileName = "temp.upload.csv"
          
      )
  } else {
      biologicSeqTools2::upload.datatable.to.database(
        host = host,
        user = db.user,
        password = db.pwd,
        prim.data.db = dbname,
        dbTableName = expDbTable,
        df.data = dfExpr,
        db.col.parameter.list = colCatList,
        new.table = T,
        cols2Index = c("gene"),
        #indexName = c("idx_gene_exp"),
        mode = dataMode  # Options: "MySQL" and "SQLite"
      )
  }
  
  
  
  ###############################################################################
  ## Create Metadata table                                                     ##
  print("Rendering Metadata...")
  
  
  dfCoord <- biologicViewerSC::createDfCoord(OsC)
  
  
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
  
  ## In case of factor problems write to file and read back in without factors
  # write.table(
  #     dfdbTable,
  #     "temp.txt",
  #     row.names = F,
  #     sep = "\t"
  # )
  # 
  # dfdbTable <- read.delim(
  #     "temp.txt",
  #     sep = "\t",
  #     stringsAsFactors = F
  # )
  
  if (file.exists("temp.txt")){
    unlink("temp.txt")
  }
  
  dfdbTable[is.na(dfdbTable)] <- ""
  
  
  
  columnDBcategoryList <- biologicSeqTools2::inferDBcategories(dfData=dfdbTable)
  
  
  biologicSeqTools2::upload.datatable.to.database(
    host = host,
    user = db.user,
    password = db.pwd,
    prim.data.db = dbname,
    dbTableName = PCAdbTableName,
    df.data = dfdbTable,
    db.col.parameter.list = columnDBcategoryList,
    new.table = TRUE,
    mode = dataMode
  )
  biologicSeqTools2::killDbConnections()
  ##                                                                           ##
  ###############################################################################
  
  ###############################################################################
  ## Create params if they don't exist                                         ##
  if (is.null(params)){
    params <- scanObjParams(OsC)
  }
  ##                                                                           ##
  ###############################################################################
  
  ###############################################################################
  ## Select default splitBy options                                            ##
  
  
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
  
  
  if (!(exists("params"))){
    paramList <- list()
  } else {
    paramList <- params
  }
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
  
  
  colorDisplayOptions <- c(
      "lg10Expr",
      colorByOptions
  )
  
  # colorByOptionsPart1 <- c(
  #   grep("clusterName", colorByOptions),
  #   grep("seurat_clusters", colorByOptions),
  #   grep("sampleName", colorByOptions)
  # )
  # 
  # 
  # if (length(colorByOptionsPart1) > 0){
  #   colorByOptionsPart2 <- colorByOptions[-colorByOptionsPart1]
  #   colorDisplayOptions <- c(
  #     "lg10Expr",
  #     colorByOptions[colorByOptionsPart1],
  #     colorByOptionsPart2
  #   )
  # } else {
  #   colorDisplayOptions <- c(
  #     "lg10Expr",
  #     colorByOptions[colorByOptionsPart1]
  #   )
  # }
  
  # if (length(colorByOptionsPart1) > 0){
  #   colorByOptionsPart2 <- colorByOptions[-colorByOptionsPart1]
  # }
  # 
  # 
  # 
  # colorDisplayOptions <- c(
  #   "lg10Expr",
  #   colorByOptions[colorByOptionsPart1],
  #   sort(colorByOptionsPart2)
  # )
  
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
  ## Create sample order and color specification files                         ##
  
  writeAppParameterFiles(
    project_id = project_id,
    projectPath = projectPath,
    params = paramList,
    menuParametersFN = "menuParameters.txt",
    colorParametersFN = "colorParameters.txt"
  )
  
  ## Done 
  ###########################################################################
  
  
  ###########################################################################
  ## Temp section                                                          ##
  
  copyVec <- paste0(
    system.file("scAppTemplate/",package="biologicViewerSC"),
    "/",
    list.files(system.file("scAppTemplate/",package="biologicViewerSC"))
  )
  
  a <- file.copy(copyVec, shinyProjectPath,  recursive=TRUE)
  
  ##                                                                       ##
  ###########################################################################
  
  
  
  ###########################################################################
  ## Create app user and credentials                                       ##
  
  biologicSeqTools2::assignDbUsersAndPrivileges(
      accessFilePath = shinyDataPath,
      hostDbUrl = host,
      appUserName = substr(paste0(project_id, "_aUser"), 1, 30),
      domains = appDomains,
      dbname = dbname,
      tables = c("coordTb" = PCAdbTableName,"exprTb" = expDbTable,"geneTb" = geneTb),
      #recreateProjectUser = TRUE,
      dbAdminUser = db.user,
      dbAdminPwd = db.pwd,
      geneDefault = geneDefault
  ) 
  
  ## Done app user and credentials                                         ##
  ###########################################################################
  
  
  
  
  print(paste0("App generated in folder ", shinyProjectPath))
  
  
}

## End                                                                       ##
###############################################################################


