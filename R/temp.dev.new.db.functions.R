###############################################################################
## Create user and assign privileges                                         ##

#' @title assignDbUsersAndPrivileges
#'
#' @description This function allows you to express your love for the superior furry animal.
#' @param agree Do you agree dogs are the best pet? Defaults to TRUE.
#' @keywords database utils
#' @export
#' @import DBI RMySQL
#' @import methods

assignDbUsersAndPrivileges <- function(
    accessFilePath = shinyDataPath,
    hostDbUrl = "10.27.241.82",
    appUserName = substr(paste0(project_id, "_aUser"), 1, 30),
    geneDefault = "SLPI",
    domains = c("shiny-bioinformatics.crick.ac.uk", "10.%"),
    dbname = "prim.data.db",
    tables = c("coordTb" = PCAdbTableName,"exprTb" = expDbTable,"geneTb" = geneTb),
    recreateProjectUser = TRUE,
    dbAdminUser = "boeings",
    dbAdminPwd = "db.pwd"
) {
    
    ############################
    ## Helper function 
    doQuery <- function(
        user = "db.user", 
        password = "db.upload.pwd", 
        host = "host",
        dbname = "primDataDB",
        query = "mysql db query",
        #existingAccessFileName = existingAccessFileName
        resOut = FALSE
    ){
        library(RMySQL)
        dbDB <- dbConnect(
            drv = RMySQL::MySQL(), 
            user = db.user, 
            password = db.pwd, 
            host = host,
            dbname = primDataDB
        ) 
        
        tryCatch(res <- DBI::dbGetQuery(dbDB, query), error = function(c) {
            c$message <- stop(paste0("Error in ", query, "."))
        })
        
        
        DBI::dbDisconnect(dbDB)
        if (resOut){
            return(res)
        }
        
    }
    
    
    if (file.exists(paste0(accessFilePath, "db.txt"))){
        df <- read.delim(paste0(accessFilePath, "db.txt"), header = T, sep="\t", stringsAsFactors = F)
        sPwd <- as.vector(df$id2)
        sUser <- as.vector(df$id)
    } else {
        sUser <- substr(appUserName, 1, 30)
        sPwd <-c(2:9,letters,LETTERS)
        sPwd <- paste(sample(sPwd,8),collapse="")
    }
        
    if (!file.exists(paste0(accessFilePath, "db.txt"))){
       ## Create user in db
        for (k in 1:length(domains)) {
            query0 <- paste0("SELECT User, Host FROM mysql.user WHERE User = '",sUser,"' AND Host = '",domains[k],"';")
            
            res <- doQuery(
                user = dbAdminUser, 
                password = dbAdminPwd, 
                host = hostDbUrl,
                dbname = dbname,
                query = query0,
                #existingAccessFileName = existingAccessFileName
                resOut = TRUE
            )
            
            if (nrow(res) > 0){
                query0a <- paste0("DROP USER '",sUser,"'@'",domains[k],"';")   
                #doQuery(Obio, query = query0a)
                
                doQuery(
                    user = dbAdminUser, 
                    password = dbAdminPwd, 
                    host = hostDbUrl,
                    dbname = dbname,
                    query = query0a,
                    #existingAccessFileName = existingAccessFileName
                    resOut = FALSE
                )
            }
            
            query1 <- paste0(
                "CREATE USER '",
                sUser, 
                "'@'",domains[k],"' IDENTIFIED BY '",
                sPwd,
                "';"
            )
            
            doQuery(
                user = dbAdminUser, 
                password = dbAdminPwd, 
                host = hostDbUrl,
                dbname = dbname,
                query = query1,
                #existingAccessFileName = existingAccessFileName
                resOut = FALSE
            )
        } # End for k-loop user fix
        
        
        ## Make password file 
        ## Create log-in file ##
        
          
        dfID <- data.frame(
            type = "main",
            dataMode = dataMode,
            url = hostDbUrl,
            id = sUser,
            id2 = sPwd,
            db = dbname,
            coordTb = tables["coordTb"],
            exprTb = tables["exprTb"],
            geneTb = tables["geneTb"],
            default = "SILI"
        )
        
        if (!file.exists(accessFilePath)){
            dir.create(accessFilePath, recursive = T)
        }
        
        write.table(dfID, paste0(accessFilePath, "db.txt"), row.names = F, sep="\t")
        ## Done making password file 
    }
    
    ## Done creating users for all domains                                   ##
    ###########################################################################
    
    ###########################################################################
    ## GRANT access to the app user to all relevant tables                   ##
    for (i in 1:length(tables)){
        for (k in 1:length(domains)){
            query7 <- paste0(
                "GRANT SELECT on ", dbname,".",tables[i], " TO '",sUser,"'@'",domains[k],"';"
            )
            
            doQuery(
                user = dbAdminUser, 
                password = dbAdminPwd, 
                host = hostDbUrl,
                dbname = dbname,
                query = query7,
                #existingAccessFileName = existingAccessFileName
                resOut = FALSE
            )
        }
    }
    ## Done                                                                  ##
    ###########################################################################
        
}    
## End of function assignDbUsersAndPrivileges                                ##
###############################################################################

###############################################################################
## Upload datatable infile                                                   ##

#' @title uploadDbTablInfile
#'
#' @description This function allows you to express your love for the superior furry animal.
#' @param agree Do you agree dogs are the best pet? Defaults to TRUE.
#' @keywords dogs
#' @export
#' @import DBI RMySQL
#' @import methods


uploadDbTableInfile <- function(
    data = "dataframe.to.upload",
    tempFileName = "temp.upload.csv",
    host = ,
    uploadUser = ,
    dbname = "",
    dbTableName = ""
){
    doQuery <- function(
        Obio, 
        query
    ){
        library(RMySQL)
        dbDB <- dbConnect(
            drv = RMySQL::MySQL(), 
            user = Obio@dbDetailList$db.user, 
            password = db.pwd, 
            host = Obio@dbDetailList$host,
            dbname = Obio@dbDetailList$primDataDB
        ) 
        
        tryCatch(dbGetQuery(dbDB, query), error = function(c) {
            c$message <- stop(paste0("Error in ", query, "."))
        })
        
        
        
        
        
        dbDisconnect(dbDB)
        
    }
    
    ## save table locally for in file upload ##
    
    ## Drop existing table if exists
    
    dbtable <- paste0(Obio@parameterList$project_id, "_gene_expr_tb")
    expDbTable <- dbtable
    
    
    query1 <- paste0("DROP TABLE IF EXISTS ", dbtable, ";\n")
    query <- query1
    
    doQuery(Obio, query = query1)
    
    query2 <- paste0(
        #query,
        "CREATE TABLE IF NOT EXISTS ",
        dbtable,
        " (gene VARCHAR(100) CHARACTER SET latin1 COLLATE latin1_swedish_ci, cellID VARCHAR(100) CHARACTER SET latin1 COLLATE latin1_swedish_ci, lg10Expr DECIMAL(6,3) NULL DEFAULT NULL, row_names INT(10) NOT NULL AUTO_INCREMENT,PRIMARY KEY (row_names)); "
    )
    
    doQuery(Obio, query = query2)
    query <- c(query, query2)
    
    
    ## infile upload
    query3 <- paste0(
        #query,
        "LOAD DATA LOCAL INFILE '",
        Obio@parameterList$workdir,
        "temp.",Obio@parameterList$project_id,".csv",
        "' INTO TABLE ",
        dbtable, 
        " FIELDS TERMINATED BY ',' ENCLOSED BY '\"' LINES TERMINATED BY '\n' IGNORE 1 LINES (gene, cellID, lg10Expr, row_names);"
    )
    
    
    doQuery(Obio, query = query3)
    query <- c(query, query3)
    
    ## alter and index 
    query4 <- paste0(
        #query,
        "ALTER TABLE ", dbtable, " ADD UNIQUE(row_names);"
    )
    
    doQuery(Obio, query = query4)
    query <- c(query, query4)
    
    query5 <- paste0(
        #query,
        "CREATE INDEX idx_gene ON ", dbtable, " (gene);"
    )
    
    doQuery(Obio, query = query5)
    query <- c(query, query5)
    ## Add user managment ##
    
}


#CREATE TABLE IF NOT EXISTS test_expr (`row_names` bigint(10) NOT NULL, `gene_id` INT(5), `gene` VARCHAR(100) CHARACTER SET latin1 COLLATE latin1_swedish_ci, `lg10Expr` DECIMAL(6,3) NULL DEFAULT NULL, PRIMARY KEY (`row_names`));




#LOAD DATA LOCAL INFILE '/camp/stp/babs/working/boeings/Projects//schaefera/tobias.ackels/319_scRNAseq_mm_olfactory_bulb_proj_neuron_analysis_SC19135/workdir/temp.csv' INTO TABLE test_expr FIELDS TERMINATED BY ',' ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 LINES (id, mycol1, mycol2);




## Done                                                                      ##
###############################################################################

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
    geneDefault = NULL
    
    
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
    
    
    primDataDB = dbname
    
    
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
    
    dfExpr <-  biologicViewerSC::createDfExpr(
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
        host = host,
        user = db.user,
        password = db.pwd,
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
        host = host,
        user = db.user,
        password = db.pwd,
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
    ## Create app user and credentials                                       ##
    
    assignDbUsersAndPrivileges(
        accessFilePath = shinyDataPath,
        hostDbUrl = host,
        appUserName = substr(paste0(project_id, "_aUser"), 1, 30),
        geneDefault = geneDefault,
        domains = appDomains,
        dbname = dbname,
        tables = c("coordTb" = PCAdbTableName,"exprTb" = expDbTable,"geneTb" = geneTb),
        #recreateProjectUser = TRUE,
        dbAdminUser = db.user,
        dbAdminPwd = db.pwd
    ) 
    
    ## Done app user and credentials                                         ##
    ###########################################################################
    
    
    pFolder <- paste0(getwd(), gsub("./", "/", shinyProjectPath))
    
    print(paste0("App generated in folder ", pFolder))
    
    
}

## End                                                                       ##
####################################################
###########################
