###############################################################################
# (19) upload.datatable.to.database()
###############################################################################

## Indexing of gene name column
# CREATE INDEX idx_mgi_symbol ON p268_rna_seq_table(mgi_symbol)
#' @title upload.datatable.to.database
#'
#' @param host URL or IP address of the database server. NULL if mode is SQLite
#' @param user database user name. Needs privileges for INSERT, DELETE, DROP and SELECT
#' @param password database password
#' @param prim.data.db primary database name
#' @param dbTableName Name of the database table to upload
#' @param df.data data.frame to upload to the database
#' @param db.col.parameter.list This list specifies the category for a database column. This ideally is filled using the function inferDbColumns(dbTableName)
#' @param increment = 5000,
#' @param new.table = FALSE,
#' @param first.row.name.index = 1,
#' @param startOnlyWithConnectionCount1 = FALSE,
#' @param cols2Index = NULL,
#' @param mode = "MySQL"
#'
#' @description mode options: SQLite, MySQL
#' @keywords MySQL SQLite Upload Database
#' @import RMySQL RSQLite DBI
#' @export
#'
#'

# 2021-08-1 Added sqlite options #

upload.datatable.to.database <- function(
    host = NULL,
    user = NULL,
    password = NULL,
    prim.data.db = "project.database",
    dbTableName = "rnaseqdbTableName",
    df.data = "df.data.to.upload",
    db.col.parameter.list = list(
        "VARCHAR(255) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("gene_description"),
        "VARCHAR(50) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("ENSG", "ENSMUSG", "hgnc_symbol", "mgi_symbol", "uniprot", "entrezgene","display_ptm", "^sequence_window", "p_site_env","for_GSEA_gene_chip","associated_gene_name", "gene_type"),
        "VARCHAR(1) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("ppos", "amino_acid", "charge","known_site"),
        "BIGINT(8) NULL DEFAULT NULL" = c("row_names"),
        "INT(8) NULL DEFAULT NULL" = c("row_id", "cluster_order","cluster_id", "count_cut_off", "^position$", "raw_counts"),
        "DECIMAL(6,3) NULL DEFAULT NULL" = c("norm_counts", "NES", "logFC", "lg2_avg", "intensity", "^int", "iBAQ","^localization_prob$"),
        "DECIMAL(6,5) NULL DEFAULT NULL" = c("padj", "pvalue","^pep$")
    ),
    increment = 5000,
    new.table = FALSE,
    first.row.name.index = 1,
    startOnlyWithConnectionCount1 = FALSE,
    cols2Index = NULL,
    indexName = NULL,
    mode = "MySQL",
    addRowNamesColumn = TRUE
){
    
    if (sum( nchar(names(df.data)) > 64) > 0){
        print("Table names clipped to 64 characters.")
        names(df.data) <- substr(names(df.data), 1, 64)
    }
    
    
    if (startOnlyWithConnectionCount1){
        
        ## helper function ##
        getConnectonCount <- function(
            user= "user",
            password = "password",
            dbname = "prim.data.db",
            host = "host"){
            
            if (mode == "SQLite"){
                
                dbDB <- DBI::dbConnect(
                    drv = RSQLite::SQLite()
                )
                
            } else {
                
                dbDB <- DBI::dbConnect(
                    drv = RMySQL::MySQL(),
                    user = user,
                    password = password,
                    host = host
                    
                )
                
            }
            
            
            
            connectionCount <- as.numeric(
                DBI::dbGetQuery(
                    dbDB,
                    "SELECT COUNT(1) ConnectionCount, user FROM information_schema.processlist WHERE user <> 'system user' AND user = 'boeingS' GROUP BY user;"
                )$ConnectionCount
            )
            
            DBI::dbDisconnect(dbDB)
            
            return(connectionCount)
        }
        
        connectionCount <- getConnectonCount(
            user= user,
            password = password,
            dbname = prim.data.db,
            host = host
        )
        
        while (connectionCount > 2){
            print(paste(connectionCount, "connections open. Sleep 30 seconds and try again."))
            
            Sys.sleep(30)
            
            
            connectionCount <- getConnectonCount(
                user= user,
                password = password,
                dbname = prim.data.db,
                host = host
            )
            
            
        }
        
    }
    
    ## Uploading of data frame to database. Happens only if all columns are defined ##
    # library(RMySQL)
    ## Connect to MySQL to check existence of database ##
    if (mode == "SQLite"){
        
        # create path to database. To do that, split everything that's not the database file name. 
        pos <- grep("/", prim.data.db)
        if (length(pos) > 0){
            prim.data.dbPath <- unlist(strsplit(prim.data.db, "/"))[1:(length(unlist(strsplit(prim.data.db, "/")))-1)]
            prim.data.dbPath <- paste0(prim.data.dbPath, collapse = "/")    
        } else {
            prim.data.dbPath <- "./"
        }
        
        
        if (!dir.exists(prim.data.dbPath)){
            dir.create(prim.data.dbPath, recursive = T)
        }
        
        
        dbDB <- DBI::dbConnect(
            drv = RSQLite::SQLite(),
            dbname=paste0(prim.data.db)
        )
        
        
        
    } else {
        
        dbDB <- DBI::dbConnect(
            drv = RMySQL::MySQL(),
            user = user,
            password = password,
            host = host,
            #dbname=prim.data.db,
            new.table = TRUE
        )
        
        ## Create the database if it does not exist already##
        res <- DBI::dbGetQuery(
            dbDB,
            paste(
                "CREATE DATABASE IF NOT EXISTS ",
                prim.data.db,
                sep = ""
            )
        )
        
    }
    
    
    
    
    
    RMySQL::dbDisconnect(dbDB)
    
    ## Ensure that df.data has a row_names column ##
    if (addRowNamesColumn){
        df.data[["row_names"]] <- first.row.name.index:(first.row.name.index+nrow(df.data)-1)
    }
    
    ## Check if all columns are assigned in db.col.parameter.list ##
    all.col.string.vec <- as.vector(do.call('c', db.col.parameter.list))
    
    ## Create a vector that contains all col names that contain at least in part the string in all.cols.vec
    
    ###############################################################################
    ## Function start                                                            ##
    get.all.col.names.with.these.strings <- function(all.col.string.vec){
        all.assigned.cols <- vector(mode="character", length=0)
        for (i in 1:length(all.col.string.vec)){
            pos <- grep(all.col.string.vec[i], names(df.data))
            if (length(pos) > 0){
                all.assigned.cols <- append(all.assigned.cols, names(df.data)[pos])
            }
        }
        return(all.assigned.cols)
    }
    
    ## End of function                                                           ##
    ###############################################################################
    all.assigned.cols <- get.all.col.names.with.these.strings(all.col.string.vec)
    
    ## Ensure that all database columns are assigned ##
    not.assigned <- names(df.data)[!(names(df.data) %in% all.assigned.cols)]
    
    if (length(not.assigned) == 0){
        print("All database columns are defined. Uploading to database...")
    } else {
        print(
            paste0(
                "The following database columns have not been defined: ",
                paste(not.assigned,
                      collapse = ', '
                ),
                ". Datatable completed.")
        )
        stop(not.assigned)
    }
    
    
    
    ## Establish connection ##
    if (mode == "SQLite"){
        
        dbDB <- DBI::dbConnect(
            drv = RSQLite::SQLite(),
            dbname=paste0(prim.data.db)
        )
        
    } else {
        
        dbDB <- DBI::dbConnect(
            drv = RMySQL::MySQL(),
            user = user,
            password = password,
            host = host,
            dbname=prim.data.db
        )
        
    }
    
    
    
    ## Remove all tables with the same name from db ##
    if (new.table){
        res <- DBI::dbExecute(
            dbDB,
            paste(
                "DROP TABLE IF EXISTS ",
                dbTableName,
                sep = ""
            )
        )
        RMySQL::dbDisconnect(dbDB)
    }
    
    ## Upload up to increment rows in one go ##
    iter <- nrow(df.data)%/%increment
    if (nrow(df.data)%%increment != 0){
        iter <- iter + 1
    }
    
    totalRows <- nrow(df.data)
    
    ## Begin of loop
    #res <- purrr::map(1:iter, function(i){
    for (i in 1:iter){
        if (nrow(df.data) > increment){
            limit <- increment
        } else {
            limit <- nrow(df.data)
        }
        df.temp <- df.data[1:limit,]
        df.data <- df.data[(increment+1):nrow(df.data),]
        
        uploaded = FALSE
        while (!uploaded){
            tryCatch({
                biologicSeqTools2::killDbConnections()
                
                ## Establish connection ##
                if (mode == "SQLite"){
                    
                    dbDB <- DBI::dbConnect(
                        drv = RSQLite::SQLite(),
                        dbname=paste0(prim.data.db)
                    )
                    
                } else {
                    
                    dbDB <- DBI::dbConnect(
                        drv = RMySQL::MySQL(),
                        user = user,
                        password = password,
                        host = host,
                        dbname=prim.data.db
                    )
                    
                }
                
                
                ## Upload new dataframe to database ##
                res <- DBI::dbWriteTable(
                    dbDB,
                    dbTableName,
                    df.temp,
                    row.names = FALSE,
                    append = TRUE,
                    overwrite = FALSE
                )
                
                RMySQL::dbDisconnect(dbDB)
                uploaded = TRUE
                #dbDisconnect(dbDB)
            }, error=function(e){cat("Upload errror :",conditionMessage(e), "\n")})
        }
        
        print(paste0(i * increment, " rows out of ",totalRows," processed..."))
        ## Connect to database for dbtable upload  ##
    }  #End of upload loop
    
    print("Processing successfully completed")
    ####################################################
    ## Function alterDBtable
    alterDBtable <- function(
        cmd.string = "mysql command",
        user = "user",
        password = "password",
        dbname = "prim.data.db",
        host = "host"
    ){
        ## Establish connection ##
        if (mode == "SQLite"){
            
            dbDB <- DBI::dbConnect(
                drv = RSQLite::SQLite(),
                dbname=paste0(prim.data.db)
            )
            
        } else {
            
            dbDB <- DBI::dbConnect(
                drv = RMySQL::MySQL(),
                user = user,
                password = password,
                host = host,
                dbname=prim.data.db
            )
            
        }
        
        
        tryCatch({
            DBI::dbExecute(
                dbDB,
                cmd.string
            )}, error=function(e) {paste0("Alter not executed. cmd.vector[", i, "]")})
        
        RMySQL::dbDisconnect(dbDB)
        
        
    }
    
    ## End of function ##
    ######################
    mysql.cmd = ""
    if (new.table){
        alterDBtable(
            cmd.string = paste(
                "ALTER TABLE `",
                dbTableName,
                "` ADD UNIQUE(`row_names`)",
                sep = ""
            ),
            user = user,
            password = password,
            dbname = dbname,
            host = host
        )
        
        ## Describe key columns in database table ##
        mysql.cmd <- paste(
            "ALTER TABLE `",
            dbTableName,
            "` ADD UNIQUE(`row_names`)",
            sep = ""
        )
        
        
        alterDBtable(
            cmd.string = paste(
                "ALTER TABLE `",
                dbTableName,
                "` ADD PRIMARY KEY(`row_names`)",
                sep = ""
            ),
            user = user,
            password = password,
            dbname = dbname,
            host = host
        )
        
        
        
        
        ###############################################################################
        ## Characterize and define secondary database columns                        ##
        for (i in 1:length(db.col.parameter.list)) {
            descriptor <- names(db.col.parameter.list[i])
            cols.in.class <-
                get.all.col.names.with.these.strings(db.col.parameter.list[[i]])
            
            if (length(cols.in.class) > 0) {
                print(
                    paste0(
                        "Assigned ",
                        paste0(cols.in.class, collapse = ', '),
                        " as ",
                        descriptor, "."
                    )
                )
                
                ## Assign column names to MySQL class ##
                alteration.string <-
                    paste0("ALTER TABLE ", dbTableName, " ")
                for (j in 1:length(cols.in.class)) {
                    alteration.string <- paste0(
                        alteration.string,
                        paste0(
                            "CHANGE `", cols.in.class[j], "` `", cols.in.class[j], "` ", descriptor, ", "
                        )
                    )
                }
                ## Remove last comma from string
                alteration.string <-
                    substr(alteration.string, 1, (nchar(alteration.string) - 2))
                
                ## Carry out alteration
                ## Connect to database for dbtable upload  ##
                ## Connection is repated to avoid loss of a short lived connection.
                alterDBtable(
                    cmd.string = alteration.string,
                    user = user,
                    password = password,
                    dbname = dbname,
                    host = host
                )
            }
            #print(alteration.string)
            mysql.cmd <- append(mysql.cmd,
                                alteration.string)
            
            #dbGetQuery(dbDB,
            #           alteration.string)
            
            
        }
        
    }
    ## End characterize and define secondary database columns                  ##
    #############################################################################
    ## Add index based on row namems ##
    
    if (length(cols2Index) > 0){
        for (i in 1:length(cols2Index)){
            print("...indexing...")
            
            
            ## Get all existing indixes
            if (mode == "SQLite"){
                
                dbDB <- DBI::dbConnect(
                    drv = RSQLite::SQLite(),
                    dbname=prim.data.db
                )
                
                ## get all existing indeces
                listCmd <- paste0(
                    "SELECT name FROM sqlite_master WHERE type = 'index';"
                )
                
                indexVec <- as.vector(DBI::dbGetQuery(
                    dbDB,
                    listCmd
                )[,1])
                
            }
            
            if (!is.null(indexName)){
                if (!is.na(indexName[i])){
                    indexCmdName <- indexName[i]
                } else {
                    indexCmdName <- paste0("idx_",cols2Index[i])
                }
            } else {
                indexCmdName <- paste0("idx_",cols2Index[i])
            }
            
            
            if (mode == "SQLite"){
                while (sum(indexCmdName %in% indexVec) > 0){
                    indexCmdName <- paste0(indexCmdName, sample(9, 1))
                }
            }
            
            cmd.string <- paste0("CREATE INDEX ",indexCmdName," ON ",dbTableName," (",cols2Index[i],")")
            
            
            ## Establish connection ##
            if (mode == "SQLite"){
                
                dbDB <- DBI::dbConnect(
                    drv = RSQLite::SQLite(),
                    dbname=prim.data.db
                )
                
                ## get all existing indeces
                listCmd <- paste0(
                    "SELECT name FROM sqlite_master WHERE type = 'index';"
                )
                
                indexVec <- DBI::dbGetQuery(
                    dbDB,
                    listCmd
                )
                
            } else {
                
                dbDB <- DBI::dbConnect(
                    drv = RMySQL::MySQL(),
                    user = user,
                    password = password,
                    host = host,
                    dbname=prim.data.db
                )
                
            }
            
            
            
            tryCatch({
                DBI::dbExecute(
                    dbDB,
                    cmd.string
                )}, error=function(e) {stop(paste0("Database table not uploaded. Problem adding index ",cols2Index[i],"."))})
            
            DBI::dbDisconnect(dbDB)
            
            print(paste0("Datatable ", dbTableName, " successfully uploaded and column(s) ",paste(cols2Index, collapse = " ")," indexed."))
        }
    }
    return(mysql.cmd)
}



## End of function                                                           ##
###############################################################################

#' @title inferDBcategories
#'
#' @description Method description
#' @param agree TBD
#' @keywords TBD
#' @export
#'
###############################################################################
## Function inferDBcategories                                                ##

inferDBcategories <- function(
    dfData
){
    dbCatList <- list()
    
    for (i in 1:ncol(dfData)){
        classLabel <- ""
        maxStringLength <- max(nchar(as.character(dfData[,i])), na.rm = T) + 2
        
        if (is.numeric(dfData[,i])){
            if (is.integer(dfData[,i])){
                if (maxStringLength <= 8) {
                    classLabel <- "INT(8) NULL DEFAULT NULL"
                } else {
                    classLabel <- "BIGINT(8) NULL DEFAULT NULL"
                }
            } else {
                if (max(dfData[,i], na.rm = T) <= 1){
                    classLabel <- "DECIMAL(6,5) NULL DEFAULT NULL"
                } else if (max(dfData[,i], na.rm = T) <= 1000){
                    classLabel <- "DECIMAL(6,3) NULL DEFAULT NULL"
                } else if (max(dfData[,i], na.rm = T) <= 10000){
                    classLabel <- "DECIMAL(6,1) NULL DEFAULT NULL"
                } else {
                    classLabel <- "DECIMAL(8,0) NULL DEFAULT NULL"
                }
            }
        } else {
            ## Running as character
            classLabel <- "VARCHAR(255) CHARACTER SET utf8 COLLATE utf8_general_ci"
            if (maxStringLength < 100){
                classLabel <- "VARCHAR(100) CHARACTER SET utf8 COLLATE utf8_general_ci"
            } else if (maxStringLength < 50){
                classLabel <- "VARCHAR(50) CHARACTER SET utf8 COLLATE utf8_general_ci"
            } else if (maxStringLength < 10){
                classLabel <- "VARCHAR(10) CHARACTER SET utf8 COLLATE utf8_general_ci"
            }
        }
        
        pos <- grep(classLabel, names(dbCatList), fixed=TRUE)
        if (length(pos) > 0 ){
            dbCatList[[classLabel]] <- c(dbCatList[[classLabel]], paste0("^", names(dfData)[i], "$"))
        } else {
            dbCatList[[classLabel]] <- paste0("^", names(dfData)[i], "$")
        }
        
    }
    
    
    
    ## Make sure row_names is prsent ##
    classLabel <- "BIGINT(8) NULL DEFAULT NULL"
    pos <- grep(classLabel, names(dbCatList), fixed=TRUE)
    if (length(pos) > 0){
        dbCatList[[classLabel]] <- c(dbCatList[[classLabel]], paste0("^row_names$"))
    } else {
        dbCatList[[classLabel]] <- paste0("^row_names$")
    }
    
    return(dbCatList)
}

## EOF inferDB category                                                      ##
###############################################################################

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

# Maximum length for maria db username is 16 #

assignDbUsersAndPrivileges <- function(
    accessFilePath = shinyDataPath,
    hostDbUrl = "10.27.241.82",
    appUserName = substr(paste0(project_id, "_aUser"), 1, 15),
    geneDefault = NULL,
    domains = c("shiny-bioinformatics.crick.ac.uk", "10.%"),
    dbname = "prim.data.db",
    tables = c("coordTb" = PCAdbTableName,"exprTb" = expDbTable,"geneTb" = geneTb),
    recreateProjectUser = TRUE,
    dbAdminUser = "boeings",
    dbAdminPwd = "db.pwd",
    dataMode = "MySQL"
) {
    
    ## Maximum length for app user name is 16
    appUserName <- substr(appUserName, 1, 15)
    
    ############################
    ## Helper function
    doQuery <- function(
        user = "db.user",
        password = "db.upload.pwd",
        host = "host",
        dbname = "primDataDB",
        query = "mysql db query",
        #existingAccessFileName = existingAccessFileName
        resOut = FALSE,
        geneDefault = NULL
    ){
        library(RMySQL)
        dbDB <- dbConnect(
            drv = RMySQL::MySQL(),
            user = user,
            password = password,
            host = host,
            dbname = dbname
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
            default = geneDefault
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
## Kill connections                                                          ##


#' @title A method
#'
#' @description Method description
#' @param agree TBD
#' @keywords TBD
#' @export
#'
#'
killDbConnections <- function () {
    library(RMySQL)
    all_cons <- dbListConnections(MySQL())
    print(all_cons)
    for(con in all_cons)
        res <- dbDisconnect(con)
    #print(paste(length(all_cons), " connections killed."))
}
##                                                                           ##
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
    host = NULL,
    user = NULL,
    password = NULL,
    prim.data.db = "project.database",
    dbTableName = "rnaseqdbTableName",
    df.data = "df.data.to.upload",
    db.col.parameter.list = list(
        "VARCHAR(255) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("gene_description"),
        "VARCHAR(50) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("ENSG", "ENSMUSG", "hgnc_symbol", "mgi_symbol", "uniprot", "entrezgene","display_ptm", "^sequence_window", "p_site_env","for_GSEA_gene_chip","associated_gene_name", "gene_type"),
        "VARCHAR(1) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("ppos", "amino_acid", "charge","known_site"),
        "BIGINT(8) NULL DEFAULT NULL" = c("row_names"),
        "INT(8) NULL DEFAULT NULL" = c("row_id", "cluster_order","cluster_id", "count_cut_off", "^position$", "raw_counts"),
        "DECIMAL(6,3) NULL DEFAULT NULL" = c("norm_counts", "NES", "logFC", "lg2_avg", "intensity", "^int", "iBAQ","^localization_prob$"),
        "DECIMAL(6,5) NULL DEFAULT NULL" = c("padj", "pvalue","^pep$")
    ),
    new.table = TRUE,
    cols2Index = NULL,
    indexName = NULL,
    mode = "MySQL",
    tempFileName = "temp.upload.csv"
    
){
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
            user = user,
            password = password,
            host = host,
            dbname = dbname
        )
        
        tryCatch(res <- DBI::dbGetQuery(dbDB, query), error = function(c) {
            c$message <- stop(paste0("Error in ", query, "."))
        })
        
        
        DBI::dbDisconnect(dbDB)
        if (resOut){
            return(res)
        }
        
    }
    
    ###########################################################################
    ## save table locally for in file upload                                 ##
    
    write.csv(df.data, tempFileName, row.names = F)
    rm(df.data)
    
    ##                                                                       ##
    ###########################################################################
    
    ## Drop existing table if exists
    
    query1 <- paste0("DROP TABLE IF EXISTS ", dbTableName, ";\n")
    
    doQuery(
        user = user,
        password = password,
        host = host,
        dbname = prim.data.db,
        query = query1,
        #existingAccessFileName = existingAccessFileName
        resOut = FALSE
    )
    
    
    query2 <- paste0(
        #query,
        "CREATE TABLE IF NOT EXISTS ",
        dbTableName,
        " (gene VARCHAR(100) CHARACTER SET latin1 COLLATE latin1_swedish_ci, cellID VARCHAR(100) CHARACTER SET latin1 COLLATE latin1_swedish_ci, lg10Expr DECIMAL(6,3) NULL DEFAULT NULL, row_names INT(10) NOT NULL AUTO_INCREMENT,PRIMARY KEY (row_names)); "
    )
    
    doQuery(
        user = user,
        password = password,
        host = host,
        dbname = prim.data.db,
        query = query2,
        #existingAccessFileName = existingAccessFileName
        resOut = FALSE
    )
    
    print("Data is being formated. This may take a few minutes for larger datasets.")
    
    ## infile upload
    query3 <- paste0(
        #query,
        "LOAD DATA LOCAL INFILE '",
        tempFileName,
        "' INTO TABLE ",
        dbTableName,
        " FIELDS TERMINATED BY ',' ENCLOSED BY '\"' LINES TERMINATED BY '\n' IGNORE 1 LINES (gene, cellID, lg10Expr, row_names);"
    )
    
    
    doQuery(
        user = user,
        password = password,
        host = host,
        dbname = prim.data.db,
        query = query3,
        #existingAccessFileName = existingAccessFileName
        resOut = FALSE
    )
    
    ## alter and index
    query4 <- paste0(
        #query,
        "ALTER TABLE ", dbTableName, " ADD UNIQUE(row_names);"
    )
    
    doQuery(
        user = user,
        password = password,
        host = host,
        dbname = prim.data.db,
        query = query4,
        #existingAccessFileName = existingAccessFileName
        resOut = FALSE
    )
    
    
    query5 <- paste0(
        #query,
        "CREATE INDEX idx_gene ON ", dbTableName, " (gene);"
    )
    
    doQuery(
        user = user,
        password = password,
        host = host,
        dbname = prim.data.db,
        query = query5,
        #existingAccessFileName = existingAccessFileName
        resOut = FALSE
    )
    unlink(tempFileName)
    print("Data loaded infile.")
}


## Done                                                                      ##
###############################################################################


