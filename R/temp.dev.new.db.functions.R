

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
    host = "",
    uploadUser = "",
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

