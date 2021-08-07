###############################################################################
## Upload Anna's samples                                                     ##


# Load Seurat Object from rds file

rdsFN <- "/camp/lab/pachnisv/working/SingleCellData/Internal/SmallIntestineDevelopmentalData/SmallIntestineDevelopmentalData.rds"
rdsFN <- "/camp/lab/pachnisv/working/SingleCellData/Internal/MouseGliaIfngrKO/MouseGliaIfngrKO.rds"
rdsFN <- "/camp/lab/pachnisv/working/SingleCellData/Internal/EGCsHPoly/EGCsHPoly.rds"
#rdsFN <- "/camp/lab/pachnisv/working/SingleCellData/Internal/InVitroData/InVitroData.rds"
rdsFN <- "/camp/lab/pachnisv/working/SingleCellData/Internal/InVitroData/InVitroFoxd3RAData.rds"

rdsFN <- "/camp/lab/pachnisv/working/Michael/FranzeColouredSeuratObject/ColoredSeuratObject.rds"
## name to be used PubMouseGliaIfngrKO


## Done                                                                      ##
###############################################################################

###############################################################################
##                                                                           ##

if (dir.exists("/Volumes/babs/working/boeings/")){
  hpc.mount <- "/Volumes/babs/working/boeings/"
} else if (dir.exists("Y:/working/boeings/")){
  hpc.mount <- "Y:/working/boeings/"
} else if (dir.exists("/camp/stp/babs/working/boeings/")){
  hpc.mount <- "/camp/stp/babs/working/boeings/"
} else {
  hpc.mount <- ""
}
#Create the environment and load a suitable version of R, e.g. so:
FN <- paste0(hpc.mount, "Projects/reference_data/documentation/BC.parameters.txt")
dbTable <- read.delim(
  FN, 
  sep = "\t",
  stringsAsFactors = F
)
db.pwd <- as.vector(dbTable[1,1])


##                                                                           ##
###############################################################################

###############################################################################
## Set Project Parameters                                                    ##
h <- unlist(strsplit(rdsFN, "\\/"))[length(unlist(strsplit(rdsFN, "\\/")))]
project_id <- gsub(".rds", "", h)
#project_id <- "InVitroFoxd3RAData"
project_id <- "PubMouseGliaIfngrKO"

host = "10.27.241.82"
db.user = "boeingS"
password = db.pwd
primDataDB = "vpl_data"


PCAdbTableName <- paste0(
    project_id, 
    "_PCA"
)


##                                                                           ##
###############################################################################






###############################################################################
##
devtools::install_github("decusInLabore/biologicSeqTools")
library(Seurat)
library(biologicSeqTools)

##
###############################################################################

###############################################################################
## Step 1: Load Seurat object                                                ##

OsC <- readRDS(rdsFN)

##                                                                           ##
###############################################################################

###############################################################################
## Create Expr table                                                         ##
dfExpr <-  createDfExpr(
    obj = OsC,
    assay = "RNA",
    #slot = "data",
    geneSel = NULL
) 




dfIDTable <- dfExpr
dfIDTable[["gene_id"]] <- 0
dfIDTable <- unique(dfIDTable[,c("gene", "gene_id")])
dfIDTable <- dfIDTable[order(dfIDTable$gene, decreasing = F), ]
dfIDTable[["gene_id"]] <- 1:nrow(dfIDTable)
upload.datatable.to.database(
    host = host,
    user = db.user,
    password = db.pwd,
    prim.data.db = primDataDB,
    dbTableName = paste0(project_id, "_geneID_tb"),
    df.data = dfIDTable,
    db.col.parameter.list = list(
        "BIGINT(8) NULL DEFAULT NULL" = c("row_names"),
        "VARCHAR(100) CHARACTER SET utf8 COLLATE utf8_general_ci" = c("gene"),
        "INT(8) NULL DEFAULT NULL" = c("gene_id")
    ),
    new.table = T,
    cols2Index = c("gene")
)


## Part2
names(dfExpr) <- gsub("condition", "cellID", names(dfExpr))
names(dfExpr) <- gsub("expr", "lg10Expr", names(dfExpr))
dfExpr$lg10Expr <- round(dfExpr$lg10Expr, 3)
## Order by gene_id
dfExpr <- dfExpr[order(dfExpr$gene, decreasing = F),]
dfExpr[["row_names"]] <- 1:nrow(dfExpr)
#dfExpr.store <- dfExpr
#dfExpr <- dfExpr[sample(dfExpr$row_names, 10000),]
tempFN <- paste0("temp.",project_id,".csv")
write.csv(dfExpr, tempFN, row.names = F)
#rm(dfExpr)
dbtable <- paste0(project_id, "_gene_expr_tb")
expDbTable <- dbtable
doQuery <- function(
    #Obio, 
    query,
    resOut = FALSE
){
    #library(RMySQL)
    res <- NULL
    dbDB <- DBI::dbConnect(
        drv = RMySQL::MySQL(), 
        user = db.user, 
        password = db.pwd, 
        host = host,
        dbname = primDataDB
    ) 
    tryCatch(res <- dbGetQuery(dbDB, query), error = function(c) {
        c$message <- stop(paste0("Error in ", query, "."))
    })
    
    dbDisconnect(dbDB)
    
    if (resOut){
        return(res)
    }
    
}

## Done expression
###############################################################################

###############################################################################
## Upload Metadata table                                                     ##

###############################################################################
## Create Metadata table                                                     ##
dfCoord <- createDfCoord(OsC)
##                                                                           ##
###############################################################################

dupTest <- duplicated(toupper(names(dfCoord)))

if (sum(dupTest) > 0 ){
  names(dfCoord)[duplicated(toupper(names(dfCoord)))] <- paste0( names(dfCoord)[duplicated(toupper(names(dfCoord)))], "_B")
}

dfCoord <- dfCoord [,!(duplicated(names(dfCoord)))]


names(dfCoord) <- gsub("\\.", "_", names(dfCoord))
names(dfCoord) <- gsub("^umap", "UMAP", names(dfCoord))
names(dfCoord) <- gsub("phase", "Phase", names(dfCoord))
names(dfCoord) <- gsub("time_point", "meta_Timepoint", names(dfCoord))
names(dfCoord) <- gsub("index", "cellID", names(dfCoord))


# Mouse project
if (project_id == "MouseGliaIfngrKO"){
  names(dfCoord)[grep("genotype", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("genotype", names(dfCoord))])
  names(dfCoord)[grep("shortName", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("shortName", names(dfCoord))])
  
  names(dfCoord)[grep("treatment", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("treatment", names(dfCoord))])
  names(dfCoord)[grep("seurat_clusters", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("seurat_clusters", names(dfCoord))])
  names(dfCoord)[grep("Seurat_Clusters_B", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("Seurat_Clusters_B", names(dfCoord))])
  names(dfCoord)[grep("quadrant", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("quadrant", names(dfCoord))])
  names(dfCoord)[grep("longTag", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("longTag", names(dfCoord))])
  names(dfCoord)[grep("isOld", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("isOld", names(dfCoord))])
  names(dfCoord)[grep("integrated_snn_res_0_8", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("integrated_snn_res_0_8", names(dfCoord))])
  
  dfCoord$meta_seurat_clusters <- paste0("Cluster_", dfCoord$meta_seurat_clusters)
  dfCoord$meta_Seurat_Clusters_B <- paste0("Cluster_", dfCoord$meta_Seurat_Clusters_B)
  
  dfdbTable <- dfCoord
  #dfdbTable[["sample_group"]] <- paste0("Cluster_r08_", dfdbTable$louvain_r0_8)
  dfdbTable[["sample_group"]] <- dfdbTable$meta_seurat_clusters
  
}


if (project_id == "PubMouseGliaIfngrKO"){
  names(dfCoord)[grep("genotype", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("genotype", names(dfCoord))])
  names(dfCoord)[grep("shortName", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("shortName", names(dfCoord))])
  
  names(dfCoord)[grep("treatment", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("treatment", names(dfCoord))])
  names(dfCoord)[grep("seurat_clusters", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("seurat_clusters", names(dfCoord))])
  names(dfCoord)[grep("Seurat_Clusters_B", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("Seurat_Clusters_B", names(dfCoord))])
  names(dfCoord)[grep("quadrant", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("quadrant", names(dfCoord))])
  names(dfCoord)[grep("longTag", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("longTag", names(dfCoord))])
  names(dfCoord)[grep("isOld", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("isOld", names(dfCoord))])
  names(dfCoord)[grep("integrated_snn_res_0_8", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("integrated_snn_res_0_8", names(dfCoord))])
  
  dfCoord$meta_seurat_clusters <- paste0("Cluster_", dfCoord$meta_seurat_clusters)
  dfCoord$meta_Seurat_Clusters_B <- paste0("Cluster_", dfCoord$meta_Seurat_Clusters_B)
  
  dfdbTable <- dfCoord
  #dfdbTable[["sample_group"]] <- paste0("Cluster_r08_", dfdbTable$louvain_r0_8)
  dfdbTable[["sample_group"]] <- dfdbTable$meta_shortName
  dfdbTable[["sample_group_color"]] <- dfdbTable$color
  
}
# Ifn project

## Rename columns ##
if (project_id == "EGCsHPoly"){
  names(dfCoord) <- gsub("index", "cellID", names(dfCoord))
  
  names(dfCoord)[grep("condition$", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("condition$", names(dfCoord))])
  names(dfCoord)[grep("batch", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("batch", names(dfCoord))])
  names(dfCoord)[grep("condition_label", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("condition_label", names(dfCoord))])
  names(dfCoord)[grep("RNA_snn_res_0_1", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("RNA_snn_res_0_1", names(dfCoord))])
  names(dfCoord)[grep("RNA_snn_res_0_8", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("RNA_snn_res_0_8", names(dfCoord))])
  names(dfCoord)[grep("seurat_clusters", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("seurat_clusters", names(dfCoord))])
  names(dfCoord)[grep("cluster$", names(dfCoord))] <- paste0("meta_", names(dfCoord)[grep("cluster$", names(dfCoord))])
  
  
  dfdbTable <- dfCoord
  #dfdbTable[["sample_group"]] <- paste0("Cluster_r08_", dfdbTable$louvain_r0_8)
  dfdbTable[["sample_group"]] <- dfdbTable$meta_cluster
}

pos <- grep("sample_group_colors", names(dfdbTable))

if (length(pos) == 0){
  dfCol <- data.frame(
    sample_group = unique(sort(dfdbTable$sample_group)),
    sample_group_colors = hue_pal()(length(unique(sort(dfdbTable$sample_group)))
    ))
  
  dfdbTable <- merge(
    dfdbTable,
    dfCol,
    by.x = "sample_group",
    by.y = "sample_group"
  )
  
  
}


dfdbTable$clusterColor <- dfdbTable[["sample_group_colors"]]



###############################################################################
## Add Sample and G2M colors if available                                    ##
if (length(grep("sampleID", names(dfdbTable))) > 0){
  identities <- levels(factor(dfdbTable[,"sampleID"]))
  sample_group_colors <- hue_pal()(length(identities))
  
  dfdbTable[["sampleID_colors"]] <- dfdbTable[["sampleID"]]
  
  for (i in 1:length(identities)){
    dfdbTable$sampleID_colors <- gsub(identities[i], sample_group_colors[i], dfdbTable$sampleID_colors)
  }
}
## G2M colors ##
if (length(grep("Phase", names(dfdbTable))) > 0){
  identities <- levels(factor(dfdbTable[,"Phase"]))
  sample_group_colors <- hue_pal()(length(identities))
  
  dfdbTable[["Phase_colors"]] <- dfdbTable[["Phase"]]
  
  for (i in 1:length(identities)){
    dfdbTable$Phase_colors <- gsub(identities[i], sample_group_colors[i], dfdbTable$Phase_colors)
  }
}
if (length(grep("DF_Classification", names(dfdbTable)) > 0)){
  #######################
  ## Edit singlet and doublet
  dfdbTable$DF_Classification <- gsub(
    "Singlet", 1, dfdbTable$DF_Classification)
  
  dfdbTable$DF_Classification <- gsub(
    "Doublet", 2, dfdbTable$DF_Classification)
  
  dfdbTable$DF_Classification <- as.numeric(dfdbTable$DF_Classification)
  ## done
  ########################
}
##
###############################################################################  

###############################################################################
## Add additonal classifications                                             ##
# dfdbTable[["Sample_Group"]] <- ""
# groupLength <- nchar(dfdbTable$Sample_Group)
# groupLength <- groupLength -1
# dfdbTable[["Sample_Group"]] <- substr(dfdbTable$sampleID, 1, groupLength)
##
###############################################################################
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



columnDBcategoryList <- inferDBcategories(dfData=dfdbTable)


upload.datatable.to.database(
  host = host,
  user = db.user,
  password = db.pwd,
  prim.data.db = primDataDB,
  dbTableName = PCAdbTableName,
  df.data = dfdbTable,
  db.col.parameter.list = columnDBcategoryList,
  new.table = TRUE
)
killDbConnections()
##
#########################

## Done 
###############################################################################

## Create user in db
if (host == "10.27.241.82"){
    serverString <- "shiny-bioinformatics.crick.ac.uk"
    shinyBasePath <- "/camp/stp/babs/www/shiny/boeings/external/"  
    FvDir <- paste0("/camp/stp/babs/www/boeings/bioLOGIC/data/", project_id)
} else if (host == "10.27.241.234"){
    serverString <- "shiny.thecrick.org"
    shinyBasePath <- "/camp/stp/babs/www/shiny/boeings/"
    FvDir <- paste0("/camp/stp/babs/www/boeings/bioLOGIC/data/", project_id)
} else {
    shinyBasePath <- "../"
    FvDir <- paste0("/camp/stp/babs/www/boeings/bioLOGIC/data/", project_id)
    #stop("no server specified")
}  

## Remove if it exists
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
shinyParamPath <-paste0(shinyProjectPath, "/parameters/")
if (!dir.exists(shinyParamPath)){
  dir.create(shinyParamPath)
}
aFN <- paste0(shinyDataPath, "db.txt")

if (file.exists(aFN)){
  df <- read.delim(aFN, header = T, sep="\t", stringsAsFactors = F)
  sPwd <- as.vector(df$id2)
  sUser <- as.vector(df$id)
} else {
    sUser <- substr(paste0(project_id, "_scData"), 1, 15)
    sPwd <-c(2:9,letters,LETTERS)
    sPwd <- paste(sample(sPwd,8),collapse="")
}

query0 <- paste0("SELECT User, Host FROM mysql.user WHERE User = '",sUser,"' AND Host = '10.%';")
res <- doQuery(query = query0, resOut=TRUE)
if (nrow(res) > 0){
    query0a <- paste0("DROP USER '",sUser,"'@'10.%';")   
    doQuery(query = query0a)
    
    
} 

query6 <- paste0(
  "CREATE USER '",
  sUser, 
  "'@'","10.%","' IDENTIFIED BY '",
  sPwd,
  "';"
)



doQuery(query = query6)
query <- c(query, query6)

query0 <- paste0("SELECT User, Host FROM mysql.user WHERE User = '",sUser,"' AND Host = '",serverString,"';")
res <- doQuery(query = query0, resOut=TRUE)
if (nrow(res) > 0){
    query0a <- paste0("DROP USER '",sUser,"'@'",serverString,"';")   
    doQuery(query = query0a)
    
    
} 

## Add user for debugging ##
query6a <- paste0(
  "CREATE USER '",
  sUser, 
  "'@'",serverString,"' IDENTIFIED BY '",
  sPwd,
  "';"
)

doQuery(query = query6a)
query <- c(query, query6a)


query1 <- paste0("DROP TABLE IF EXISTS ", dbtable, ";\n")
query <- query1
doQuery(query = query1)
query2 <- paste0(
    #query,
    "CREATE TABLE IF NOT EXISTS ",
    dbtable,
    " (gene VARCHAR(100) CHARACTER SET latin1 COLLATE latin1_swedish_ci, cellID VARCHAR(100) CHARACTER SET latin1 COLLATE latin1_swedish_ci, lg10Expr DECIMAL(6,3) NULL DEFAULT NULL, row_names INT(10) NOT NULL AUTO_INCREMENT,PRIMARY KEY (row_names)); "
)
doQuery(query = query2)
query <- c(query, query2)
#CREATE TABLE IF NOT EXISTS test_expr (`row_names` bigint(10) NOT NULL, `gene_id` INT(5), `gene` VARCHAR(100) CHARACTER SET latin1 COLLATE latin1_swedish_ci, `lg10Expr` DECIMAL(6,3) NULL DEFAULT NULL, PRIMARY KEY (`row_names`));
query3 <- paste0(
    #query,
    "LOAD DATA LOCAL INFILE '",
    
    "temp.",project_id,".csv",
    "' INTO TABLE ",
    dbtable, 
    " FIELDS TERMINATED BY ',' ENCLOSED BY '\"' LINES TERMINATED BY '\n' IGNORE 1 LINES (gene, cellID, lg10Expr, row_names);"
)
doQuery(query = query3)
query <- c(query, query3)
#LOAD DATA LOCAL INFILE '/camp/stp/babs/working/boeings/Projects//schaefera/tobias.ackels/319_scRNAseq_mm_olfactory_bulb_proj_neuron_analysis_SC19135/workdir/temp.csv' INTO TABLE test_expr FIELDS TERMINATED BY ',' ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 LINES (id, mycol1, mycol2);
query4 <- paste0(
    #query,
    "ALTER TABLE ", dbtable, " ADD UNIQUE(row_names);"
)
doQuery(query = query4)
query <- c(query, query4)
query5 <- paste0(
    #query,
    "CREATE INDEX idx_gene ON ", dbtable, " (gene);"
)
doQuery(query = query5)
query <- c(query, query5)


query7 <- paste0(
    "GRANT SELECT on ", primDataDB,".",PCAdbTableName, " TO ",sUser,"@'",serverString,"';"
)
doQuery(query = query7)
query <- c(query, query7)
query7a <- paste0(
    "GRANT SELECT on ", primDataDB,".",PCAdbTableName, " TO ",sUser,"@'10.%';"
)
doQuery(query = query7a)
query <- c(query, query7a)
query8 <- paste0(
    "GRANT SELECT on ", primDataDB,".",expDbTable, " TO ",sUser,"@'",serverString,"';"
)
doQuery(query = query8)
query <- c(query, query8)
query8a <- paste0(
    "GRANT SELECT on ", primDataDB,".",expDbTable, " TO ",sUser,"@'10.%';"
)
doQuery(query = query8a)
query <- c(query, query8a)
geneTb <- paste0(project_id, "_geneID_tb")
query9 <- paste0(
    "GRANT SELECT on ", primDataDB,".",geneTb, " TO ",sUser,"@'",serverString,"';"
)
doQuery(query = query9)
query <- c(query, query9)
geneTb <- paste0(project_id, "_geneID_tb")
query9a <- paste0(
    "GRANT SELECT on ", primDataDB,".",geneTb, " TO ",sUser,"@'10.%';"
)
doQuery(query = query9a)
query <- c(query, query9a)
# 
# GRANT SELECT on csl_data.p316_rna_seq_table TO csl316data@'shiny-bioinformatics.crick.ac.uk';
## Create relevant directory ##
## Create log-in file ##
url <- host
id <- sUser
id2 <- sPwd
db <- primDataDB
coordTb <- PCAdbTableName
exprTb <- paste0(project_id, "_gene_expr_tb")
geneTb <- paste0(project_id, "_geneID_tb")
#dfPercCellsExpr <- Obio@dataTableList$dfPercCellsExpr
#dfPercCellsExpr <- dfPercCellsExpr[dfPercCellsExpr$gene %in% Obio@dataTableList$referenceList$integrated_top30var, ]
#dfPercCellsExpr <- dfPercCellsExpr[order(dfPercCellsExpr$count_cut_off, decreasing = T),]
default <- as.vector(dfExpr[1,"gene"])
default <- "S100b"

if (!(default %in% dfExpr$gene)){
    stop("Set default gene correctly")
    
}
dfID <- data.frame(
    url,
    id,
    id2,
    db,
    coordTb,
    exprTb,
    geneTb,
    default
)
## Maintain old password if exists
# FN <- paste0(shinyDataPath, "db.txt")
# if (file.exists(FN)){
#     dfCont <- read.delim(
#         FN, 
#         header = T,
#         sep = "\t",
#         stringsAsFactors = F
#     )
#     dfID$id2 <- dfCont$id2
# }
#write.table(dfID, FN, row.names = F, sep="\t")
## Copy rest ##

cpString1 <- paste0("cp -r /camp/stp/babs/www/shiny/boeings/external/scAppTemplate/* ", shinyProjectPath)


#cpString1 <- paste0("cp assets/template_sc_app/ui.r ", shinyProjectPath)
system(cpString1)

FN <- paste0(shinyDataPath, "db.txt")
write.table(dfID, FN, row.names = F, sep="\t")
# cpString2 <- paste0("cp assets/template_sc_app/server.r ", shinyProjectPath)
# system(cpString2)
## Done creating app ##
## Make featureView entry ##
FvDir <- paste0("/camp/stp/babs/www/boeings/bioLOGIC_external/data/", project_id)
if (!dir.exists(FvDir)){
    dir.create(FvDir)
}
FvDir <- paste0(FvDir, "/html/") #html_local_dir
#[1] "/camp/stp/babs/working/boeings/P
#paste0("/camp/stp/babs/www/boeings/bioLOGIC_external/data/", project_id, "/html")
if (!dir.exists(FvDir)){
    dir.create(FvDir)
}
## Feature View Created ##
## Export upload string ##
write.table(
    data.frame(query),
    "query.string.txt",
    row.names = F,
    sep = "\t"
)


##                                                                           ##
###############################################################################







################
##

createHTMLoutput <- function(FvFN = NULL, project_id = NULL){
    sink(FvFN)
    cat('<!DOCTYPE html>');cat("\n");
    cat('<html lang="en">');cat("\n");
    cat('<head>');cat("\n");
    cat('    <meta charset="UTF-8">');cat("\n");
    cat(paste0('<title> ',project_id,' FeatureView</title>    <style>'));cat("\n");
    cat('        body {');cat("\n");
    cat('             margin: 0;            /* Reset default margin */');cat("\n");
    cat('         }');cat("\n");
    cat('        iframe {');cat("\n");
    cat('            display: block;       /* iframes are inline by default */');cat("\n");
    cat('            background: #000;');cat("\n");
    cat('            border: none;         /* Reset default border */');cat("\n");
    cat('            height: 100vh;        /* Viewport-relative units */');cat("\n");
    cat('            width: 100vw;');cat("\n");
    cat('        }');cat("\n");
    cat('    </style>');cat("\n");
    cat('</head>');cat("\n");
    cat('<body>');cat("\n");
    cat(paste0('<iframe src="https://shiny-bioinformatics.crick.ac.uk/shiny/boeings/',project_id,'_app/" title="FeatureView"></iframe>'));cat("\n");
    cat('</body>');cat("\n");
    cat('</html>');cat("\n");
    sink()
    
}
createHTMLoutput(
    FvFN = paste0(FvDir, "FeatureView_", project_id,".html"), 
    project_id = project_id
)

