# Create an Interactive Single-cell App Based on a Seurat single-cell Object

## Run App Locally
```
# Install R-package for app creation
devtools::install_github("decusinlabore/biologicViewerSC")
devtools::install_github("decusInLabore/biologicSeqTools")

# When prompted to update packages, try first the 'None' option. 

```
Once the installation is complete run

```
library(biologicViewerSC)
library(biologicSeqTools)
library(Seurat)

project_id = "test2_PBMC"





###############################################################################
## Step 1: Load Seurat object                                                ##

library(Seurat)
testObj <- pbmc_small
all.genes <- rownames(testObj)
testObj <- ScaleData(testObj, features = all.genes)
testObj <- RunPCA(testObj, npcs = 3)
testObj <- RunUMAP(testObj, dims = 1:3)

testObj <- FindNeighbors(testObj, dims = 1:3)
testObj <- FindClusters(testObj, resolution = 0.5)

testObj@meta.data[["all"]] <- "all"

testObj@meta.data[["sampleName"]] <- "SampleA"
testObj@meta.data[["sampleColor"]] <- "#009900"

testObj@meta.data[["clusterName"]] <- paste0("Cluster_", as.vector(testObj@meta.data$seurat_clusters))

clusterVec <- unique(testObj@meta.data$clusterName)
clusterCols <- scales::hue_pal()(length(clusterVec))
dfCol <- data.frame(clusterVec, clusterCols)

dfCell <- data.frame(cellID = row.names(testObj@meta.data), clusterName = testObj@meta.data$clusterName) 
dfCell <- merge(dfCell, dfCol, by.x = "clusterName", by.y = "clusterVec")

addVec <- dfCell$clusterCols
names(addVec) <- row.names(dfCell)

testObj <- Seurat::AddMetaData(
                object = testObj, 
                metadata = addVec, 
                "clusterColor"
)

testObj@meta.data[["meta_Region"]] <- "Tumor"
testObj@meta.data[1:20,"meta_Region"] <- "Normal"

OsC <- testObj
names(OsC@meta.data) <- gsub("\\.", "_", names(OsC@meta.data))

  

##                                                                           ##
###############################################################################

###############################################################################
## Create Single-cell Application                                            ##

shinyBasePath <- "./"
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
    mode = "SQLite"  # Options: "MySQL" and "SQLite"
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
    mode = "SQLite"  # Options: "MySQL" and "SQLite"
)


###############################################################################
## Create Metadata table                                                     ##

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
   mode = "SQLite"
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
    dataMode = "SQLite",
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







##                                                                           ##
###############################################################################



###############################################################################
## Customise app                                                             ##

params <- biologicSeqTools::scanObjParams(OsC)

## Done                                                                      ##
###############################################################################


## run your app locally

shiny::runApp(paste0(project_id, "_app"))


```



## Option A: Database Version
```
devtools::install_github(
    "FrancisCrickInstitute/hc",
    ref = "main", 
    auth_token = "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"
)
```
Once the installation is complete, run:

```
library(biologicViewerSC)

biologicViewerSC::run_app()
```



## Working on the code of this app on your local computer

### Get the present app from github

Go to github repo for this project https://github.com/FrancisCrickInstitute/hc and select the branch you'd like to work on. 

Get all files to your computer using the path given in the green code button like so:

```
git clone git@github.com:FrancisCrickInstitute/hc.git
```

Create a new branch for your improvement work

```
git checkout -b app_improvements
```

To render and run the app locally from the downloaded code, e.g. after you've made changes, do the following 
```
setwd("path/to/project/folder/hc")

devtools::check()

# Set options here
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode

# Detach all loaded packages and clean your environment
golem::detach_all_attached()
# rm(list=ls(all.names = TRUE))

# Document and reload your package
golem::document_and_reload()

# Run the application
run_app()
```

To check for the name of your current github branch you can do
```
git branch
```

## Deploying this app on a shiny server
To deploy this app on a shiny server, download the content of the repo using the git clone approach described above, make sure cytoscape is installed and copy the PPI folder onto the shiny server. Once all relevant R-packages are installed, the app should run. 

Once you are happy with your improvements, you can push your development branch to the remote github repo by doing (inside the PPI folder)

```
git add -A

git commit -m "Description of the changes you've made"

git push origin [name_of_your_new_branch]
```
