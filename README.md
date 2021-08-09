# Create an Interactive Single-cell App Based on a Seurat single-cell Object

This R-package will create a shiny viewer for your (Seurat-) single-cell RNA-Sequencing dataset. 
As input we need only your analysed Seurat object. The app can be used in two modes: as a 'local' version that stores all required data in files locally, or as a 'remote database' version, that relies on a MySQL/MariaDB database. This database can be on your local machine as well as hosted on a remote server. Using the local version wil be easier in terms of setup, the database version has significant performance advantages particularly when it comes to larger datasets. Both versions of the app can be readily deployed on a shiny server. 

## Option 1 Create App with local data storage. 

Let's get started and install the required R-packges. 
```
devtools::install_github("decusInLabore/biologicSeqTools")
devtools::install_github("decusinlabore/biologicViewerSC")


library(biologicViewerSC)
library(biologicSeqTools)
library(Seurat)

```


```
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
```

If you have a Seurat object containing a single-cell experiment already, just rename it to OsC and create two meta data columns by the name of 'clusterName' and 'sampleName'. The presence of these columns with these names is required. For example you could create them by following the example below. 

```
###############################################################################
## Prepare your Seurat Object                                                ##

# OsC <- YourSeuratObjectHere

OsC@meta.data[["clusterName"]] <- paste0("C", OsC@meta.data$seurat_clusters)
OsC@meta.data[["sampleName"]] <- paste0("C", OsC@meta.data$orig.ident)

##                                                                           ##
###############################################################################
```

And now we create the app in two steps: In a first step we create a list with project parameters (which you can customise before proceding) and in a second step the single-cell shiny app. 

<b>Step 1 Create Parameter List</b>
```

params <- biologicSC::scanObjParams(OsC)

## Review the default app settings and colors:
params

```

<b>Step2 Create app 
```
project_id <- "test_PBMC"

projectPath = "./"


seuratObjectToLocalViewer(
    #params = params,
    project_id = project_id,
    projectPath = projectPath,
    OsC = OsC,
    dataMode = "SQLite"
)



```

Run app locally, e.g. in Rstudio

```
shiny::runApp(paste0(projectPath, project_id, "_app"))
```

If you want to deploy the app on a shiny server, simply transfer the project folder, in the example case names test_PBMC_app onto the shiny server. 


## Option B: MySQL/MariaDB Database Version






