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
project_id <- "testApp"

params <- biologicSC::scanObjParams(OsC)

seuratObjectToLocalViewer(
    project_id = project_id,
    projectPath = "./",
    OsC = OsC,
    dataMode = "SQLite"
)

Run app locally

sewd("testApp")
shiny::runApp(paste0("../", project_id))



```



## Option B: MySQL/MariaDB Database Version






