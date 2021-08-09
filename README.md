# Create an Interactive Single-cell App Based on a Seurat single-cell Object

This R-package will create a shiny viewer for your (Seurat-) single-cell RNA-Sequencing dataset. 
As input we need only your analysed Seurat object. The app can be used in two modes: as a 'local' version that stores all required data in files locally, or as a 'remote database' version, that relies on a MySQL/MariaDB database. This database can be on your local machine as well as hosted on a remote server. Using the local version wil be easier in terms of setup, the database version has significant performance advantages particularly when it comes to larger datasets. Both versions of the app can be readily deployed on a shiny server. 

## Option 1 Create App with local data storage. 

Let's get started and install the required R-packges. 
```
install.packages("devtools")
devtools::install_github("decusInLabore/biologicSeqTools")
devtools::install_github("decusinlabore/biologicViewerSC")

```


```
###############################################################################
## Step 1: Load Exampe Seurat object                                         ##

library(Seurat)
library(dplyr)
library(biologicViewerSC)
library(biologicSeqTools)


all.genes <- rownames(pbmc_small)

testObj <- ScaleData(pbmc_small, testObj, features = all.genes)

testObj <- pbmc_small %>% 
    ScaleData(features = all.genes) %>% 
    RunPCA(npcs = 3) %>%
    RunUMAP(dims = 1:3)  %>%
    FindNeighbors(dims = 1:3) %>%
    FindClusters(resolution = 0.5)



testObj@meta.data[["meta_Region"]] <- "Tumor"
testObj@meta.data[1:20,"meta_Region"] <- "Normal"

OsC <- testObj


  

##                                                                           ##
###############################################################################
```

If you have a Seurat object containing a single-cell experiment already, just rename it to OsC and create two meta data columns by the name of 'clusterName' and 'sampleName'. The presence of these columns with these names is required. For example you could create them by following the example below. 

```
###############################################################################
## Prepare your Seurat Object                                                ##

# OsC <- YourSeuratObjectHere

OsC@meta.data[["clusterName"]] <- paste0("C", OsC@meta.data$seurat_clusters)
OsC@meta.data[["sampleName"]] <- paste0("S", OsC@meta.data$orig.ident)
names(OsC@meta.data) <- gsub("\\.", "_", names(OsC@meta.data))

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
Depending on the size of your Seurat object, it might take a couple of minutes for the seuratObjectToLocalViewer function to run. Very large single cell objects might have to be processed on a high-performance computing system. 

```
project_id <- "test_zf"

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






