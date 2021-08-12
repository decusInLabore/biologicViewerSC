# Create an Interactive Single-cell App Based on a Seurat single-cell Object

This R-package will create a shiny viewer for your (Seurat-) single-cell RNA-Sequencing dataset. 
As input we need only your analysed Seurat object. The app can be used in two modes: as a 'local' version that stores all required data in files locally, or as a 'remote database' version, that relies on a MySQL/MariaDB database. This database can be on your local machine as well as hosted on a remote server. Using the local version wil be easier in terms of setup, the database version has significant performance advantages particularly when it comes to larger datasets. Both versions of the app can be readily deployed on a shiny server. 

## Option 1 Create App with local data storage. 

Let's start by installing required R-packages:
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

testObj <- pbmc_small %>% 
    ScaleData(features = all.genes) %>% 
    RunPCA(npcs = 3) %>%
    RunUMAP(dims = 1:3)  %>%
    FindNeighbors(dims = 1:3) %>%
    FindClusters(resolution = 0.5)



testObj@meta.data[["meta_Region"]] <- "RandomAcat"
testObj@meta.data[1:20,"meta_Region"] <- "RandomBcat"


##                                                                           ##
###############################################################################
```

If you have a Seurat object containing a single-cell experiment already, just rename it to OsC and create two meta data columns by the name of 'clusterName' and 'sampleName'. The presence of these columns with these names is required. For example you could create them by following the example below. 

<b>
```
###############################################################################
## Prepare your Seurat Object                                                ##

## We are inserting here as an example the small Seurat test object we've created above.


OsC <- testObj


## To use your analysed Seurat object, and rename it to OsC. 


# OsC <- YourAnalysedSeuratObject

names(OsC@meta.data) <- gsub("\\.", "_", names(OsC@meta.data))


##                                                                           ##
###############################################################################
```
</b>

And now we create the app in two steps: In a first step we create a list with project parameters (which you can customise before proceding) and in a second step the single-cell shiny app. You also can later on customise display options in the parameter files in the parameter folder. 

<b>Step 1 Create Parameter List</b>
```

params <- biologicViewerSC::scanObjParams(OsC)

## Review the default app settings and colors:
params

```

<b>Step2 Create app </b>
Depending on the size of your Seurat object, it <b>might take a couple of minutes for the seuratObjectToLocalViewer function to run</b>. Very large single cell objects might have to be processed on a high-performance computing system. A dataset with 5000 cells should take less than a minute to render on your local system. 
<b>

```
project_id <- "testExperiment"

projectPath = "./"


biologicViewerSC::seuratObjectToLocalViewer(
    params = params,
    project_id = project_id,
    projectPath = projectPath,
    OsC = OsC,
    dataMode = "SQLite"
)



```
</b>
Run app locally, e.g. in Rstudio

<b>
```
shiny::runApp(paste0(projectPath, project_id, "_app"))
```
</b>
If you want to deploy the app on a shiny server, simply transfer the project folder, in the example case names test_PBMC_app onto the shiny server. 


## Option B: MySQL/MariaDB Database Version






