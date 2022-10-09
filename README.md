# Create an Interactive Single-cell App Based on a Seurat single-cell Object

This R-package will create a shiny viewer for your (Seurat-) single-cell RNA-Sequencing dataset. 
As input we need only your Seurat object. The app can be used in two modes: as a 'local' version that stores all required data in files locally, or as a 'remote database' version, that relies on a MySQL/MariaDB database. This database can be on your local machine as well as hosted on a remote server. Using the local version wil be easier in terms of setup, the database version has significant performance advantages particularly when it comes to larger datasets. Both versions of the app can be readily deployed on a shiny server. 

You will get a local Shiny app, that can be used locally or copied onto a shiny server. If you need help with the database setup or the shiny server at the Crick, please do get in touch with the Software development team, the database team or email me (stefan.boeing@crick.ac.uk).

## Example App Neuroblastoma (comparison single-cell / single-nuclei Workflow)
In the section below you will be able to create a single-cell data viewer using your own Seurat single-cell object. As an example you can review a comparison of a single-cell and single-nuclei workflow on a Neuroblastoma sample <a href="https://bioinformatics.crick.ac.uk/shiny/users/boeings/Neuroblastoma_app/" target="_blank">here</a>. The About this dataviewer section will have additional details on the experiment.  

**Features**
* Backend
  * Two-step creation procedure from Seurat object
  * Custom-selection of display options via the params list
  * Convenient resetting of default category colors in the app's parameters/colorParameters.txt file at any time after app creation
  * App can be run with a datafile or remote database
* Front-end
  * Free Category Color Selection
  * All figures can be exported as high-quality PDFs (Download Button)
  * Any snapshot state of the app can be captured and shared as a url link (Bookmark button)
  * User-directed figure creation
  * User-directed data slicing
  * FeatureView plots
  * Histograms
  * Ridgeplots
  * Densityplots
* Performance
  * Small webserver memory footprint
  
  

## Create single-cell Shiny Data Viewer from your Seurat Object

In the example below we use create a small Seurat single-cell object from a single-cell dataset build into Seurat. In order to use your Seurat object, start at the OsC <- YourAnalysedSeuratObject step.


Let's start by installing required R-packages:
```

if (!require("Seurat")){
    install.packages("Seurat")
}

if (!require("remotes")){
    install.packages("remotes")
}

remotes::install_github("decusinlabore/biologicViewerSC")

```


```
###############################################################################
## Step 1: Load Exampe Seurat object                                         ##

library(Seurat)
library(dplyr)
library(biologicViewerSC)

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


And now we create the app in two steps: In a first step we create a list with project parameters (which you can customise before proceding) and in a second step the single-cell shiny app. You also can later on customise display options in the parameter files in the parameter folder. 

<b>Step 1 Create Parameter List</b>
```

params <- biologicViewerSC::scanObjParams(OsC)

## Review the default app settings and colors:
params

```

<b>Step2 Create app </b>
Depending on the size of your Seurat object, it <b>might take a couple of minutes for the seuratObjectToLocalViewer function to run</b>. Very large single cell objects might have to be processed on a high-performance computing system. A dataset with 5000 cells should take less than a minute to render on your local system. 

## Option A: Fully Local Version
```
project_id <- "testExperiment"

projectPath = "./"


biologicViewerSC::seuratObjectToLocalViewer(
    params = params,
    project_id = project_id,
    projectPath = projectPath,
    OsC = OsC,
    dataMode = "SQLite",
    geneDefault = "CD3E" # This will set the default gene in the app
)



```


Run app locally, e.g. in Rstudio


```
shiny::runApp(paste0(projectPath, project_id, "_app"))
```

If you want to deploy the app on a shiny server, simply transfer the project folder, in the example case names test_PBMC_app onto the shiny server. 


## Option 2B: MySQL/MariaDB Remote Database Version

This option is particularly for large datasets. Also it can be useful to minimise the datafootprint, as data only needs to be stored in one place for world-wide availability. 

For this option you need to have access to a MySQL/MariaDB database with permissions to GRANT access to project users. For the database upload we are using the administrator database privileges. For the app to access the database we will create a restricted user that will only be able to read tables for the project it works on. 

### Setup database access credentials

```
dbHostURL <- "mysql.database.IP.address"
dbAdminUser <- "admin.user.name"

FN <- "/camp/stp/babs/working/boeings/Projects/reference_data/documentation/BC.parameters.txt"
dbTable <- read.delim(
    FN, 
    sep = "\t",
    stringsAsFactors = F
)
dbAdminPassword <- as.vector(dbTable[1,1])


project_id <- "testDb"

projectPath = "./"


biologicViewerSC::seuratObjectToViewer(
    params = params,
    project_id = project_id,
    projectPath = projectPath,
    OsC = OsC,
    dataMode = "MySQL",
    host = "10.27.241.82",
    dbname = "test_data",
    db.pwd = dbAdminPassword,
    db.user = "boeings",
    appDomains = c("bioinformatics.crick.ac.uk","10.%"),
    geneDefault = "CD3E" # This will set the default gene in the app
)



```

## FAQs

### How do I get to a Violin Plot?
In the data viewer, select as x-axis a category column, for example, seurat_clusters in the provided example, Meta Region, as y-axis log10 pression and as colorBy Meta Region. 

### How do I get to a Histogram?
In the data viewer, select as x-axis a numeric column, e.g. nFeature_RNA, and as y-axis Histogram and as colorBy, for example, meta Region. In your experiment, splitting by sample might be helpful.

### How do I change colors in the app?
To change category colors, e.g. cluster colors, navigate in the app folder to the parameters folder and open the colorParameters.txt file in Excel and edit the color in the colSel column.

Alternatively you can edit the color HEX codes in the params list prepared in step 1.

### How do I limit the amout of parameters listed as selection options?
If you wish to limit the parameters listed in your viewer, edit the params list in step 1. 

To see all default parameters for the x-axis selection, y-axis selection, splitPlotsBy selection or colorPlotsBy selection, do after performing step 1:

```
names(params[["x_axis"]]) or
names(params[["y_axis"]]) or 
names(params[["splitPlotsBy"]]) or
nanes(params[["colorPlotsBy"]]) or

```

In order to for example remove the TSNE coordinates from being listed as x- and y-axis options, do

```
# Removing options TSNE 1 and TSNE 2 from the x-axis listing:

entriesToRemoveX <- c("TSNE 1", "TSNE 2")
itemsToRetainX <- names(params[["x_axis"]])[!(names(params[["x_axis"]]) %in% entriesToRemoveX)]

params[["x_axis"]] <- params[["x_axis"]][itemsToRetainX]

# Removing options TSNE 1 and TSNE 2 from the y-axis listing:

entriesToRemoveY <- c("TSNE 1", "TSNE 2")
itemsToRetainY <- names(params[["y_axis"]])[!(names(params[["y_axis"]]) %in% entriesToRemoveY)]

params[["y_axis"]] <- params[["y_axis"]][itemsToRetainY]

```

Once you've customised the params list as shown above, re-run the seuratObjectToLocalViewer or seuratObjectToViewer functions to create a local or database Shiny app.  