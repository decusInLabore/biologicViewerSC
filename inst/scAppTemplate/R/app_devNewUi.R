#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import ggplot2
#' @import DT
#' @import colourpicker
#' @import DBI
#' @import RMySQL
#' @import RSQLite
#' @noRd


###############################################################################
## Load parameter file if available                                          ##
FNparameters <- "parameters/menuParameters.txt"

if (file.exists(FNparameters)){
    dfParam <- read.delim(
        FNparameters, 
        header = T, 
        sep = "\t",
        stringsAsFactors = F
    )
    
    parameterFileLoaded <- TRUE
    ## Check file integrity ##
    if (!(sum(names(dfParam) %in% c("menuName", "displayName", "colSel", "displayOrder")))){
        rm(dfParam)
        parameterFileLoaded <- FALSE
    }
    
} else {
    parameterFileLoaded <- FALSE
}

## Done                                                                      ##
###############################################################################

###############################################################################
## Load category color file if available                                     ##
FNcolParameters <- "parameters/colorParameters.txt"

if (file.exists(FNparameters)){
  dfCol <- read.delim(
    FNcolParameters, 
    header = T, 
    sep = "\t",
    stringsAsFactors = F
  )
  
  colorFileLoaded <- TRUE
  ## Check file integrity ##
  if (!(sum(names(dfCol) %in% c("menuName", "displayName", "colSel", "displayOrder")))){
    rm(dfCol)
    colorFileLoaded <- FALSE
  }
  
} else {
  colorFileLoaded <- FALSE
}


##                                                                           ##
###############################################################################

###############################################################################
## Data Access Module                                                        ##

FNkey <- "data/connect/db.txt"
FNrda <- "data/dfkey.rda"

if (file.exists(FNkey)){
    dfkey <- read.delim(FNkey, stringsAsFactors = F, sep="\t")
} else if (file.exists(FNrda)){
    load(FNrda)
} else {
    data("dfkey")
} 


geneDefault = as.character(dfkey$default)
host <- as.character(dfkey$url)
user <- as.character(dfkey$id)
DBpd <- as.character(dfkey$id2)
dbname <- as.character(dfkey$db)
coordinateTbName <- as.character(dfkey$coordTb)
exprTbName <- as.character(dfkey$exprTb)
geneID_TbName <- as.character(dfkey$geneTb)

pos <- grep("dataMode", names(dfkey))
if (length(pos) == 1){
    if (dfkey$dataMode == "SQLite"){
        dataMode <- dfkey$dataMode
    } else {
        dataMode <- "MySQL"
    }
} else {
    dataMode <- "MySQL"
}

## Done Data Access Module                                                   ##
###############################################################################

###############################################################################
##  Retrieve Meta data table                                                 ##      

oldw <- getOption("warn")
options(warn = -1)

if (dataMode == "SQLite"){
  
  dbDB <- DBI::dbConnect(
    drv = RSQLite::SQLite(),
    dbname=dbname
  )
  
} else {
  
  dbDB <- DBI::dbConnect(
    drv = RMySQL::MySQL(),
    user = user, 
    password = DBpd, 
    host = host, 
    dbname=dbname
    
  )
  
}


query <- paste0("SELECT DISTINCT * FROM ", coordinateTbName)
dfCoordSel <- DBI::dbGetQuery(dbDB, query)
DBI::dbDisconnect(dbDB)

## Add column for all values
dfCoordSel[["all"]] <- "all"
##                                                                           ##
###############################################################################

###############################################################################
## Create order in which samples are displayed                               ##

if (colorFileLoaded){
  dfSO <- dfCol[grep("^SAMPLENAME$", toupper(dfCol$menuName)), ]  
  if (nrow(dfSO) > 0){
    dfSO <- unique(dfSO[,c("colOption", "displayOrder")])
    dfSO <- dfSO[order(dfSO$displayOrder, decreasing = F),]
    conditionVec <- as.vector(dfSO$colOption)
  } else {
    conditionVec <- unique(sort(dfCoordSel$sampleName))  
  }
} else {
  pos <- grep("sampleOrder", names(dfCoordSel))
  
  if (length(pos) > 0){
    dfCoordSel <- unique(dfCoordSel[,c("sampleName", "sampleOrder")])
    dfCoordSel <- dfCoordSel[order(dfCoordSel$sampleOrder, decreasing = F),]
    conditionVec <- as.vector(dfCoordSel$sampleName)
  } else {
    conditionVec <- unique(sort(dfCoordSel$sampleName))  
  }
}


Nsamples <- length(conditionVec)

##                                                                           ##
###############################################################################

###############################################################################
## Select Display options                                                    ##       

## Create x and y axis selections if no parameterfile is loaded ##
if (!parameterFileLoaded){
  allOptions <- names(dfCoordSel)
  
  ## Remove common undesirable colums ##
  rmNameVec <-c(
    "^DC",
    "uniquecellID",
    "hmIdent",
    "old_ident",
    "cellID", 
    "sample_group",
    "DF_pANN",
    "clusterColor",
    "sampleColor",
    "clustIdent",
    "G2M_Score",
    #"DM_Pseudotime",
    "^Sub_clusters_ExNeurons$",
    "sample_group_colors",
    "row_names",
    "sampleID"
  )
  
  rmVec <- as.vector(NULL, mode = "numeric")
  for (i in 1:length(rmNameVec)){
    rmVec <- c(
      rmVec,
      grep(rmNameVec[i], names(dfCoordSel))
    )
  }
  
  XYsel <- allOptions
  if (length(rmVec) > 0){
    XYsel <- XYsel[-rmVec]
  }
  
  ## Reorder
  XYsel <- c(
    XYsel[grep("UMAP_", XYsel)],
    XYsel[grep("tSNE_", XYsel)],
    XYsel[grep("sampleName", XYsel)],
    XYsel[grep("clusterName", XYsel)],
    XYsel[grep("ClusterTame", XYsel)],
    XYsel[grep("ClusterTest", XYsel)],
    XYsel[grep("PC_", XYsel)],
    XYsel[grep("DM_Pseudotime", XYsel)],
    XYsel[grep("meta_", XYsel)],
    #XYsel[grep("DF_Classification", XYsel)],
    XYsel[grep("nCount", XYsel)],
    XYsel[grep("nFeatures", XYsel)],
    XYsel[grep("nCount", XYsel)],
    XYsel[grep("percent", XYsel)],
    XYsel[grep("nCount", XYsel)],
    XYsel[grep("nCount", XYsel)]
  )
  
  Xsel <- XYsel
  Ysel <- XYsel
  xDisplayName <- "Choose a X-axis"
  yDisplayName <- "Choose a Y-axis"
} else {
  Xsel <- as.vector(dfParam[dfParam$displayName == "x_axis","colSel"])
  Ysel <- as.vector(dfParam[dfParam$displayName == "y_axis","colSel"])
  
  xDisplayName <- gsub("_", " ", unique(dfParam[dfParam$displayName == "x_axis", "menuName"]))
  yDisplayName <- gsub("_", " ", unique(dfParam[dfParam$displayName == "y_axis", "menuName"]))
}

## check if all column names are valid ##
check <- c("lg10Expr", names(dfCoordSel))
Xsel <- Xsel[Xsel %in% check]
Ysel <- Ysel[Ysel %in% check]

defaultX <- "UMAP_1"
if (length(defaultX %in% Xsel) != 1){
  defaultX <- XYsel[2]
}

defaultY <- "UMAP_2"
if (length(defaultX %in% Ysel) != 1){
  defaultX <- XYsel[3]
}


## Add to dropdownlist 
dropDownList <- list()

dropDownList[["x_axis"]] <- list(
  "displayName" = xDisplayName,
  "selOptions" = Xsel,
  "selDisplayOptions" = gsub("_", " ", Xsel),
  "default" = defaultX
)

## Add to dropdownlist 
dropDownList[["y_axis"]] <- list(
  "displayName" = yDisplayName,
  "selOptions" = Ysel,
  "selDisplayOptions" = gsub("_", " ", Ysel),
  "default" = defaultY
)


##                                                                           ##
###############################################################################


###############################################################################
## Set color options                                                         ##

if(!parameterFileLoaded){
  allColorOptions <- c(
    "log10 Expression" = "lg10Expr",
    names(dfCoordSel)
  )
  
  if (!parameterFileLoaded){
    #############################################
    ## Make Color Selections 
    ## Get color selection ##
    allColorOptions <- c(
      #"Log10 Expresson" = "lg10Expr",
      #"DM Pseudotime"  = "DM_Pseudotime",
      "Sample" = "sampleName",
      "Cluster" = "clusterName",
      "Subcluster" = "subClusterName",
      # "WT vs. IDH" = "WT_IDH",
      "Gender" = "Gender",
      #  "Norm vs Hyp" = "Norm_Hyp",
      #  "Con Prad AZ" = "Con_Prad_AZ",
      "Cells From Tumor" = "CellFromTumor",
      "Patient" = "Patient",
      "Region" = "Region",
      "Article Cell Type" = "Article_Cell_Type",
      "Doublet Classification" = "DF_Classification" ,
      "nCount_RNA" = "nCount_RNA",
      "nFeature_RNA" = "nFeature_RNA",
      "percent_mt" = "percent_mt",
      "S Phase Score" = "S_Score",
      "G2M Score" = "G2M_Score",
      "Cell Cycle Phase" = "Phase",
      "Uniform" = "all"
    )
    
    colAddvec <- c(
      XYsel[grep("meta_", XYsel)],
      XYsel[grep("ClusterTestRes", XYsel)]
    )
    
    names(colAddvec) <- colAddvec
    
    allColorOptions <- c(
      allColorOptions, 
      colAddvec
    )
    
    
    allColorOptions <- allColorOptions[allColorOptions %in% names(dfCoordSel)]
    
    c(
      "Log10 Expression" = "lg10Expr",
      allColorOptions
    )
    
  } 
  
} else {
  ## If paramsfile is loaded 
  allColorOptions <- unique(dfParam[dfParam$displayName == "Color Plots By", "colSel"])
  names(allColorOptions) <- gsub("_", " ", unique(dfParam[dfParam$displayName == "Color Plots By", "colSel"]))
}


## Organise order ##
headVec <- unique(
  c(
    grep("LG10EXPR", toupper(allColorOptions)),
    grep("CLUSTERNAME", toupper(allColorOptions)),
    grep("CLUSTER", toupper(allColorOptions)),
    grep("SAMPLENAME", toupper(allColorOptions)),
    grep("META_", toupper(allColorOptions)),
    grep("CLUSTERNAME", toupper(allColorOptions)),
    grep("subClusterName", toupper(allColorOptions))
  )
)

if (length(headVec) > 0){
  headOptions <- allColorOptions[headVec]
  restVec <- allColorOptions[-headVec]
  allColorOptions <- c(
    headOptions,
    restVec
  )
}

if (length(allColorOptions[allColorOptions == "all"]) > 0){
  names(allColorOptions[allColorOptions == "all"]) <- "Unicolor"
}

## check if all column names are valid ##
check <- c("lg10Expr", names(dfCoordSel))
allColorOptions <- allColorOptions[allColorOptions %in% check]


defaultCol <- "lg10Expr"
if (length(defaultCol %in% allColorOptions) != 1){
  defaultCol <- allColorOptions[1]
}

## Add to dropdownlist 
dropDownList[["colorBy"]] <- list(
  "displayName" = xDisplayName,
  "selOptions" = allColorOptions,
  "selDisplayOptions" = gsub("_", " ", names(allColorOptions)),
  "default" = defaultCol
)
## Done with color options                                                   ##
###############################################################################


###############################################################################
## Set split options                                                         ##
if (!parameterFileLoaded){
  splitOptions <- names(dfCoordSel)
  
  rmVec <- c(
    grep("orig_", splitOptions),
    grep("sampleID", splitOptions),
    grep("old_ident", splitOptions),
    grep("hmIdent", splitOptions),
    grep("color", tolower(splitOptions))
    
  )
  
  if (length(rmVec) > 0){
    splitOptions <- splitOptions[-rmVec]
  }
  
} else {
  ## If paramsfile is loaded 
  splitOptions <- unique(dfParam[dfParam$menuName == "splitPlotsBy", "colSel"])
  names( splitOptions) <- gsub("_", " ", unique(dfParam[dfParam$menuName == "splitPlotsBy", "colSel"]))
}


## Remove all split options with more than 20 options ##
Nopt <- apply(dfCoordSel[,splitOptions], 2, function(x) length(unique(x)))
Nopt <- sort(Nopt[Nopt < 42], decreasing = F)


splitOptions <- as.vector(names(Nopt))


Nsamples <- length(unique(dfCoordSel$sampleName))

if (Nsamples > 3 | nrow(dfCoordSel) < 5000){
  headVec <- c(
    grep("all", splitOptions),
    grep("sampleName", splitOptions),
    grep("meta_", tolower(splitOptions)),
    grep("clusterName", splitOptions),
    grep("subClusterName", splitOptions)
  )  
} else {
  headVec <- c(
    grep("sampleName", splitOptions),
    grep("meta_", tolower(splitOptions)),
    grep("all", splitOptions),
    grep("clusterName", splitOptions),
    grep("subClusterName", splitOptions)
  )  
}



if (length(headVec) > 0){
  headOptions <- splitOptions[headVec]
  restVec <- splitOptions[-headVec]
  splitOptions <- c(
    headOptions,
    restVec
  )
}

names(splitOptions) <- splitOptions
names(splitOptions) <- gsub("meta_", "", names(splitOptions) )
names(splitOptions) <- gsub("META_", "", names(splitOptions) )
names(splitOptions) <- gsub("meta_", "", names(splitOptions) )
names(splitOptions) <- gsub("sampleName", "Sample", names(splitOptions) )
names(splitOptions) <- gsub("clusterName", "Cluster", names(splitOptions) )
names(splitOptions) <- gsub("all", "None", names(splitOptions) )

numOptions <- names(dfCoordSel)[!(names(dfCoordSel)) %in% splitOptions]
numOptions <- c(
  "lg10Expr",
  numOptions
)
  
  if (length(splitOptions[splitOptions == "all"]) > 0){
    names(splitOptions[splitOptions == "all"]) <- "None"
  }
  
  ## check if all column names are valid ##
  check <- c(names(dfCoordSel))
  splitOptions <- splitOptions[splitOptions %in% check]
  
  
  defaultS <- splitOptions[1]
  
  sDisplayName <- gsub("_", " ", unique(dfParam[dfParam$menuName == "colorPlotsBy", "displayName"]))
  
  ## Add to dropdownlist 
  dropDownList[["colorBy"]] <- list(
    "displayName" = sDisplayName,
    "selOptions" = splitOptions,
    "selDisplayOptions" = gsub("_", " ", names(splitOptions)),
    "default" = defaultS
  )
  ## Done setting split options                                                ##
  ###############################################################################
  
  
  
  numOptions <- names(dfCoordSel)[!(names(dfCoordSel)) %in% splitOptions]
  numOptions <- c(
    "lg10Expr",
    numOptions
  )
  
  ##                                                                           ##
  ###############################################################################
  

###############################################################################
## Main App - server                                                         ##       


app_ui <- function(request) {
  
  # theme <- bslib::bs_theme(
  #   bg = "#d3d3d3", 
  #   fg = "#EEE8D5", 
  #   primary = "#2AA198" #,
  #   # bslib also makes it easy to import CSS fonts
  #   #base_font = bslib::font_google("Pacifico")
  # )
  
  #theme = bs_theme(version = 4, bootswatch = "cosmo")
  
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
                                                                         ##       
    ## Start new
    navbarPage(
        "biologic SC",   
        
        tabPanel("FeatureView", 
            sidebarPanel(
                mod_FeatureViewSidebar_ui("FeatureViewSidebar_ui_1")
                
            ),
            mainPanel(
                fluidRow(
                    textOutput("dev_text")
                ),
                fluidRow(
                    column(8,
                        uiOutput("multi_plot_ui")
                    )
                )
            )
        ),
        # tabPanel("Downloads", 
        #         sidebarPanel("sidebar panel 2"),
        #         mainPanel("main panel",
        #             downloadButton('plotDL2')
        #         )
        # ),
        # tabPanel("CategoryView", "three"),
        navbarMenu("Settings", 
                   tabPanel("Set Colors", "four-a"),
                   tabPanel("About", "four-b")
        )
    ) ## End Navbar Parge
    ## End new
    
    
    
    ## old below ##
    
    )
    
    ##                                                                           ##
    ###############################################################################
    
  
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
 
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'biologicSC'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

