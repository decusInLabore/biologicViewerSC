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
#' @noRd


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

## Done Data Access Module                                                   ##
###############################################################################

###############################################################################
##  Retrieve dfCoordSel                                                      ##       
oldw <- getOption("warn")
options(warn = -1)
dbDB <- DBI::dbConnect(RMySQL::MySQL(), user = user, password = DBpd, host = host, dbname=dbname)
query <- paste0("SELECT DISTINCT * FROM ", coordinateTbName)
dfCoordSel <- DBI::dbGetQuery(dbDB, query)
DBI::dbDisconnect(dbDB)

dfCoordSel[["all"]] <- "all"
##                                                                           ##
###############################################################################


###############################################################################
## Select Display options                                                    ##       

pos <- grep("sampleOrder", names(dfCoordSel))

if (length(pos) > 0){
    dfOrder <- unique(dfCoordSel[,c("sampleName", "sampleOrder")])
    dfOrder <- dfOrder[order(dfOrder$sampleOrder, decreasing = F),]
    conditionVec <- as.vector(dfOrder$sampleName)
} else {
    conditionVec <- unique(sort(dfCoordSel$sampleName))  
}

Nsamples <- length(conditionVec)


allOptions <- names(dfCoordSel)

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
allColorOptions <- 
  c(
    "Log10 Expression" = "lg10Expr",
    allColorOptions
  )


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


## Remove all split options with more than 20 options ##
Nopt <- apply(dfCoordSel[,splitOptions], 2, function(x) length(unique(x)))
Nopt <- sort(Nopt[Nopt < 25], decreasing = F)


splitOptions <- as.vector(names(Nopt))

headVec <- c(
    grep("meta_", tolower(splitOptions)),
    grep("sampleName", splitOptions),
    grep("clusterName", splitOptions),
    grep("subClusterName", splitOptions)
)

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

##                                                                           ##
###############################################################################

###############################################################################
##  Get full gene list                                                       ##       

oldw <- getOption("wafrn")
options(warn = -1)

dbDB <- DBI::dbConnect(MySQL(), user = user, password = DBpd, host = host, dbname=dbname)
query <- paste0("SELECT DISTINCT gene FROM ", geneID_TbName)
allGenes <- as.vector(DBI::dbGetQuery(dbDB, query)[,"gene"])
dbDisconnect(dbDB)

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

