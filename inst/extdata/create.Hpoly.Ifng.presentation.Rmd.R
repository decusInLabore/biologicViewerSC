###############################################################################
## The following R-script should allow you to create powerpoint presentations## 
##from your single-cell data.                                                ##

## For each slide, you can freely choose a title (Slidetitle). To make a plot, 
## go to the online single-cell Viewer, choose the plot you'd like to 
## incorporate in your presentation and use the bookmark button to generate 
## the required plotlink (figureLink). Each viewer requires the specification 
## of a few parameters which are given in the app-specific appParameterFile and
## the pathToAppColorParameters file

## Set the working directory to where you've copied the four files below or 
## adjust the paths to the correct location of those files, e.g. by doing

## !!! You need to be inside the Crick or VPN connected to the Crick for the 
## following procedure to work !!!

## setwd("/folder/where/you/have/copied/the/files/to")

## Then we install required R-packages:
install.packages("remotes")
remotes::install_github("decusInLabore/biologicSeqTools")

remotes::install_github("decusinlabore/biologicViewerSC")

## Here are build-in example files for the Hpoly dataset:
## here: https://shiny-bioinformatics.crick.ac.uk/shiny/boeings/EGCsHPoly_app
appParameterFile_Hpoly <-   system.file("extdata/examples/Hpoly.db.txt",package = "biologicViewerSC")
pathToAppColorParameters_Hpoly <- system.file("extdata/examples/Hpoly.colorParameters.txt",package = "biologicViewerSC")   

## Here are build-in example files for the Ifng dataset:
## here: https://shiny-bioinformatics.crick.ac.uk/shiny/boeings/GliaIfngrKO_app
appParameterFile_WT_MT_infected <- system.file("extdata/examples/WT_MT_infected.db.txt",package = "biologicViewerSC") 
pathToAppColorParameters_WT_MT_infected <- system.file("extdata/examples/WT_MT_infected.colorParameters.txt",package = "biologicViewerSC") 

## Add as many items to the Slidelist as you'd like to have slides in your 
## presentation

slideList <- list(
    "Slide1" = list(
        "Slidetitle" = "H Poly Expression UMAP Ifngr1",
        "figureLink" = "https://shiny-bioinformatics.crick.ac.uk/shiny/boeings/EGCsHPoly_app/?_inputs_&gene=%22Ifngr1%22&lowColor=%22%23D3D3D3%22&y_axis=%22UMAP_2%22&dotcolor=%22%2300008B%22&showPlotLegend=true&splitByColumn=%22all%22&x_axis=%22UMAP_1%22&colorBy=%22lg10Expr%22&background=%22white%22&dotsize=1",
        "appParameterFile" = appParameterFile_Hpoly,
        "pathToAppColorParameters" = pathToAppColorParameters_Hpoly
    ),
    "Slide2" = list(
        "Slidetitle" = "A UMAP Test Plot",
        "figureLink" = "https://shiny-bioinformatics.crick.ac.uk/shiny/boeings/GliaIfngrKO_app/?_inputs_&gene=%22S100b%22&lowColor=%22%23D3D3D3%22&y_axis=%22UMAP_2%22&dotcolor=%22%2300008B%22&showPlotLegend=true&splitByColumn=%22all%22&x_axis=%22UMAP_1%22&colorBy=%22lg10Expr%22&background=%22white%22&dotsize=1",
        "appParameterFile" = appParameterFile_WT_MT_infected,
        "pathToAppColorParameters" = pathToAppColorParameters_WT_MT_infected 
    ),
    "Slide3" = list(
        "Slidetitle" = "Another UMAP Test Plot",
        "figureLink" = "https://shiny-bioinformatics.crick.ac.uk/shiny/boeings/GliaIfngrKO_app/?_inputs_&MuscularisMacrophages3=%22%23FD7A63%22&Fibroblasts2=%22%23A0B915%22&lowColor=%22%23D3D3D3%22&SmoothMuscleCells=%22%2333D3A6%22&showPlotLegend=false&EndothelialCells2=%22%23E82F55%22&DendriticCells1=%22%2336DCB6%22&EpithelialCells=%22%23C48779%22&MonocyteMacrophages=%22%23EE7A6F%22&DendriticCells2=%22%23A67066%22&y_axis=%22UMAP_2%22&MuscularisMacrophages2=%22%23B706DC%22&dotcolor=%22%2300008B%22&MuscularisMacrophages4=%22%2361BBC6%22&NaturalKillerCells=%22%23621934%22&TCells1=%22%230351B5%22&dotsize=1&colorBy=%22ClusterName%22&CyclingMacrophages=%22%23C6B1F5%22&gene=%22S100b%22&EndothelialCells1=%22%23927E75%22&MesothelialCells2=%22%2335D156%22&EntericGliaCells=%22%23447BE9%22&Fibroblasts1=%22%23A91CC3%22&MesothelialCells1=%22%231069DF%22&MesothelialCells5=%22%23628D55%22&TCells3=%22%23DEA93F%22&OtherLymphoidCells=%22%234412C1%22&InterstitialCellsOfCajal=%22%23F1CB7A%22&Pericytes=%22%23E5AA42%22&MuscularisMacrophages1=%22%237FC9E4%22&splitByColumn=%22all%22&Bcells=%22%23576825%22&MesothelialCells4=%22%23AFC858%22&x_axis=%22UMAP_1%22&background=%22white%22&LymphaticEndothelialCells=%22%235F994B%22&MesothelialCells3=%22%232DA5D6%22&Neutrophils=%22%23203E53%22&TCells2=%22%2362841F%22",
        "appParameterFile" = appParameterFile_WT_MT_infected,
        "pathToAppColorParameters" = pathToAppColorParameters_WT_MT_infected 
    ),
    "Slide4" = list(
        "Slidetitle" = "A UMAP Plot 3",
        "figureLink" = "https://shiny-bioinformatics.crick.ac.uk/shiny/boeings/GliaIfngrKO_app/?_inputs_&MuscularisMacrophages3=%22%23FD7A63%22&Fibroblasts2=%22%23A0B915%22&lowColor=%22%23D3D3D3%22&SmoothMuscleCells=%22%2333D3A6%22&showPlotLegend=false&EndothelialCells2=%22%23E82F55%22&DendriticCells1=%22%2336DCB6%22&EpithelialCells=%22%23C48779%22&MonocyteMacrophages=%22%23EE7A6F%22&DendriticCells2=%22%23A67066%22&y_axis=%22tSNE_2%22&MuscularisMacrophages2=%22%23B706DC%22&dotcolor=%22%2300008B%22&MuscularisMacrophages4=%22%2361BBC6%22&NaturalKillerCells=%22%23621934%22&TCells1=%22%230351B5%22&dotsize=1&colorBy=%22ClusterName%22&CyclingMacrophages=%22%23C6B1F5%22&gene=%22S100b%22&EndothelialCells1=%22%23927E75%22&MesothelialCells2=%22%2335D156%22&EntericGliaCells=%22%23447BE9%22&Fibroblasts1=%22%23A91CC3%22&MesothelialCells1=%22%231069DF%22&MesothelialCells5=%22%23628D55%22&TCells3=%22%23DEA93F%22&OtherLymphoidCells=%22%234412C1%22&InterstitialCellsOfCajal=%22%23F1CB7A%22&Pericytes=%22%23E5AA42%22&MuscularisMacrophages1=%22%237FC9E4%22&splitByColumn=%22genotype%22&Bcells=%22%23576825%22&MesothelialCells4=%22%23AFC858%22&x_axis=%22tSNE_1%22&background=%22white%22&LymphaticEndothelialCells=%22%235F994B%22&MesothelialCells3=%22%232DA5D6%22&Neutrophils=%22%23203E53%22&TCells2=%22%2362841F%22",
        "appParameterFile" = appParameterFile_WT_MT_infected,
        "pathToAppColorParameters" = pathToAppColorParameters_WT_MT_infected 
    ),
    "Slide5" = list(
        "Slidetitle" = "A UMAP Plot 4",
        "figureLink" = "https://shiny-bioinformatics.crick.ac.uk/shiny/boeings/GliaIfngrKO_app/?_inputs_&MuscularisMacrophages3=%22%23FD7A63%22&Fibroblasts2=%22%23A0B915%22&lowColor=%22%23D3D3D3%22&SmoothMuscleCells=%22%2333D3A6%22&showPlotLegend=false&EndothelialCells2=%22%23E82F55%22&DendriticCells1=%22%2336DCB6%22&EpithelialCells=%22%23C48779%22&MonocyteMacrophages=%22%23EE7A6F%22&DendriticCells2=%22%23A67066%22&y_axis=%22Densityplot%22&MuscularisMacrophages2=%22%23B706DC%22&dotcolor=%22%2300008B%22&MuscularisMacrophages4=%22%2361BBC6%22&NaturalKillerCells=%22%23621934%22&TCells1=%22%230351B5%22&dotsize=1&colorBy=%22ClusterName%22&CyclingMacrophages=%22%23C6B1F5%22&gene=%22Malat1%22&EndothelialCells1=%22%23927E75%22&MesothelialCells2=%22%2335D156%22&EntericGliaCells=%22%23447BE9%22&Fibroblasts1=%22%23A91CC3%22&MesothelialCells1=%22%231069DF%22&MesothelialCells5=%22%23628D55%22&TCells3=%22%23DEA93F%22&OtherLymphoidCells=%22%234412C1%22&InterstitialCellsOfCajal=%22%23F1CB7A%22&Pericytes=%22%23E5AA42%22&MuscularisMacrophages1=%22%237FC9E4%22&splitByColumn=%22genotype%22&Bcells=%22%23576825%22&MesothelialCells4=%22%23AFC858%22&x_axis=%22lg10Expr%22&background=%22white%22&LymphaticEndothelialCells=%22%235F994B%22&MesothelialCells3=%22%232DA5D6%22&Neutrophils=%22%23203E53%22&TCells2=%22%2362841F%22",
        "appParameterFile" = appParameterFile_WT_MT_infected,
        "pathToAppColorParameters" = pathToAppColorParameters_WT_MT_infected 
    )
)


## The function below will generate a powerpoint presentation under the name 
## specified in the outputfile in your working directory. 

biologicViewerSC::createPowerpointPresentation(
    slideList = slideList, 
    outputfile="singleCell.powerpoint.presentation.pptx"
)





