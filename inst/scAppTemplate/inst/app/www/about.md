## About this app

### How to explore your single-cell dataset?
You can explore this single-cell dataset in various ways. 

#### Looking up cell-level gene expression (FeaturePlots)
In a first place you may just to review gene expression levels per cell. On many occassions it is helpful to split the display (using the splityBy drop-down menu) according to sampleIDs, cell cycle phase, condition or treatment. 

#### Looking up cell-level gene-category expression (category FeaturePlots)
Similar to gene level expression levels, many experiments will have gene category level gene expression values available. Get a listing for those by typing cat_ into the gene name box. 

#### Looking up category (e.g. cluster) level gene Expression (Violin Plots)

When dealing with lowly expressed genes, it is often helpful to switch the display to Violin plot by selecting a category option, such as cluster name, sample name or cell cycle phase as a-axis selection and, for example, log10 Expression as y-Axis selection. 

#### Reviewing quality control parameters
To review quality parameters of the single-cell experiment, review measures for the number of RNA features as well as percentages of mitochondrial genes in your experiment by selecting the corresponding columns on the x- and y-axis. 

## Creating a data viewer for your single-cell dataset
If you have an R-Seurat object with your single-cell analysis, you can create a dataviewer similar to this one following the instructions in this [github repository](https://github.com/decusInLabore/biologicViewerSC).


