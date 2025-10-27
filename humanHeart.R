
# 1.    I download the data ( Heart Human Cell Atlas - anger Institute) from here:
#       "      https://explore.data.humancellatlas.org/projects/ad98d3cd-26fb-4ee3-99c9-8a2ab085e737/project-matrices      "
#       Downloaded data are in .loom format
#       They are on my PC
#       I put my wd directory to look at the directory I have the downloaded data

# I have moved the R script and files to the directory right down.
getwd()
setwd("/Users/nikolasgiannak/heartHumanCellAtlas/")
getwd()

# 2.    Install packages/ Load libraries/ Build environment

library(reticulate)
# Create a fresh environment with all compatible packages
reticulate::conda_create("loom_r_env", python_version = "3.9")

# Install all needed packages with compatible versions
reticulate::conda_install(
  "loom_r_env", 
  packages = c("numpy<2.0", "scipy=1.11.0", "pandas", "anndata", "loompy"),
  pip = FALSE
)

# Use this new environment
use_condaenv("loom_r_env", required = TRUE)

# Import loompy
loompy <- reticulate::import("loompy")

# Restart R


library(reticulate)

# Use the new environment BEFORE loading any other packages or importing Python modules
use_condaenv("loom_r_env", required = TRUE)

# Verify it's using the right Python
py_config()

# Now import loompy
loompy <- reticulate::import("loompy")
class(loompy)

# Load other R packages AFTER setting up Python
library(LoomExperiment)
library(SingleCellExperiment)
library(sceasy)
library(loomR)


# Convert between files
# From loom --> h5ad --> rds



#  I will use sceasy package to convert the .loom file to .sce (Single Cell Experiment)
sceasy::convertFormat('human-heart-10XV2.loom', from="loom", to="sce",
                      outFile='human-heart-10XV2.rds')


sceasy::convertFormat('human-heart-10XV2.loom', from="loom", to="anndata",
                      outFile='humanHeart10XV2.h5ad')



sceasy::convertFormat("humanHeart10XV2.h5ad", from="anndata", to="seurat",
                      outFile='humanHeart10XV2.rds')

library(Seurat)

humanHeart <- readRDS(file = "humanHeart10XV2.rds")

#Recreate Seurat object from scratch
# Extract the count matrix
counts <- GetAssayData(humanHeart, layer = "counts")

# Create a fresh Seurat object
humanHeart <- CreateSeuratObject(counts = counts, 
                                 project = "humanHeart",
                                 min.cells = 3, 
                                 min.features = 200)

# Data normalize Normalization
humanHeart <- NormalizeData(humanHeart, 
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000)


humanHeart <- NormalizeData(humanHeart, normalization.method = "LogNormalize", scale.factor = 10000)

dim(humanHeart)


humanHeart <- FindVariableFeatures(humanHeart, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(humanHeart), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(humanHeart)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Data scaling
all.genes <- rownames(humanHeart)
humanHeart <- ScaleData(humanHeart, features = all.genes)

# Get all genes
all_genes <- rownames(humanHeart)

# Split into chunks
chunk_size <- 5000
chunks <- split(all_genes, ceiling(seq_along(all_genes) / chunk_size))

# Scale in chunks
for (i in seq_along(chunks)) {
  message(paste("Processing chunk", i, "of", length(chunks)))
  humanHeart <- ScaleData(humanHeart, 
                          features = chunks[[i]], 
                          verbose = FALSE)
}


humanHeart <- RunPCA(humanHeart, features = VariableFeatures(object = humanHeart))


DimPlot(humanHeart, reduction = "pca") + NoLegend()

ElbowPlot(humanHeart)

humanHeart <- FindNeighbors(humanHeart, dims = 1:20)
humanHeart <- FindClusters(humanHeart, resolution = 0.5)

DimPlot(humanHeart, reduction = "pca")

humanHeart <- RunUMAP(humanHeart, dims = 1:20)
humanHeart2 <- RunUMAP(humanHeart, dims = 1:10)

DimPlot(humanHeart, reduction = "umap")
DimPlot(humanHeart2, reduction = "umap")



# Now I will use a publicly available atlas to annotate the cells in the humnaHeart seurat object
# install.packages("BiocManager")
#BiocManager::install("celldex")
library(celldex)
surveyReferences()
#ref <- fetchReference("hpca", "2024-02-26")

# This package includes also datasets that can be fetched (downloaded). For studies that generate multiple datasets, 
# the dataset of interest must be explicitly requested via the path= argument:
BiocManager::install("scRNAseq")
library(scRNAseq)
all.ds <- surveyDatasets()
all.ds
# By default, array data is loaded as a file-backed DelayedArray from the HDF5Array package.
# Setting realize.assays=TRUE and/or realize.reduced.dims=TRUE will coerce these to more conventional
# in-memory representations like ordinary arrays or dgCMatrix objects.
BiocManager::install("HDF5Array")
library(HDF5Array)
# THIS IS MY REFERENCE ATLAS
ref <- celldex::HumanPrimaryCellAtlasData()

# Step 1 - extratc the normalized matrix -   till seurat v4 arg was GetAssayData(seuratObj, slot = "data), but for Seurat v5:
 matrix <- LayerData(humanHeart, layer = "data")

# Step 2 - run singleR
 #BiocManager::install("SingleR")
 library(SingleR)
 res.singleR.main <- SingleR( test = matrix, ref = ref, labels = ref$label.main)

 plotScoreHeatmap(res.singleR.main)
 plotDeltaDistribution(res.singleR.main, ncol =4, size = 0.1)
 install.packages("hexbin")
 library(hexbin)
 BiocManager::install("scrapper")
 library(scrapper)
 commoncells <- intersect(colnames(ref), rownames(res.singleR.main))
 
 # Subset the objects to keep only the common cells
 ref_filtered <- ref[, commoncells]
 groups_filtered <- res.singleR.main[commoncells, ]
 
 # Now, use the filtered objects in your plotMarkerHeatmap call
 plotMarkerHeatmap(groups_filtered, ref_filtered, label = "Hepatocytes")
 
 
 
 
 
 
 plotMarkerHeatmap(res.singleR.main, ref,label = "Hepatocytes")
 # Step 3 - add the labels to the seurat object 
 humanHeart$SingleR.labels = res.singleR.main$labels
 
 DimPlot(humanHeart, group.by = "SingleR.labels", label = TRUE, repel = TRUE)
 
 
 DimPlot(humanHeart, reduction = "umap", group.by = "SingleR.labels", label = TRUE, repel = TRUE) 
 #
 
 # check fine labels instead of main
 res.singleR.fine <- SingleR( test = matrix, ref = ref, labels = ref$label.fine)
 
 humanHeart$SingleR.labelsfine = res.singleR.fine$labels
 
 DimPlot(humanHeart, group.by = "SingleR.labels", label = TRUE, repel = TRUE)
 
 
 DimPlot(humanHeart, reduction = "umap", group.by = "SingleR.labels", label = TRUE, repel = TRUE) 
 
 
 
 
# Check if there's a mismatch
length(rownames(humanHeart@assays$RNA@data))
nrow(humanHeart@assays$RNA@meta.features)
