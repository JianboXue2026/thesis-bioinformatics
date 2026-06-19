
# Install and load Seurat
install.packages('Seurat')

library('Seurat')

# Install and load monocle3

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

install.packages("devtools")
install.packages("sf") # sf installed separately because monocle3 direct install may fail to download it
devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)

# Install and load clusterProfiler
BiocManager::install("clusterProfiler")

library(clusterProfiler)


# Install and load CellChat
devtools::install_github("jinworks/CellChat")

install.packages('NMF')
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")

library(CellChat)

# Install and load SingleR
library(BiocManager)
BiocManager::install("SingleR")
BiocManager::install("celldex")
library(SingleR)
