library(Seurat)
library(dplyr)
library(SingleR)
library(celldex)
library(ggplot2)

# Own dataset
PKU.combined.COVID <- readRDS("path/to/your/project/All_cells-20251204.rds")
# Public dataset
seurat_merged <- readRDS("path/to/your/project/all_merged_seurat-20251204.rds")

# Own data
PKU.combined.COVID$celltype <- Idents(PKU.combined.COVID)

merge_map_mine <- c("T cell" = "T Cells",
                    "CD14 Monocyte" = "Myeloid Cells",
                    "NK" = "NK Cells",
                    "Low-quality cell" = "Other Cells",
                    "CD14CD16 Monocyte" = "Myeloid Cells",
                    "B cell" = "B Cells",
                    "Platelet" = "Platelets",
                    "Plasma" = "Plasma Cells",
                    "CD16 Monocyte" = "Myeloid Cells",
                    "DC" = "Myeloid Cells")
# Use recode() to update idents
PKU.combined.COVID$com_cell_type <- recode(PKU.combined.COVID$celltype, !!!merge_map_mine)
# Update idents
Idents(PKU.combined.COVID) <- PKU.combined.COVID$com_cell_type
head(Idents(PKU.combined.COVID))

saveRDS(PKU.combined.COVID,
        file = "path/to/your/project/All_cells-20251204.rds",
        compress = "gzip")

# Public data
# For public dataset
merge_map_pub <- c("T cells" = "T Cells",
                   "Monocytes" = "Myeloid Cells",
                   "Dendritic cells" = "Myeloid Cells",
                   "NK cells" = "NK Cells",
                   "Platelet" = "Platelets",
                   "B cells" = "B Cells",
                   "Neutrophil" = "Other Cells",
                   "NK-T cells" = "NK Cells",
                   "Plasma cells" = "Plasma Cells")
# Use recode() to update idents
seurat_merged$com_cell_type <- recode(seurat_merged$manual_annotation, !!!merge_map_pub)
# Update idents
Idents(seurat_merged) <- seurat_merged$com_cell_type
head(Idents(seurat_merged))
# Save the modified data
saveRDS(seurat_merged, "path/to/your/project/all_merged_seurat-20251204.rds")


# Common cell-type markers with the public database
# Reorder cell types on the x-axis
cell_type_order <- c("T Cells", "B Cells", "Plasma Cells", "NK Cells", "Myeloid Cells", "Other Cells", "Platelets")

# Own dataset dot plot
DotPlot(PKU.combined.COVID, features = c(
  "PPBP", "PF4", #Platelet
  "CD3D", "CD3E", "CD4", "CD8A", "AQP3", "LEF1", "TCF7", "CXCR3", "CCL4", "CCR6", #T-cells
  "CD38", "MZB1", "IGKC", "SDC1", #Plasma cells
  "MS4A1", "CD19", "CD27", "TCL1A", #B-cells
  "NCAM1", "FCGR3A", "NKG7", "GNLY", #NK cells
  "CD14", "CD33", "LYZ" #Myeloid cells
)) + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(limits = cell_type_order)

# Public dataset dot plot
DotPlot(seurat_merged, features = c(
  "PPBP", "PF4", #Platelet
  "CD3D", "CD3E", "CD4", "CD8A", "AQP3", "LEF1", "TCF7", "CXCR3", "CCL4", "CCR6", #T-cells
  "CD38", "MZB1", "IGKC", "SDC1", #Plasma cells
  "MS4A1", "CD19", "CD27", "TCL1A", #B-cells
  "NCAM1", "FCGR3A", "NKG7", "GNLY", #NK cells
  "CD14", "CD33", "LYZ" #Myeloid cells
)) + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(limits = cell_type_order)

# Own dataset UMAP plot
DimPlot(PKU.combined.COVID, reduction = "umap", label = TRUE, cols = c(
  "T Cells" = "#f8766d",
  "Myeloid Cells" = "#cd9600",
  "NK Cells" = "#7cae00",
  "Platelets" = "#00be67",
  "B Cells" = "#00bfc4",
  "Other Cells" = "#00a9ff",
  "Plasma Cells" = "#ff61cc"
))

# Public dataset UMAP plot
DimPlot(seurat_merged, reduction = "umap", label = TRUE, cols = c(
  "T Cells" = "#f8766d",
  "Myeloid Cells" = "#cd9600",
  "NK Cells" = "#7cae00",
  "Platelets" = "#00be67",
  "B Cells" = "#00bfc4",
  "Other Cells" = "#00a9ff",
  "Plasma Cells" = "#ff61cc"
))

# Attempt to apply threshold: only plot cells with expression > 0
# Public data: dimensionality reduction results must be removed before subsetting
temp_seurat_obj <- seurat_merged
temp_seurat_obj@reductions <- list()
filtered_cells_CARD19_pub <- subset(temp_seurat_obj, subset = CARD19 > 0)
filtered_cells_IFNAR2_pub <- subset(temp_seurat_obj, subset = IFNAR2 > 0)
filtered_cells_IL15RA_pub <- subset(temp_seurat_obj, subset = IL15RA > 0)
filtered_cells_LY6E_pub <- subset(temp_seurat_obj, subset = LY6E > 0)
filtered_cells_FHL1_pub <- subset(temp_seurat_obj, subset = FHL1 > 0)

# Own data: direct cell selection
filtered_cells_CARD19_self <- subset(PKU.combined.COVID, subset = CARD19 > 0)
filtered_cells_IFNAR2_self <- subset(PKU.combined.COVID, subset = IFNAR2 > 0)
filtered_cells_IL15RA_self <- subset(PKU.combined.COVID, subset = IL15RA > 0)
filtered_cells_LY6E_self <- subset(PKU.combined.COVID, subset = LY6E > 0)
filtered_cells_FHL1_self <- subset(PKU.combined.COVID, subset = FHL1 > 0)

VlnPlot(filtered_cells_IFNAR2_self, features = "IFNAR2",
        #group.by = "com_cell_type",
        #split.by = "compub",
        group.by = "compub",
        log = TRUE, pt.size = 0,
        ncol = 3, combine = FALSE,
        #idents = c("T Cells", "Myeloid Cells", "NK Cells", "Platelets", "B Cells", "NK-T Cells", "Plasma Cells"),
        #cols = c(
        # "T Cells" = "#f8766d",
        #"Myeloid Cells" = "#cd9600",
        #"#NK Cells" = "#7cae00",
        #"Platelets" = "#00be67",
        #"B Cells" = "#00bfc4",
        #"NK-T Cells" = "#c77cff",
        #"Plasma Cells" = "#ff61cc"
        #)
)
