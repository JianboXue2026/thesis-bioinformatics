library(Seurat)
library(Matrix)
library(stringr)

install.packages('devtools')
devtools::install_github('immunogenomics/presto')

# Set data directory
data_dir <- "path/to/your/public/data/Data_Renamed"   # Modify to your actual path

# Search for all matrix files; sample name = prefix before "-matrix.mtx"
matrix_files <- list.files(data_dir, pattern = "-matrix.mtx$", full.names = TRUE)
sample_names <- str_replace(basename(matrix_files), "-matrix.mtx", "")

seurat_list <- list()

for (i in seq_along(sample_names)) {
  sample <- sample_names[i]
  message("Processing: ", sample)

  matrix_file   <- file.path(data_dir, paste0(sample, "-matrix.mtx"))
  features_file <- file.path(data_dir, paste0(sample, "-features.tsv"))
  barcodes_file <- file.path(data_dir, paste0(sample, "-barcodes.tsv"))

  # Read data
  mat <- readMM(matrix_file)
  genes <- read.table(features_file, header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.table(barcodes_file, header = FALSE, stringsAsFactors = FALSE)

  rownames(mat) <- genes$V1
  colnames(mat) <- barcodes$V1

  # Prevent duplicate cell barcodes: add sample prefix
  colnames(mat) <- paste(sample, colnames(mat), sep = "_")

  # Create Seurat object and add sample info
  seu <- CreateSeuratObject(counts = mat, project = sample, min.cells = 0, min.features = 0)
  seu$sample <- sample

  seurat_list[[sample]] <- seu
}

# Merge all Seurat objects
# Keep all features
seurat_merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))

# Inspect
table(seurat_merged$sample)
head(seurat_merged@meta.data)

# Save (optional)
# saveRDS(seurat_merged, "all_merged_seurat.rds")
seurat_merged <- readRDS("path/to/your/public/data/all_merged_seurat.rds")
saveRDS(seurat_merged, "path/to/your/public/data/all_merged_seurat-20251004.rds")

seurat_merged$group <- sub("^([A-Za-z]+)\\d+-(\\d+)$", "\\1-\\2", seurat_merged$sample)
table(seurat_merged$group)
head(seurat_merged$group)
# Normalization
seurat_merged <- NormalizeData(seurat_merged, normalization.method = "LogNormalize", scale.factor = 10000)
# Variable features
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 2000)
# Data scaling
seurat_merged <- ScaleData(seurat_merged, features = VariableFeatures(seurat_merged))

seurat_merged <- RunPCA(seurat_merged, features = VariableFeatures(seurat_merged))
ElbowPlot(seurat_merged)  # Inspect PC selection
plot1 <- DimPlot(seurat_merged, reduction = "pca", group.by="orig.ident")
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:10)   # Typically use first 20 PCs
seurat_merged <- FindClusters(seurat_merged, resolution = 0.3) # resolution is tunable

seurat_merged <- RunUMAP(seurat_merged, dims = 1:10)
DimPlot(seurat_merged, reduction = "umap", group.by = "group")  # View custom grouping
DimPlot(seurat_merged, reduction = "umap", label = TRUE)

seurat_merged <- JoinLayers(seurat_merged)
markers <- FindAllMarkers(seurat_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

FeaturePlot(seurat_merged, features = c("PPBP", "PF4"))

install.packages("org.Hs.eg.db")
library(org.Hs.eg.db)
# Assume current genes are Ensembl Gene IDs
current_genes <- rownames(seurat_merged)
# Convert to gene symbols
gene_names <- mapIds(org.Hs.eg.db, keys = current_genes, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
# Update gene names in the Seurat object
rownames(seurat_merged) <- gene_names

expr_matrix <- GetAssayData(seurat_merged, assay = "RNA", slot = "data")
expr_matrix <- as.matrix(expr_matrix)
ref_data <- MonacoImmuneData()
main_predictions <- SingleR(
  test = expr_matrix,        # Matrix to annotate (from Seurat)
  ref = ref_data,            # Reference dataset
  labels = ref_data$label.main,  # Reference cell-type labels (e.g., "T cell", "B cell")
  clusters = seurat_merged$seurat_clusters  # Optional: annotate by Seurat cluster (improves stability)
)
fine_predictions <- SingleR(
  test = expr_matrix,        # Matrix to annotate (from Seurat)
  ref = ref_data,            # Reference dataset
  labels = ref_data$label.fine,  # Reference cell-type labels (e.g., "T cell", "B cell")
  clusters = seurat_merged$seurat_clusters  # Optional: annotate by Seurat cluster (improves stability)
)
# View annotation results
table(main_predictions$labels)
table(fine_predictions$labels)

# Add annotation results as new columns in the Seurat object
# At cell level
seurat_merged$singleR_main_labels <- main_predictions$labels
seurat_merged$singleR_fine_labels <- fine_predictions$labels
# At cluster level
seurat_merged$singleR_cluster_main_labels <- main_predictions$labels[seurat_merged$seurat_clusters]
seurat_merged$singleR_cluster_fine_labels <- fine_predictions$labels[seurat_merged$seurat_clusters]

# Set annotation results as cell Idents
seurat_merged <- SetIdent(seurat_merged, value = "singleR_cluster_main_labels")
seurat_merged <- SetIdent(seurat_merged, value = "singleR_cluster_fine_labels")

table(seurat_merged$seurat_clusters, seurat_merged$singleR_cluster_main_labels)

cluster_annotation <- c("T cells", "B cells", "Monocytes", "NK cells", "Monocytes",
                        "Platelet", "Dendritic cells", "Neutrophil", "Platelet", "Dendritic cells",
                        "NK-T cells", "Monocytes", "Platelet", "Plasma cells", "Dendritic cells")
# cluster_annotation is a character vector with one annotation per cluster
# Its length must match the number of seurat_clusters categories
seurat_merged$manual_annotation <- cluster_annotation[as.numeric(seurat_merged$seurat_clusters)]

seurat_merged <- SetIdent(seurat_merged, value = "manual_annotation")

Idents(seurat_merged) <- "manual_annotation"

DotPlot(seurat_merged, features = c(
  "PPBP", "PF4", #Platelet
  "CD3D", "CD3E", "CD4", "CD8A", "AQP3", "LEF1", "TCF7", "CXCR3", "CCL4", "CCR6", #T-cells
  "CD38", "MZB1", "IGKC", "SDC1", #Plasma cells
  "MS4A1", "CD19", "CD27", "TCL1A", #B-cells
  "NCAM1", "FCGR3A", "NKG7", "GNLY", #NK cells
  "CD14", "CD33", "LYZ" #Myeloid cells
)) + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

DimPlot(seurat_merged, reduction = "umap", label = TRUE)

DotPlot(seurat_merged, features = c("IL15RA", "CARD19", "IFNAR2", "LY6E", "TRBV3-1", "FHL1"
),
group.by = "group"
) + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(seurat_merged, features = c("IL15RA", "CARD19", "IFNAR2", "LY6E", "TRBV3-1", "FHL1"
)) + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Platelet-only analysis

allcells <- seurat_merged

print(allcells@reductions)

allcells@reductions$pca <- NULL
allcells@reductions$umap <- NULL

Plt_cells_public <- subset(allcells, subset = manual_annotation == "Platelet")
saveRDS(Plt_cells_public,
        file = "path/to/your/public/data/plt-public.rds",
        compress = "gzip")
Plt_cells_public_fixed <- readRDS("path/to/your/public/data/plt-public.rds")

# -- 1. Extract raw count data from your existing object
#    The default Assay is typically "RNA"
raw_counts <- GetAssayData(Plt_cells_public, assay = "RNA", slot = "counts")

# -- 2. Extract metadata
metadata <- Plt_cells_public@meta.data

# -- 3. Get current (possibly duplicated) row names
original_rownames <- rownames(raw_counts)

# -- 4. Use make.unique() to create new, unique row names
#    E.g., two "MALAT1" entries become "MALAT1" and "MALAT1.1"
unique_rownames <- make.unique(original_rownames)

# -- 5. Apply the new unique row names to the raw count matrix
rownames(raw_counts) <- unique_rownames

# -- 6. Rebuild a fully functional Seurat object with the fixed data
#    This is the safest approach, ensuring complete internal consistency
Plt_cells_public_fixed <- CreateSeuratObject(counts = raw_counts, meta.data = metadata)

# -- 7. Verify no duplicate row names remain (should return FALSE)
any(duplicated(rownames(Plt_cells_public_fixed)))


print(Plt_cells_public_fixed)

Plt_cells_public_fixed <- NormalizeData(Plt_cells_public_fixed)
Plt_cells_public_fixed <- FindVariableFeatures(Plt_cells_public_fixed)

Plt_cells_public_fixed <- ScaleData(Plt_cells_public_fixed)
Plt_cells_public_fixed <- RunPCA(Plt_cells_public_fixed, npcs = 30)

Plt_cells_public_fixed <- FindNeighbors(Plt_cells_public_fixed, reduction = "pca", dims = 1:12, k.param = 10)
Plt_cells_public_fixed <- FindClusters(Plt_cells_public_fixed, resolution = 0.2, algorithm = 3)

Plt_cells_public_fixed <- RunUMAP(Plt_cells_public_fixed, dims = 1:10) # Adjust dims for your data

#Plt_cells_public_fixed <- SetIdent(Plt_cells_public_fixed, value = "seurat_clusters")


DimPlot(Plt_cells_public_fixed, reduction = "umap", label = TRUE)

FeaturePlot(Plt_cells_public_fixed, features = c("PPBP", "PF4", "ITGA2B", "ANO6", "HIST1H2BO", "HIST1H4H"))

all.markers <- FindAllMarkers(Plt_cells_public_fixed, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
all_markers_top5 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(Plt_cells_public_fixed, features = all_markers_top5$gene) + NoLegend()

cluster_to_plt_type <- c(
  "Type A",
  "Type A",
  "Type B",
  "Type B",
  "Type A",
  "Type A",
  "Type C",
  "Type C",
  "Type D")
Plt_cells_public_fixed$Manual_plt_type <- cluster_to_plt_type[as.numeric(as.character(Plt_cells_public_fixed$seurat_clusters)) + 1]
Plt_cells_public_fixed <- SetIdent(Plt_cells_public_fixed, value = "Manual_plt_type")
DimPlot(Plt_cells_public_fixed, reduction = "umap", label = TRUE)

all.markers_new <- FindAllMarkers(Plt_cells_public_fixed, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers_new %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
all_markers_top5_new <- all.markers_new %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(Plt_cells_public_fixed, features = all_markers_top5_new$gene) + NoLegend()


# FHL1-related plots
FeaturePlot(Plt_cells_public_fixed, features = c("FHL1"))
VlnPlot(Plt_cells_public_fixed, features = c("FHL1"), split.by = "group")
VlnPlot(Plt_cells_public_fixed, features = c("FHL1"))
DimPlot(Plt_cells_public_fixed, reduction = "umap", label = TRUE, split.by = "group")

# Count cells per group and per cell type
cell_type_counts <- as.data.frame(table(Plt_cells_public_fixed$group, Plt_cells_public_fixed$Manual_plt_type))
# Calculate proportion of each cell type within each group
cell_type_counts <- cell_type_counts %>%
  group_by(Var1) %>%  # Var1 is group_column
  mutate(Percent = Freq / sum(Freq) * 100)
# Stacked bar plot
ggplot(cell_type_counts, aes(x = Var1, y = Percent, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar chart
  scale_fill_manual(values = c("brown", "orange", "green", "purple")) +  # Custom colors
  labs(x = "Group", y = "Cell type proportion (%)", fill = "Cell type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Build contingency table for different cell types
contingency_table <- table(Plt_cells_public_fixed$Manual_plt_type %in% c("Type A", "Type B"), Plt_cells_public_fixed$group)
rownames(contingency_table) <- c("Activated Platelets", "Other Platelets")
#colnames(contingency_table) <- c("Control", "Treatment")
contingency_table
chisq_test <- chisq.test(contingency_table)
chisq_test
fisher_test <- fisher.test(contingency_table)
fisher_test


DotPlot(seurat_merged, features = c("PPBP", "PF4", "ITGA2B", "ANO6", "HIST1H2BO", "HIST1H4H"
),
group.by = "group"
) + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(seurat_merged, features = c("IL15RA", "CARD19", "IFNAR2", "LY6E", "TRBV3-1", "FHL1"
)) + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))


