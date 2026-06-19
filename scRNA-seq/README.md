# scRNA-seq Joint Analysis Toolkit

R scripts for single-cell RNA sequencing data processing, cell-type annotation, differential
expression analysis, intercellular communication, and publication-ready visualization.

---

## Dependencies

| Package         | Version (tested) | Purpose                                   |
| --------------- | ---------------- | ----------------------------------------- |
| R               | >= 4.2           | Runtime                                   |
| Seurat          | >= 5.0           | scRNA-seq data handling, clustering, DEG  |
| monocle3        | >= 1.3           | Trajectory / pseudotime analysis          |
| clusterProfiler | >= 4.0           | GO / KEGG enrichment                      |
| CellChat        | >= 1.6           | Ligand–receptor interaction inference     |
| SingleR         | >= 2.0           | Reference-based cell-type annotation      |
| celldex         | —                | Reference expression datasets for SingleR |
| ggplot2         | —                | Plotting engine                           |
| dplyr           | —                | Data manipulation                         |
| reshape2        | —                | Data reshaping (melt / cast)              |
| NMF             | —                | Non-negative matrix factorization         |
| circlize        | —                | Circular visualization (CellChat dep.)    |
| ComplexHeatmap  | —                | Advanced heatmaps (CellChat dep.)         |
| sf              | —                | Spatial features (monocle3 dep.)          |
| devtools        | —                | GitHub package installation               |
| BiocManager     | —                | Bioconductor package management           |

### Install

Run `SingcellAnalysis-packages_install.R` once to set up the environment:

```r
source("SingcellAnalysis-packages_install.R")
```

Alternatively, install step-by-step:

```r
# Core
install.packages('Seurat')

# Bioconductor ecosystem
BiocManager::install(c('clusterProfiler', 'SingleR', 'celldex',
                        'BiocGenerics', 'DelayedArray', 'SingleCellExperiment',
                        'SummarizedExperiment'))

# From GitHub
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("jinworks/CellChat")

# Utilities
install.packages(c('NMF', 'sf'))
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
```

---

## Project Structure

```
project/
│
├── data/                              # Input Seurat objects (.rds) and CSV expression tables
│   ├── own_dataset.rds                # Your in-house scRNA-seq data
│   └── public_dataset.rds             # Downloaded public dataset
│
├── results/                           # Output figures and CSV tables (auto-generated)
│
├── exp-OUT.R                          # Expression output — compute mean, proportion & significance
├── Box-bar plotting.R                 # Box–bar plots for selected genes across cell types
├── Dataset Modification.R             # Harmonise idents, reorder groups & produce dot/UMAP plots
├── pub-data-processing.R              # Process public dataset (merge, annotate, subset)
├── SingcellAnalysis-packages_install.R# One-shot package installation
└── README.md
```

---

## Scripts

### 1. `SingcellAnalysis-packages_install.R` — Environment Setup

Installs and loads all required R packages in one go. Run this first on a fresh machine.

```r
source("SingcellAnalysis-packages_install.R")
```

- Installs Seurat, monocle3, clusterProfiler, CellChat, SingleR and their dependencies
- Uses `BiocManager` for Bioconductor packages and `devtools` for GitHub sources
- Loads each package after installation

### 2. `exp-OUT.R` — Expression Statistics & Significance Testing

Computes per-gene expression summaries across cell types and groups, then runs
pairwise significance tests.

```r
source("exp-OUT.R")
```

- **User-configurable section** at the top: set your Seurat object, gene list, group column and output prefix
- Computes mean, median, standard deviation, and quantiles (Q5, Q25, Q75, Q95) per cell-type × group
- Calculates expression proportion (fraction of cells with expression > 0)
- Performs Wilcoxon rank-sum and Fisher's exact tests for every group pair
- Applies Benjamini–Hochberg correction to all p-values
- Saves three CSV files to the specified `data_dir`:
  - `*_average_expression_YYYYMMDD.csv`
  - `*_expression_proportion_YYYYMMDD.csv`
  - `*_significance_test_YYYYMMDD.csv`

### 3. `Box-bar plotting.R` — Gene Expression Box–Bar Plots

Generates combined box-and-bar charts showing per-cell-type expression of one
or more genes, faceted by group.

```r
source("Box-bar plotting.R")
```

- Reads a pre-computed average-expression CSV (output of `exp-OUT.R`)
- Dynamically filters gene–cell-type combinations with mean expression above a threshold
- Aligns cell-type order across facets
- Uses `ggplot2` with `geom_boxplot` + `geom_bar` (mean ± SD)
- Auto-assigns colour palette based on the number of groups in the data

### 4. `Dataset Modification.R` — Identity Harmonisation & Visualisation

Aligns cell-type identities (`Idents`) between your own dataset and a public
reference, then generates comparative dot plots and UMAP projections.

```r
source("Dataset Modification.R")
```

- Loads both own and public Seurat objects
- **Recode idents** — maps categorical labels to a unified ontology (e.g. B cell → B, CD8+ T → CD8 T, etc.)
- **Reorder x-axis** — ensures consistent cell-type order in all downstream plots
- Generates four figures:
  - Own dataset dot plot (marker expression × cell type)
  - Public dataset dot plot
  - Own dataset UMAP coloured by cell type
  - Public dataset UMAP coloured by cell type
- Supports threshold-based cell filtering (e.g. keep only cell types with ≥ N cells)

### 5. `pub-data-processing.R` — Public Dataset Pipeline

Full processing chain for a public scRNA-seq dataset: read multiple samples,
merge, normalise, annotate cell types, and extract subpopulations.

```r
source("pub-data-processing.R")
```

- Reads individual sample directories (10X-style `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`)
- Creates Seurat objects and merges them into a single combined object
- Performs QC filtering (`nFeature_RNA`, `percent.mt`), normalisation (`LogNormalize`),
  variable-feature selection, scaling, PCA, and UMAP
- Runs **SingleR** for reference-based cell-type annotation using `celldex` references
  (e.g. `HumanPrimaryCellAtlasData`, `BlueprintEncodeData`)
- Generates a cell-type × group table and a UMAP coloured by annotation
- Extracts platelet subpopulation for downstream focused analysis
- Saves the annotated Seurat object as `plt-public.rds`

---

## Typical Workflow

1. **Install** — run `SingcellAnalysis-packages_install.R` once.
2. **Prepare data** — place your `.rds` files in `data/`.
3. **Process public data** — run `pub-data-processing.R` to merge, annotate and subset.
4. **Harmonise idents** — run `Dataset Modification.R` to align cell-type labels.
5. **Expression analysis** — run `exp-OUT.R` to compute statistics and significance.
6. **Visualise** — run `Box-bar plotting.R` to produce publication figures.

---

## Notes

- All paths containing personal identifiers have been replaced with `path/to/your/…`
  placeholders. Update them to match your local directory structure before execution.
- The `exp-OUT.R` script has a clearly marked **user-modifiable section** at the top;
  the remainder is designed to run without further editing.
- `Dataset Modification.R` assumes a one-to-one cell-type mapping between your
  own and public datasets. Adjust the `recode` mappings if the ontologies differ.
- SingleR annotation in `pub-data-processing.R` uses multiple references in cascade;
  change the reference list in the script to suit your tissue or species.
- All scripts use base R or Bioconductor conventions — no hard-coded absolute paths
  remain in the distributed version.
