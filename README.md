# Bioinformatics Analysis Code for Graduation Thesis

This repository contains the data analysis code used in the bioinformatics section of my graduation thesis.

**Research on Susceptibility Genes for Severe Community Acquired Pneumonia in the Chinese Population Based on the Human Whole Genome**

## Overview

The scripts and related files in this repository are associated with the data analysis part of my thesis project. They are intended to document the computational workflow, improve reproducibility, and provide a public record of the analyses performed during the study.

Please note that this repository is currently under organization. I will continue to clean, annotate, and upload the relevant scripts over the coming period.

## Data Availability

The data used in this study are available from the following sources:

- Part of the original data has been submitted to the university archives.
- **Note:** Data files in this repository (under `data/`) are example datasets only. The complete clinical data are archived at Peking University and are not publicly included here.
- Public database data were obtained from: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5102902`.
- Single-cell sequencing data have been uploaded to GEO database: `GSE262512`.
- Whole-genome sequencing data have been archived at Peking University and are available through the university's data access procedures.

## Repository Structure

```
.
├── WGS/                            # Whole-genome sequencing analysis ★
│   ├── preprocessing/              # Data preprocessing & quality control
│   ├── mapping/                    # Mapping, deduplication & QC statistics
│   ├── variant-calling/            # Variant calling (GATK & DeepVariant)
│   ├── annotation/                 # Variant annotation (SnpEff / VEP / ANNOVAR)
│   ├── variant-comparison/         # Two-cohort variant comparison (chi-square / Fisher)
│   ├── association/                # PLINK association analysis
│   ├── visualization/              # Manhattan / QQ / circular plots, IGV inspection
│   ├── figures/                    # Downstream thesis figure scripts (to be added)
│   └── README.md
├── scRNA-seq/                      # Single-cell RNA sequencing analysis ★
│   ├── data/                       # Input Seurat objects (.rds) and CSV expression tables
│   ├── results/                    # Output figures and CSV tables (auto-generated)
│   ├── exp-OUT.R                   # Expression output — compute mean, proportion & significance
│   ├── Box-bar plotting.R          # Box–bar plots for selected genes across cell types
│   ├── Dataset Modification.R      # Harmonise idents, reorder groups & produce dot/UMAP plots
│   ├── pub-data-processing.R       # Process public dataset (merge, annotate, subset)
│   ├── SingcellAnalysis-packages_install.R  # One-shot package installation
│   └── README.md
├── SCAP-prediction/                # Genotype-panel screening & SCAP prediction models ★
│   ├── convert_to_onehot.py        # Genotype table (-1/0/1/2) → one-hot features
│   ├── lasso_selection_glmnet.R    # LASSO (glmnet, lambda.1se) feature screening + figures
│   ├── split_train_test.py         # Stratified 80/20 train/test split
│   ├── consensus_feature_selection.py  # RF / SVM-RFE / XGBoost+SHAP consensus screening
│   ├── train_final_model.py        # Final XGBoost + logistic-regression panel (ROC / calibration / DCA / nomogram / SHAP)
│   ├── predict_new_samples.py      # Inference on new patients from saved model assets
│   └── plot_confusion_corr.py      # Confusion matrix & feature-correlation heatmaps
├── LPCAT-Machine_Learning/         # DNN classification on clinical data ★
│   ├── data/                       # Example CSV datasets (full data archived at Peking University)
│   ├── models/                     # Saved Keras models (auto-generated)
│   ├── results/roc/                # ROC curve output images
│   ├── Machine_Learning.py         # Train a model from scratch
│   ├── model_performance.py        # Evaluate a saved model + single ROC
│   ├── ROCs.py                     # Compare ROC across dataset variants
│   ├── Self_Defining_Function.py   # Shared prediction utilities
│   └── README.md
└── README.md
```

The final structure may be adjusted as additional scripts are organized.

---

## Component Details

### 1. WGS — Whole-Genome Sequencing Analysis ★

Shell / Python / R scripts for the complete WGS pipeline of the SCAP cohort: starting from the raw data (.fq.gz) provided by the sequencing company, through alignment against **three reference genomes in parallel** (hg38 / GRCh38, T2T-CHM13, T2T-YAO), variant calling with two independent callers (GATK and DeepVariant), functional annotation with three annotators (SnpEff / VEP / ANNOVAR), case–control variant comparison, PLINK association analysis, and publication-ready visualization.

All scripts were executed on an offline Linux compute server. Server paths, usernames, project identifiers, and sample IDs in the released scripts have been replaced with placeholders (`/path/to/...`, `<sample_id>`, `CAP001` as an example) — edit the user-configurable section at the top of each script before running.

#### Dependencies

| Tool / Package | Version (used) | Purpose                                          |
|----------------|----------------|--------------------------------------------------|
| fastp          | —              | Raw read filtering (parameters set by vendor)    |
| bwa (bwa-mem)  | —              | Read alignment                                   |
| samtools       | 1.x            | Sorting, fixmate, markdup, indexing, stats       |
| GATK           | 4.5.0.0        | HaplotypeCaller / CombineGVCFs / GenotypeGVCFs / VariantFiltration / SelectVariants / MergeVcfs |
| DeepVariant    | 1.6.1 (GPU Docker) | Alternative variant caller                   |
| GLnexus        | —              | Joint-calling of DeepVariant GVCFs               |
| bcftools       | —              | VCF filtering / merging / stats                  |
| vcftools       | —              | SNP / INDEL separation and counting (early stage)|
| bgzip / tabix  | —              | Compression and indexing of (g)VCF              |
| SnpEff         | —              | Variant annotation (ANN field)                   |
| VEP            | Docker         | Variant annotation (CSQ field)                   |
| ANNOVAR        | —              | Variant annotation (refGene databases)           |
| gff3ToGenePred | UCSC binary    | GFF3 → GenePred conversion for ANNOVAR databases |
| PLINK / PLINK2 | 1.9 / 2.x      | Association analysis and QC                      |
| bedtools       | —              | Reference sequence extraction                    |
| Python         | ≥ 3.8          | Runtime for extraction / comparison scripts      |
| pysam          | —              | VCF parsing                                       |
| numpy / scipy  | —              | Contingency-table tests                           |
| pandas         | —              | Table manipulation                                |
| Biopython      | —              | FASTA parsing (chromosome lengths)                |
| R (qqman, CMplot) | —          | Manhattan / QQ / circular plots                   |
| IGV            | 2.17.2         | Manual BAM / VCF inspection                       |

Environment setup on the offline server: conda environments were built in a local Linux VM (`wgs` for pysam/numpy/scipy/pandas; `wgs_plot` for plotting; an R environment for karyoploteR), then packaged and uploaded to the server.

#### Directory Layout & Scripts

**`preprocessing/` — Data Preprocessing & Quality Control**

- **`fastp_clean_data.sh`** — filter vendor-provided paired-end raw reads with fastp; filtering parameters (`-g -q 5 -u 50 -n 15 -l 150`, overlap limits) follow the sequencing company's specification and must stay consistent across all samples. Outputs cleaned `.fq.gz` plus an HTML/JSON QC report per sample.
- **`filter_hg38_primary_chroms.sh`** — keep only primary chromosomes (chr1–22, chrX, chrY, chrM) from the hg38 FASTA and GFF3; the UCSC hg38 contains hundreds of unplaced/alt contigs that otherwise break downstream per-chromosome steps.
- **`extract_bamstat_info.py`** — pull key metrics (total reads, mapped, properly paired, …) from `samtools stats` reports of every sample × reference combination into one CSV for mapping-quality comparison. Reads fixed line numbers of the samtools stats layout — verify against your samtools version.
- **`extract_vcfstat_counts.py`** — extract SNP / INDEL counts from bcftools-stats-style log files into a two-column CSV.

**`mapping/` — Alignment & BAM Processing**

- **`bwa_index.sh`** — build bwa indexes for the three reference genomes (once).
- **`bwa_alignment.sh`** — align paired-end clean reads with `bwa mem` (RG mandatory — GATK rejects BAMs without it), piped into `samtools sort -n` (name order required by the next step).
- **`samtools_fixmate.sh`** — fill mate coordinates / ISIZE / flags (prerequisite for markdup). Supports both name-sorted input (direct call) and coordinate-sorted input (re-sort pipe).
- **`samtools_markdup.sh`** — coordinate-sort then remove PCR duplicates with `samtools markdup`.
- **`samtools_index_stats.sh`** — index the dedup BAMs and generate `samtools stats` reports consumed by `extract_bamstat_info.py`.
- **`samtools_add_rg.sh`** — rescue script re-adding RG to BAMs aligned without `-R`.
- **`workflow_fastq_to_dedup_bam.sh`** — fully automated FASTQ → sort → fixmate → markdup → index → stats workflow, one call per sample, three references in parallel; supports sample-ID standardisation (`old_name` → `new_name`).
- **`md5sum_bam_archive.sh`** — generate md5 checksums of the pre-variant-calling BAM archive for transfer verification.

**`variant-calling/` — Variant Calling (GATK + DeepVariant)**

- **`gatk_haplotypecaller_gvcf.sh`** — per-sample `HaplotypeCaller -ERC GVCF` + bgzip + tabix, three references in parallel. GVCF mode is used because merging plain VCFs cannot distinguish `./.` (not called) from `0/0` (reference) genotypes.
- **`gatk_combine_gvcfs.sh`** — hierarchical `CombineGVCFs` merging (samples → lv1 parts → lv2 parts → cohort final); whole-cohort GVCFs are 60–70 GB each, so merging is levelled. Reference genome chosen via indirect variable expansion.
- **`gatk_gvcf_to_vcf.sh`** — `GenotypeGVCFs` joint genotyping of a combined cohort GVCF into a multi-sample VCF.
- **`gatk_hardfilter_vcf.sh`** — hard filtering: split SNP / INDEL (SelectVariants) → filter with GATK-recommended thresholds (SNP: QD<2, QUAL<30, SOR>3, FS>60, MQ<40, MQRankSum<-12.5, ReadPosRankSum<-8; INDEL: QD<2, QUAL<30, FS>200, ReadPosRankSum<-20) → keep PASS only → MergeVcfs. Separate single-expression filters (not compound `||`) so sites lacking MQRankSum/ReadPosRankSum annotations are not silently dropped.
- **`gatk_split_gvcf_by_chrom.sh`** — merge two cohort GVCFs chromosome by chromosome (SelectVariants -L → CombineGVCFs → GenotypeGVCFs per chromosome, all chromosomes in parallel) to avoid memory/storage blow-up; edit the chromosome list per reference (CHM13: no chrM; hg38: extract contig list from FASTA header first).
- **`gatk_merge_chr_vcfs.sh`** — merge the per-chromosome VCFs into one cohort VCF.
- **`deepvariant_batch.sh`** — batch DeepVariant 1.6.1 (GPU Docker) calling driven by a sample CSV (`SAMPLE_ID,REF_GENOME,SEX`; male samples get `--haploid_contigs=chrX,chrY`). Pins one container per GPU.
- **`glnexus_merge_gvcfs.sh`** — joint-call DeepVariant GVCFs with GLnexus (equivalent to CombineGVCFs + GenotypeGVCFs, much faster); supports `DeepVariantWGS` and `DeepVariant_unfiltered` configs; each run needs its own temp directory.

**`annotation/` — Functional Annotation**

- **`snpeff_anno_clean.sh`** — SnpEff annotation of filtered VCFs, then remove WARNING lines SnpEff mixes into the redirected VCF (they corrupt tabix/bcftools parsing), then bgzip + tabix.
- **`vep_anno.sh`** — VEP annotation via Docker with local FASTA + GFF3 (fully offline). With `--fasta`/`--gff` set, VEP is already offline — do **not** add `--offline` (it triggers a cache check that fails).
- **`annovar_build_db.sh`** — build custom ANNOVAR refGene databases from reference GFF3 + FASTA (gff3ToGenePred → GenePred; retrieve_seq_from_fasta.pl → transcript FASTA).
- **`annovar_annotate.sh`** — ANNOVAR `table_annovar.pl` gene-based annotation of filtered VCFs; reference database inferred from the VCF file name.
- **`fix_gff3_order.sh`** — re-sort a GFF3 by chromosome/start (tabix refuses unsorted GFF3); used for the T2T-CHM13 annotation.
- **`fix_gff3_parent_refs.py`** — repair invalid `Parent` references in the raw T2T-YAO GFF3 (records kept as documentation of the repair attempt; the final YAO database was built from an AGAT-cleaned GFF3).

**`variant-comparison/` — Case–Control Variant Comparison ★**

- **`vcf_info_extraction_snpeff.py`** — compare two annotated cohort VCFs (NSCP vs SCAP, same reference) variant by variant: build a 2×2 contingency table from AC/AN allele counts (variants unique to one group included, missing group AC=0), run a chi-square test per variant and switch to Fisher's exact test when any observed count in the 2×2 table is < 5, all variants processed in parallel (ThreadPoolExecutor). Outputs a combined CSV, per-chromosome CSVs (primary chromosomes only), HIGH/MODERATE-impact CSVs with gene names parsed from the SnpEff ANN field, and two VCFs of group-unique variants.
- **`vcf_info_extraction_annovar.py`** — ANNOVAR counterpart; instead of impact tiers it extracts variants with exonic function changes (ExonicFunc.refGene) and amino-acid changes (AAChange.refGene) into dedicated CSVs.
- **`vcf_info_extraction_vep.py`** — VEP counterpart; parses the CSQ field and additionally reports the MODIFIER impact tier in its own CSV.
- **`count_significant_variants.sh`** — count rows with P < 0.05 in a comparison output CSV.
- **`bedtools_extract_refseq.sh`** — extract reference sequences (±200/±500 bp) around candidate variants with `bedtools getfasta`.

**`association/` — PLINK Association Analysis ★**

- **`vcf_to_plink_bed.sh`** — prepare cohort VCFs for PLINK: keep biallelic SNPs only (`bcftools view -m2 -M2 -v snps`), merge NSCP + SCAP cohorts per reference (`bcftools merge`), convert to bed/bim/fam with `plink2 --make-bed` (sex supplied via `--psam`; chrX PAR boundaries via `--split-par`; `--allow-extra-chr` for hg38; `--vcf-half-call m` for DeepVariant VCFs).
- **`csv_to_psam.py`** — convert a two-column sample/sex CSV into a PLINK2 `.psam` file.
- **`bim_add_variant_id.py`** — assign each variant a unique deterministic ID (`chr:pos:REF:ALT`) by rewriting the .bim column 2, so association results are unambiguously addressable.
- **`plink_qc.sh`** — PLINK QC (`--mind 0.05 --geno 0.05 --maf 0.01 --hwe 1e-6`; thresholds relaxed to 0.60 for very-low-frequency variant sets — document the choice per dataset).
- **`plink_assoc.sh`** — PLINK 1.9 case/control association (`--assoc` / `--logistic`, each with `--adjust` for Bonferroni / Holm / Sidak / FDR_BH / FDR_BY corrected P values), covering GATK- and DeepVariant-derived bed files, QC-filtered and unfiltered.
- **`extract_significant_variants.sh`** — extract P < 0.05 variants from `.assoc` results → CSV → BED → pull the matching records out of the annotated VCF with `bcftools view -R`.
- **`extract_impact_info.py`** — extract HIGH / MODERATE impact variants (with gene name and variant type) from a SnpEff-annotated significant-variant VCF into readable CSVs.

**`visualization/` — Visualization & Manual Inspection**

- **`manhattan_qq.R`** — Manhattan and QQ plots from PLINK `.assoc` results with qqman (cleans chromosome codes, removes PAR rows and invalid P values).
- **`circular_manhattan.R`** — circular Manhattan plot with CMplot (≤ 1,000,000 variants — downsample larger result sets first).
- **`extract_chromosome_lengths.py`** — chromosome lengths from a reference FASTA into a CSV (input for karyoploteR / RIdeogram chromosome ideograms).
- **`extract_gene_positions.py`** — look up candidate-gene coordinates across the three reference GFF3s into one wide CSV for cross-reference gene marking.
- **`extract_regions_for_igv.sh`** — slice BAMs (`samtools view -L`) and the annotated VCF (`bcftools view -R`) around candidate sites, merge per cohort, for manual read-coverage / genotype inspection in IGV.

**`figures/` — Downstream Thesis Figures (to be added)**

Reserved for the figure-generation scripts used to produce the final WGS figures of the thesis (cohort comparison summaries, per-chromosome statistics, combined multi-panel figures). Scripts will be uploaded as the repository is updated.

#### Typical Workflow

1. **Preprocess** — `fastp_clean_data.sh` on vendor raw data; run `filter_hg38_primary_chroms.sh` once for hg38.
2. **Align & dedup** — `workflow_fastq_to_dedup_bam.sh <old_name> <new_name>` per sample (or the individual mapping scripts); verify mapping stats with `extract_bamstat_info.py`.
3. **Call variants (GATK)** — `gatk_haplotypecaller_gvcf.sh <sample>` per sample → `gatk_combine_gvcfs.sh` level by level → `gatk_gvcf_to_vcf.sh <cohort> <ref>` → `gatk_hardfilter_vcf.sh <cohort> <ref>` (or the chromosome-split route: `gatk_split_gvcf_by_chrom.sh` + `gatk_merge_chr_vcfs.sh`).
4. **Call variants (DeepVariant)** — `deepvariant_batch.sh samples.csv` per GPU → `glnexus_merge_gvcfs.sh`.
5. **Annotate** — `snpeff_anno_clean.sh` / `vep_anno.sh` / `annovar_annotate.sh` (build the SnpEff / ANNOVAR databases first if needed).
6. **Compare cohorts** — `vcf_info_extraction_snpeff.py` (or the ANNOVAR / VEP variants) for per-variant chi-square / Fisher tests; `count_significant_variants.sh` for a quick tally.
7. **Associate** — `vcf_to_plink_bed.sh` → `plink_qc.sh` → `plink_assoc.sh` → `extract_significant_variants.sh` → `extract_impact_info.py`.
8. **Visualise & verify** — `manhattan_qq.R` / `circular_manhattan.R`; `extract_regions_for_igv.sh` for manual IGV inspection of candidate sites.

#### Technical Notes

- **Three references, two callers, three annotators.** Every analysis step was run in parallel for hg38 / T2T-CHM13 / T2T-YAO, and variant calling was performed independently with GATK and DeepVariant; results were cross-validated (see the PLINK section for the parallel result sets).
- **GVCF over VCF for merging** — plain VCF merging cannot distinguish `./.` from `0/0` at multi-sample level, biasing cohort comparisons; hence the GVCF route.
- **Hard-filter expressions are single, not compound** — sites without MQRankSum / ReadPosRankSum annotations (non-heterozygous sites) would be dropped by a compound `||` filter.
- **SnpEff stdout pollution** — redirecting SnpEff output mixes WARNING lines into the VCF, breaking tabix/bcftools; the `awk '!/^WARNING/'` cleanup in `snpeff_anno_clean.sh` is mandatory.
- **VEP offline mode** — with `--fasta` + `--gff` supplied, adding `--offline` triggers a failing cache check; omit it.
- **DeepVariant GPU pinning** — parallel containers without `--gpus '"device=N"'` caused GPU out-of-memory crashes; run one container per GPU, one script per device.
- **GLnexus temp directories** — concurrent `glnexus_cli` runs collide in `.GLnexus.DB`; give each run its own working directory.
- **Long-running jobs** — run scripts inside `screen` or under `nohup`; SSH disconnects otherwise kill GATK workflows mid-run.
- **Version consistency** — GATK was upgraded from 4.0.5.1 to 4.5.0.0 mid-study; earlier gVCFs were re-generated so the whole cohort used one GATK version. Keep caller versions uniform across a cohort.
- **Line-number-based extraction scripts** (`extract_bamstat_info.py`, `extract_vcfstat_counts.py`) assume fixed report layouts; re-verify the offsets against your tool versions before running.

### 2. scRNA-seq — Single-Cell RNA Sequencing Analysis ★

R scripts for single-cell RNA sequencing data processing, cell-type annotation, differential expression analysis, intercellular communication, and publication-ready visualization.

#### Dependencies

| Package         | Version (tested) | Purpose                                   |
|-----------------|------------------|-------------------------------------------|
| R               | ≥ 4.2            | Runtime                                   |
| Seurat          | ≥ 5.0            | scRNA-seq data handling, clustering, DEG  |
| monocle3        | ≥ 1.3            | Trajectory / pseudotime analysis          |
| clusterProfiler | ≥ 4.0            | GO / KEGG enrichment                      |
| CellChat        | ≥ 1.6            | Ligand–receptor interaction inference     |
| SingleR         | ≥ 2.0            | Reference-based cell-type annotation      |
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

Run `SingcellAnalysis-packages_install.R` once to set up the environment:

```r
source("SingcellAnalysis-packages_install.R")
```

#### Scripts

**`SingcellAnalysis-packages_install.R` — Environment Setup**

Installs and loads all required R packages in one go. Run this first on a fresh machine.

```r
source("SingcellAnalysis-packages_install.R")
```

- Installs Seurat, monocle3, clusterProfiler, CellChat, SingleR and their dependencies
- Uses `BiocManager` for Bioconductor packages and `devtools` for GitHub sources
- Loads each package after installation

**`exp-OUT.R` — Expression Statistics & Significance Testing**

Computes per-gene expression summaries across cell types and groups, then runs pairwise significance tests.

```r
source("exp-OUT.R")
```

- **User-configurable section** at the top: set your Seurat object, gene list, group column and output prefix
- Computes mean, median, standard deviation, and quantiles (Q5, Q25, Q75, Q95) per cell-type × group
- Calculates expression proportion (fraction of cells with expression > 0)
- Performs Wilcoxon rank-sum and Fisher's exact tests for every group pair
- Applies Benjamini–Hochberg correction to all p-values
- Saves three CSV files to the specified `data/` directory:
  - `*_average_expression_YYYYMMDD.csv`
  - `*_expression_proportion_YYYYMMDD.csv`
  - `*_significance_test_YYYYMMDD.csv`

**`Box-bar plotting.R` — Gene Expression Box–Bar Plots**

Generates combined box-and-bar charts showing per-cell-type expression of one or more genes, faceted by group.

```r
source("Box-bar plotting.R")
```

- Reads a pre-computed average-expression CSV (output of `exp-OUT.R`)
- Dynamically filters gene–cell-type combinations with mean expression above a threshold
- Aligns cell-type order across facets
- Uses `ggplot2` with `geom_boxplot` + `geom_bar` (mean ± SD)
- Auto-assigns colour palette based on the number of groups in the data

**`Dataset Modification.R` — Identity Harmonisation & Visualisation**

Aligns cell-type identities (`Idents`) between your own dataset and a public reference, then generates comparative dot plots and UMAP projections.

```r
source("Dataset Modification.R")
```

- Loads both own and public Seurat objects
- **Recode idents** — maps categorical labels to a unified ontology (e.g. B cell → B, CD8+ T → CD8 T)
- **Reorder x-axis** — ensures consistent cell-type order in all downstream plots
- Generates four figures:
  - Own dataset dot plot (marker expression × cell type)
  - Public dataset dot plot
  - Own dataset UMAP coloured by cell type
  - Public dataset UMAP coloured by cell type
- Supports threshold-based cell filtering (e.g. keep only cell types with ≥ N cells)

**`pub-data-processing.R` — Public Dataset Pipeline**

Full processing chain for a public scRNA-seq dataset: read multiple samples, merge, normalise, annotate cell types, and extract subpopulations.

```r
source("pub-data-processing.R")
```

- Reads individual sample directories (10X-style `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`)
- Creates Seurat objects and merges them into a single combined object
- Performs QC filtering (`nFeature_RNA`, `percent.mt`), normalisation (`LogNormalize`), variable-feature selection, scaling, PCA, and UMAP
- Runs **SingleR** for reference-based cell-type annotation using `celldex` references (e.g. `HumanPrimaryCellAtlasData`, `BlueprintEncodeData`)
- Generates a cell-type × group table and a UMAP coloured by annotation
- Extracts platelet subpopulation for downstream focused analysis
- Saves the annotated Seurat object as `plt-public.rds`

#### Typical Workflow

1. **Install** — run `SingcellAnalysis-packages_install.R` once.
2. **Prepare data** — place your `.rds` files in `data/`.
3. **Process public data** — run `pub-data-processing.R` to merge, annotate and subset.
4. **Harmonise idents** — run `Dataset Modification.R` to align cell-type labels.
5. **Expression analysis** — run `exp-OUT.R` to compute statistics and significance.
6. **Visualise** — run `Box-bar plotting.R` to produce publication figures.

#### Technical Notes

- `exp-OUT.R` has a clearly marked **user-modifiable section** at the top; the remainder is designed to run without further editing.
- `Dataset Modification.R` assumes a one-to-one cell-type mapping between your own and public datasets. Adjust the `recode` mappings if the ontologies differ.
- SingleR annotation in `pub-data-processing.R` uses multiple references in cascade; change the reference list in the script to suit your tissue or species.

### 3. SCAP-prediction — Genotype-Panel Screening & Prediction ★

Python / R scripts for building the SCAP genotype prediction panel: one-hot encoding of the cohort genotype table, LASSO dimensionality reduction, multi-model consensus feature screening (Random Forest / SVM-RFE / XGBoost + SHAP), final dual-model construction (XGBoost for discrimination, multivariable logistic regression for clinical interpretability and nomogram), held-out test evaluation (ROC / calibration / decision curve analysis / SHAP), and single-sample inference.

The panel is derived from the WGS variant table (T2T-YAO coordinates): each locus is coded 0 = wild type, 1 = heterozygous, 2 = homozygous, −1 = unknown, and expanded into binary indicator columns so penalised models see mutually exclusive dummies rather than an ordinal coding.

Scripts were executed under WSL Ubuntu 22.04 (conda Python 3.12) and Windows 11 (R 4.4.2). Paths in the released scripts have been replaced with placeholders (`/path/to/scap_genotype/...`) — edit the user-configurable section at the top of each script before running.

#### Dependencies

| Tool / Package | Version (used) | Purpose                                          |
|----------------|----------------|--------------------------------------------------|
| Python         | 3.12 (WSL Ubuntu 22.04, conda) | Runtime for the modelling scripts   |
| pandas / numpy | —              | Table & matrix manipulation                      |
| scikit-learn   | —              | Logistic regression, RF, SVM-RFE, scaling, metrics |
| xgboost        | —              | Gradient-boosted trees (screening + final model) |
| shap           | —              | SHAP values (feature importance, beeswarm plots) |
| statsmodels    | —              | Logistic regression with coefficient p-values (nomogram input) |
| matplotlib (+ matplotlib-venn, seaborn) | — | Vector PDF figures, Venn diagram, heatmaps |
| joblib         | —              | Model-asset serialisation                        |
| R              | 4.4.2 (Windows) | LASSO screening runtime                         |
| glmnet / ggplot2 / dplyr | —     | LASSO (L1) logistic regression + figures         |

#### Directory Layout & Scripts

- **`convert_to_onehot.py`** — convert the genotype table (ID, SCAP label, one column per locus coded −1/0/1/2) into one-hot indicator columns (`YAO001_0`, `YAO001_1`, …); genotype values are cast to strings first so pandas treats them as four independent categories, not ordinal numbers.
- **`lasso_selection_glmnet.R`** — LASSO (L1-penalised logistic) feature screening with `glmnet`: 10-fold CV (`type.measure = "auc"`), features extracted at `lambda.1se` (the most parsimonious model within one SE of the optimum — preferred for clinical feature screening over `lambda.min`). Exports three vector PDFs (coefficient path, CV curve, non-zero weights with risk/protective colouring), the weight table, and the reduced dataset (`Final_Data_For_Python_XGBoost.csv`) for the downstream steps.
- **`split_train_test.py`** — stratified 80/20 train/test split preserving the severe/non-severe ratio; the test set never participates in any feature screening or tuning.
- **`consensus_feature_selection.py`** — three screening models run in parallel threads on the training set only: Random Forest (Gini importance), SVM-RFE (linear-kernel recursive elimination, on standardised features), XGBoost + SHAP (mean |SHAP|). Produces full per-model rankings, a SHAP summary figure and a Venn diagram of the Top-K sets, then exports train/test subsets for several consensus sets. `selection_mode`:
  - `"pairwise"` (main pipeline) — three pairwise intersections, the strict 3-model intersection, and the **union of pairwise intersections** (features selected by ≥ 2 models); the final model is trained on `Train/Test_Subset_Union_of_Pairs.csv`;
  - `"topk_union"` (sensitivity analysis) — union of the three Top-K lists (top 20/25/30/35 …), used to probe panel-size sensitivity.
- **`train_final_model.py`** — final panel construction and evaluation: XGBoost (raw 0/1 features) for maximum discrimination + multivariable logistic regression (statsmodels on raw features for nomogram coefficients/p-values; sklearn on standardised features for prediction). Evaluates AUC, Brier score and the Youden-optimal cut-off (sensitivity / specificity / PPV / NPV) on the held-out test set; exports ROC curve, calibration curves, DCA (net benefit computed manually — no external DCA package), the **clinical nomogram** (per-locus points scaled so the largest |coefficient| = 100 points; total points → SCAP probability via the log-odds transform with the fitted intercept) plus a text version of the score system, a SHAP beeswarm plot, and all fitted objects to `Final_SCAP_Model_Assets.pkl`. Optional `remove_collinear = True` drops zero-variance and |r| > 0.95 features before modelling (used in the all-LASSO-features sensitivity analysis to avoid a singular design matrix).
- **`predict_new_samples.py`** — single/batch inference from the saved assets; validates that all panel loci are present and applies the correct preprocessing per model (XGBoost on raw features, logistic regression through the stored scaler).
- **`plot_confusion_corr.py`** — supplementary evaluation from the saved assets: confusion-matrix heatmap on the test set and a feature-correlation heatmap (top-20 features by XGBoost importance) to check residual collinearity.

#### Typical Workflow

1. **Encode** — `convert_to_onehot.py`: genotype table → one-hot feature matrix.
2. **Screen (LASSO)** — `lasso_selection_glmnet.R`: p >> n dimensionality reduction at `lambda.1se` → `Final_Data_For_Python_XGBoost.csv`.
3. **Split** — `split_train_test.py`: stratified 80/20 split.
4. **Screen (consensus)** — `consensus_feature_selection.py` (`selection_mode = "pairwise"`, `top_k = 20`): three-model ranking + Venn + consensus subsets (train aligned, test only aligned).
5. **Model & evaluate** — `train_final_model.py` on `Train/Test_Subset_Union_of_Pairs.csv`: XGBoost + logistic regression, ROC / calibration / DCA / nomogram / SHAP, model assets saved.
6. **Infer** — `predict_new_samples.py` with `Final_SCAP_Model_Assets.pkl` and a new-patient CSV; optionally `plot_confusion_corr.py` for supplementary plots.

#### Technical Notes

- **Why LASSO first** — one-hot expansion makes p far exceed n; unpenalised regression then suffers from multicollinearity (LD between neighbouring loci) and overfitting. The L1 penalty yields a sparse solution, doing variable selection and parameter estimation simultaneously.
- **Test-set isolation** — the held-out test set is never used for feature ranking or hyperparameter choice; it enters only the final evaluation, with its columns aligned to whatever subset the training pipeline selects.
- **Scaling convention (important for inference)** — XGBoost is fitted on the **raw** 0/1 features, the sklearn logistic regression on **standardised** features (the fitted `StandardScaler` is stored in the assets); `predict_new_samples.py` applies each preprocessing path to the right model. Mixing them up silently degrades the logistic predictions.
- **Nomogram maths** — points per locus = coefficient × (100 / max |coefficient|); a total-points target for probability p is `(logit(p) − intercept) × scaling`. The probability axis shares the physical extent of the total-points axis, so the two must stay strictly aligned.
- **DCA without a DCA library** — net benefit is computed directly from the confusion matrix at each threshold (`NB = TP/N − FP/N × t/(1−t)`) with treat-all / treat-none reference curves; no external package needed.
- **Sensitivity analyses covered by parameters** — Top-K union screening (`selection_mode = "topk_union"`) and collinearity wash-out before modelling (`remove_collinear = True`, |r| > 0.95 plus zero-variance removal) reproduce the panel-size and all-LASSO-features variants of the study without separate scripts.
- **Percentage formatting fix** — the working notes printed Youden-cut-off sensitivity/specificity as fractions with a "%" suffix; the released script multiplies by 100 before formatting.
- **statsmodels optimiser** — BFGS with an lbfgs fallback is used for the nomogram logistic fit; if both struggle on an ill-conditioned design matrix, enable `remove_collinear = True` first.

### 4. LPCAT-Machine Learning — DNN Classification ★

Deep neural network for binary classification on tabular clinical data. Includes model training, performance evaluation, and comparative ROC analysis.

#### Dependencies

| Package       | Version (tested) | Purpose                              |
|---------------|------------------|--------------------------------------|
| Python        | ≥ 3.8            | Runtime                              |
| tensorflow    | ≥ 2.6            | Model building / training / TF ops   |
| numpy         | —                | Array operations                     |
| pandas        | —                | CSV I/O, data manipulation           |
| scikit-learn  | —                | Train-test split, ROC/AUC metrics    |
| matplotlib    | —                | ROC curve plotting                   |

```bash
pip install tensorflow numpy pandas scikit-learn matplotlib
```

#### CSV Format

All datasets share the same 27-column layout (zero-indexed):

| Col | Name        | Role  | Description                                  |
|-----|-------------|-------|----------------------------------------------|
| 0   | num         | —     | Patient ID (not used as a feature)           |
| 1   | SCAP        | Label | Binary target: 0 = negative, 1 = positive    |
| 2   | Death       | —     | Mortality outcome                            |
| 3   | Gender      | —     | 0 = male, 1 = female                         |
| 4   | Age         | —     | Age in years                                 |
| 5   | BMI         | —     | Body mass index                              |
| 6   | Temperature | —     | Body temperature (°C)                        |
| 7   | RR          | —     | Respiratory rate                             |
| 8   | HR          | —     | Heart rate                                   |
| 9   | Conscious   | —     | 0 = alert, 1 = altered                       |
| 10  | LPCAT1      | ★     | Lysophosphatidylcholine acyltransferase 1    |
| 11  | WBC         | ★     | White blood cell count                       |
| 12  | NE          | ★     | Neutrophil percentage                        |
| 13  | LY          | ★     | Lymphocyte percentage                        |
| 14  | NLR         | ★     | Neutrophil-to-lymphocyte ratio               |
| 15  | Hb          | —     | Hemoglobin                                   |
| 16  | PLT         | —     | Platelet count                               |
| 17  | BUN         | —     | Blood urea nitrogen                          |
| 18  | Scr         | —     | Serum creatinine                             |
| 19  | ALB         | —     | Albumin                                      |
| 20  | CK          | —     | Creatine kinase                              |
| 21  | ESR         | ★     | Erythrocyte sedimentation rate               |
| 22  | CRP         | ★     | C-reactive protein                           |
| 23  | PCT         | ★     | Procalcitonin                                |
| 24  | Glu         | —     | Glucose                                      |
| 25  | CURB65      | —     | CURB-65 severity score (0–5)                 |
| 26  | PSI         | —     | Pneumonia Severity Index                     |

★ = feature columns used by the model (indices 10, 11, 12, 13, 14, 21, 22, 23).

The label column (`SCAP`, index 1) contains 0 (negative) or 1 (positive). Missing values are represented by `-1`.

#### Dataset Variants

Each variant keeps the same 27-column structure but replaces one feature column with `-1` across all rows. This tests the model's dependence on each feature without changing the input shape.

| File                          | Dropped Feature | Column | Legend Label in `ROCs.py` |
|-------------------------------|-----------------|--------|---------------------------|
| `MachineLearning.csv`         | (none)          | —      | Prediction Model          |
| `MachineLearning-nopct.csv`   | PCT             | 23     | Drop-out PCT              |
| `MachineLearning-nowbc.csv`   | WBC             | 11     | Drop-out WBC              |
| `MachineLearning-noesr.csv`   | ESR             | 21     | Drop-out ESR              |
| `MachineLearning-nocrp.csv`   | CRP             | 22     | Drop-out CRP              |
| `MachineLearning-noly.csv`    | LY              | 13     | Drop-out LY%              |
| `MachineLearning-nolpcat.csv` | LPCAT1          | 10     | Drop-out LPCAT            |
| `MachineLearning-none.csv`    | NE              | 12     | Drop-out NE%              |

#### Scripts

**`Machine_Learning.py` — Train**

Builds a 10-layer DNN and trains on the full dataset.

```bash
python Machine_Learning.py
```

- Reads `data/MachineLearning.csv`
- Splits into train / validation / test (60 / 20 / 20, stratified)
- Trains for 500 epochs with Adam + SparseCategoricalCrossentropy
- Saves model to `models/<timestamp>/`
- Prints test accuracy, loss, and confusion matrix

**`model_performance.py` — Evaluate**

Loads a saved model and evaluates it on a single dataset variant.

```bash
python model_performance.py
```

- Set `MODEL_TIMESTAMP` in the config to match a saved model folder
- Default dataset: `data/MachineLearning-nopct.csv`
- Outputs accuracy, loss, confusion matrix, AUC, and a ROC curve to `results/roc/`

**`ROCs.py` — Compare ROC**

Runs the same saved model against all 8 dataset variants and plots their ROC curves on a single figure.

```bash
python ROCs.py
```

- Expects all 8 CSV files in `data/`
- Generates `results/roc/<timestamp>_ROCs.tiff`

**`Self_Defining_Function.py` — Utilities**

Shared helper functions used by the other scripts. Not meant to be run directly.

| Function            | Used by                | Purpose                              |
|--------------------|------------------------|--------------------------------------|
| `prediction()`      | `Machine_Learning.py`  | Softmax wrapper, returns class labels |
| `prediction_v30()`  | `model_performance.py` | Prediction + metrics + ROC plot      |
| `random_test()`     | (standalone utility)   | Random subset sampling               |
| `auto_log_to_csv()` | (standalone utility)   | Append experiment results to CSV log |

#### Typical Workflow

1. **Train** — run `Machine_Learning.py`, note the model timestamp.
2. **Evaluate** — copy the timestamp into `model_performance.py` config, then run it.
3. **Compare** — copy the timestamp into `ROCs.py`, then run it to see all ROC curves.

#### Technical Notes

- The model uses `from_logits=True`; the final Dense layer has no activation. Prediction helpers in `Self_Defining_Function` prepend a Softmax layer automatically.
- All random seeds are pinned (`RANDOM_STATE = 1`) for reproducibility.
- Edit the `FEATURE_COLUMNS` list in each script if your CSV layout differs.

---

## Notes

- Some raw data files are not included in this repository due to institutional archive requirements, database submission policies, file size limitations, or privacy considerations.
- The uploaded scripts will be cleaned and annotated to improve readability and reproducibility.
- Additional documentation will be added as the repository is updated.
- Server paths, usernames, project identifiers, and sample IDs appearing in the WGS scripts are placeholders for privacy reasons; edit the user-configurable section at the top of each script before running.

## Citation

If you use or refer to this repository, please cite the corresponding graduation thesis or contact the author for further information.

## License

The license for this repository will be determined based on the final scope of the publicly released scripts and any related intellectual property considerations.
