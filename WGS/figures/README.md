# WGS Figures — Downstream Thesis Figure Scripts

This directory holds the scripts that produced the **downstream WGS figures of
the thesis**: cohort mapping-QC comparisons across the three reference genomes
(T2T-CHM13 / GRCh38 / T2T-YAO), per-sample difference plots, per-chromosome
variant-density ideograms, and UpSet comparisons of variant / gene sets.

The analysis pipeline itself (FASTQ → BAM → VCF → association) lives in the
sibling directories (`preprocessing/`, `mapping/`, `variant-calling/`, …); see
the repository README for the full workflow.

All scripts are parameterised, with a user-configurable section at the top —
edit it before running with your own data.

## Directory Layout

```
figures/
├── samtools_stats_extract.py       # samtools stats -> per-sample metric CSVs
├── box_line_plots.py               # box-line plots per metric (3 references)
├── difference_plots.py             # per-sample difference scatter plots
├── make_density_windows.py         # 1 Mb window position tables for density tracks
├── chromosome_density_ideogram.R   # RIdeogram karyotype + density + markers
├── upsetr_plot.R                   # UpSet plots of variant/gene sets
├── output/                         # generated figures (auto-created)
└── data/                           # example datasets (anonymised IDs)
    ├── boxline/                    # wide + long metric tables
    ├── difference/                 # YAO-vs-others per-sample differences
    ├── chromosome-density/         # per-1Mb density + gene marker tables
    └── upsetr/                     # tab-delimited multi-set list
```

**The data files use anonymised placeholder sample IDs** (20 samples,
`CAP001…CAP020`) and reproduce the exact file layouts so every script runs
out of the box. The per-sample QC tables (`boxline/`, `difference/`) contain
the **same aggregate sample-level values as reported in the thesis** —
disclosed here deliberately, since per-sample QC aggregates carry no
individual-level genomic or clinical data and nothing that identifies
patients. The chromosome-density, gene-marker and UpSetR tables ship with
placeholder values/IDs to document the layout. The complete clinical data are
archived at Peking University and are not publicly distributed. Replace the
files in `data/` with your own tables (identical column layout) to reproduce the
thesis figures.

> **Two separate cohorts.** The WGS cohort shown here and the clinical cohort in
> `../LPCAT-Machine-Learning/` are entirely different study cohorts with **no
> overlapping patients**; identical or similarly numbered identifiers across the
> two datasets do **not** refer to the same individual.

## Scripts

### `samtools_stats_extract.py` — Metric Extraction from samtools stats

Walks a directory of `<sample_id>_<reference>.stats` reports and extracts one
selected `SN` metric (e.g. `reads mapped and paired`, `error rate`) into a
long-format CSV, plus an optional wide-format CSV (one column per reference)
that feeds `box_line_plots.py`.

```bash
python samtools_stats_extract.py \
    --stats-dir /path/to/01.Samtools_stats \
    --metric mapped_and_paired \
    --output data/boxline/raw_mapped_and_paired.csv \
    --wide-output data/boxline/box-line_mapped_and_paired.csv
python samtools_stats_extract.py --list-metrics   # show available metrics
```

- Expects the file-naming convention produced by `../mapping/samtools_index_stats.sh`.

### `box_line_plots.py` — Box–Line Comparison Plots

For each alignment/QC metric, draws a box plot per reference genome, overlays
one line per sample (the comparison is paired — every sample was aligned
against all three references), annotates Q25 / median / Q75, and marks
pairwise Wilcoxon signed-rank tests with significance stars (`*` < 0.05,
`**` < 0.01, `***` < 0.001).

```bash
python box_line_plots.py                       # all eight metrics
python box_line_plots.py --metric error_rate   # a single metric
```

- Input: `data/boxline/box-line_<metric>.csv` with columns
  `Sample_ID, CHM13, GRCh38, YAO`.
- Samples whose value for any reference falls outside the group-wise 5%–95%
  quantile range are dropped from the line overlay (the box plot keeps the
  full data) — the whole sample drops, never a single reference value.
- Output: `output/boxline/box-line_<metric>.pdf`.
- Covered metrics: `bases_mapped_cigar`, `mismatches`, `error_rate`,
  `mapped_and_paired`, `reads_unmapped`, `Het_Mutant`, `Hom_Mutant`,
  `All_Mutant`.

### `difference_plots.py` — Per-Sample Difference Scatter Plots

Plots, per sample, how much smaller/larger each metric is on T2T-YAO than on
T2T-CHM13 and GRCh38 (YAO as baseline), sorted by the GRCh38 difference.

```bash
python difference_plots.py                          # all three metrics
python difference_plots.py --metric mapped_bases
```

- Input: `data/difference/<metric>.csv` with columns
  `Sample_ID, YAO-CHM13 (Mb), YAO-GRCh38 (Mb)`.
- Two variants per metric, as in the thesis: the original PNG style (both
  differences plus the YAO baseline at 0, inverted y-axis) and the renewed
  PDF style (reductions relative to YAO only).
- Output: `output/difference/difference_<metric>-{png,renew}.{png,pdf}`.

### `make_density_windows.py` — 1 Mb Window Position Tables

Generates the per-chromosome 1 Mb window position table (`Chr, Start, End,
Count`) used to build the variant-density tracks. Chromosome lengths come
from a reference FASTA (`.fai` is created via `samtools faidx` if absent) or
a two-column `Chr,End` CSV.

```bash
python make_density_windows.py --fasta reference.fasta --output windows.csv
```

Fill the `Count` column with per-window variant counts (counted upstream from
the filtered cohort VCF), then normalise to 0–1 into the `Value` column
consumed by `chromosome_density_ideogram.R`.

### `chromosome_density_ideogram.R` — Karyotype + Variant Density + Gene Markers

Draws a full 24-chromosome ideogram for each reference genome with RIdeogram:
the per-1Mb variant-density heat map is overlaid on the chromosomes and
candidate (high-impact) genes are marked on top. Each SVG is converted to
PDF (300 dpi) and TIFF (900 dpi).

```r
Rscript chromosome_density_ideogram.R
```

- Karyotype chromosome lengths are hardcoded in the script — these are
  public assembly properties, not study data.
- Density input: `data/chromosome-density/Chrome_<REF>_Variants_Density.csv`
  with columns `Chr, Start, End, Value` (Value normalised to 0–1, one row per
  1 Mb window).
- Gene marker input: `data/chromosome-density/ChrPlot-<REF>-HighImpact-Gene.csv`
  with columns `Type, Shape, Chr, Start, End, color` (`color` as an
  `"R,G,B"` string) — the RIdeogram `marker` label format.

### `upsetr_plot.R` — UpSet Plots of Variant / Gene Sets

Combines the three reference-specific element lists (CHM13 / GRCh38 / YAO)
of each set category into one UpSet plot per category.

```r
Rscript upsetr_plot.R
```

- Input: `data/upsetr/UpSetR_example.txt` — tab-delimited, one column per
  set named `<PANEL>-<REFERENCE>`, panels separated by an empty column, each
  column listing element IDs (one per row, ragged columns allowed).
- **Panel naming** — `<CALLER><ANALYSIS><ANNOTATION>`:
  - 1st letter — variant caller: `G` = GATK, `D` = DeepVariant;
  - 2nd letter — differential set: `C` = chi-square test, `P` = PLINK;
  - 3rd letter — annotation: `S` = SnpEff, `V` = VEP.
  - e.g. `GCS` = GATK + chi-square + SnpEff, `DPV` = DeepVariant + PLINK + VEP.
- Output: `output/upsetr/<PANEL>.pdf`.
- Note: UpSetR 1.4.0 is incompatible with ggplot2 ≥ 3.5 through its internal
  `ggplot_build` call — the panel label is therefore drawn with `grid.text`
  instead of `upset(plot.title = ...)`.

## Typical Workflow

1. **Extract metrics** — after `samtools_index_stats.sh`, run
   `samtools_stats_extract.py` once per metric to build the long + wide CSVs.
2. **Compare references** — `box_line_plots.py` for the per-metric box-line
   figures; `difference_plots.py` for the YAO-vs-others scatter figures.
3. **Density ideograms** — `make_density_windows.py` to build the 1 Mb window
   table, count variants per window from the filtered cohort VCFs upstream,
   normalise values to 0–1, then run `chromosome_density_ideogram.R`.
4. **Set comparisons** — export the significant-variant / gene ID lists per
   reference into the tab-delimited format, then run `upsetr_plot.R`.

## Privacy Notes

- All data files under `data/` use **anonymised placeholder sample IDs** and
  contain **aggregate, sample-level statistics only** — no individual-level
  genomic or clinical data, and nothing that identifies patients. The QC
  metric tables (`boxline/`, `difference/`) disclose the **same aggregate
  per-sample values as reported in the thesis** (a deliberate choice: these
  aggregates cannot be linked back to individuals). The chromosome-density,
  gene-marker and UpSetR tables use placeholder values/IDs. The complete
  clinical data are archived at Peking University and are not publicly
  distributed.
- The WGS cohort in this directory and the clinical cohort in
  `../LPCAT-Machine-Learning/` are **two separate cohorts with no overlapping
  patients**; matching or similarly numbered IDs across the two are unrelated
  individuals.
