#!/usr/bin/env Rscript
# Chromosome ideograms with variant-density overlay and candidate-gene markers
#
# For each of the three reference genomes (T2T-CHM13, GRCh38, T2T-YAO) this
# script draws a full karyotype with RIdeogram:
#   - the per-1Mb variant-density track is overlaid on the chromosomes
#     ("overlaid" heat-map layer), and
#   - candidate genes (high-impact / significant loci) are marked on top
#     ("marker" label layer).
# Output: one SVG per reference, converted to PDF (300 dpi) and TIFF (900 dpi).
#
# Inputs (example datasets shipped in data/chromosome-density/):
#   - karyotype: chromosome lengths are hardcoded below; they are public
#     reference-genome properties (assembly sizes), NOT study data.
#   - density CSV: Chr, Start, End, Value  (Value normalised to 0-1)
#   - gene marker CSV: Type, Shape, Chr, Start, End, color
#     (Type = gene label, Shape = marker shape, color = "R,G,B" string)
#
# Edit the user-configurable section before running with your own data.

# install.packages("RIdeogram")   # uncomment on a fresh machine
library(RIdeogram)

# ---------------------------------------------------------------------------
# User-configurable section
# ---------------------------------------------------------------------------

DATA_DIR <- "data/chromosome-density"
OUTPUT_DIR <- "output/chromosome-density"

# ---------------------------------------------------------------------------
# Karyotypes: chromosome lengths (bp) of the three reference assemblies.
# Public information - safe to ship.
# ---------------------------------------------------------------------------

CHM13_chrom <- data.frame(
  Chr = c(1:22, "X", "Y"),
  Start = rep(0, 24),
  End = c(248387328, 242696752, 201105948, 193574945, 182045439, 172126628,
          160567428, 146259331, 150617247, 134758134, 135127769, 133324548,
          113566686, 101161492, 99753195, 96330374, 84276897, 80542538,
          61707364, 66210255, 45090682, 51324926, 154259566, 62460029)
)

GRCh38_chrom <- data.frame(
  Chr = c(1:22, "X", "Y"),
  Start = rep(0, 24),
  End = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
          159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
          114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
          58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
)

YAO_chrom <- data.frame(
  Chr = c(1:22, "X", "Y"),
  Start = rep(0, 24),
  End = c(243724521, 242614988, 200304709, 191838926, 183358055, 171543956,
          160820566, 145070391, 132347575, 135075359, 134119163, 134464661,
          106548181, 104236899, 97185997, 88312866, 83649057, 80471738,
          60544037, 66996478, 42396787, 50500499, 155058850, 51523716)
)

# Reference -> (karyotype, density CSV, gene marker CSV)
REFS <- list(
  CHM13 = list(karyotype = CHM13_chrom,
               density = "Chrome_CHM13_Variants_Density.csv",
               genes = "ChrPlot-CHM13-HighImpact-Gene.csv"),
  GRCh38 = list(karyotype = GRCh38_chrom,
                density = "Chrome_GRCh38_Variants_Density.csv",
                genes = "ChrPlot-GRCh38-HighImpact-Gene.csv"),
  YAO = list(karyotype = YAO_chrom,
             density = "Chrome_YAO_Variants_Density.csv",
             genes = "ChrPlot-YAO-HighImpact-Gene.csv")
)

# ---------------------------------------------------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (ref in names(REFS)) {
  cfg <- REFS[[ref]]

  gene_markers <- read.csv(file.path(DATA_DIR, cfg$genes), sep = ",",
                           stringsAsFactors = FALSE)
  gene_density <- read.table(file.path(DATA_DIR, cfg$density), sep = ",",
                             header = TRUE, stringsAsFactors = FALSE)

  svg_file <- file.path(OUTPUT_DIR, paste0("Chrome_Variants_Density-", ref, ".svg"))
  ideogram(karyotype = cfg$karyotype,
           overlaid = gene_density,
           label = gene_markers,
           label_type = "marker",
           output = svg_file)

  base_name <- file.path(OUTPUT_DIR,
                         paste0("Chrome_Variants_Density&Marker-", ref))
  svg2pdf(svg = svg_file, file = base_name,
          width = 8.2677, height = 11.6929, dpi = 300)
  svg2tiff(svg = svg_file, file = base_name,
           width = 8.2677, height = 11.6929, dpi = 900)

  message("Finished reference: ", ref)
}
