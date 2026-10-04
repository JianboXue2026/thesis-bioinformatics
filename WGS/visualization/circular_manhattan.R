# circular_manhattan.R
# Circular Manhattan plot with the CMplot R package.
#
# NOTE: CMplot handles at most ~1,000,000 variants — downsample / filter the
# association results (e.g. random sampling + OR normalisation) before
# plotting when the full result set is larger.

# install.packages("CMplot")   # uncomment on a fresh machine
library(CMplot)

# --- User-configurable section ------------------------------------------------
# Input CSV: one row per variant with columns
#   SNP, Chromosome, Position, P-value
# (CMplot treats every extra numeric column as a P-value track — keep the
# table to these four columns, or select the P-value column before plotting.)
# The default points at a small SYNTHETIC example dataset
# (data/circular_manhattan_example.csv) so the script runs out of the box;
# replace it with your own association results. The complete per-variant
# comparison table is archived at Peking University and is not distributed.
input_csv <- "data/circular_manhattan_example.csv"
# -----------------------------------------------------------------------------

data_01 <- read.csv(input_csv)

# Density binning: CMplot colours 1 Mb windows by marker count, and the
# break boundaries must not exceed the densest window in the data. The
# breaks are therefore scaled to the input (with million-variant results
# this reaches roughly the 100-500 range).
bin_size <- 1e6
markers_per_bin <- table(trunc(data_01[["Position"]] / bin_size))
bin_breaks <- seq(1, max(markers_per_bin), length.out = 9)

CMplot(data_01,
       type = "p",                    # points
       plot.type = "c",               # circular Manhattan plot
       chr.labels = paste("Chr", c(1:22), sep = ""),
       chr.pos.max = TRUE,
       r = 0.4,
       bin.size = bin_size,
       bin.breaks = bin_breaks,
       outward = FALSE,
       cir.chr.h = 1.3,
       chr.den.col = c("darkgreen", "yellow", "red"),
       file = "tiff",
       dpi = 600,
       threshold = c(NULL, 1e-5),
       threshold.col = c("white", "red"),
       file.output = TRUE, verbose = TRUE,
       width = 10, height = 10)
