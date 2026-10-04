# circular_manhattan.R
# Circular Manhattan plot with the CMplot R package.
#
# NOTE: CMplot handles at most ~1,000,000 variants — downsample / filter the
# association results (e.g. random sampling + OR normalisation) before
# plotting when the full result set is larger.

install.packages("CMplot")   # once
library(CMplot)

# --- User-configurable section ------------------------------------------------
# Input CSV: one row per variant with chromosome, position and P value
input_csv <- "merged-YAO-nofil-2csv-fillup-adjusted-norm-selected.csv"
# -----------------------------------------------------------------------------

data_01 <- read.csv(input_csv)

CMplot(data_01,
       type = "p",                    # points
       plot.type = "c",               # circular Manhattan plot
       chr.labels = paste("Chr", c(1:22), sep = ""),
       chr.pos.max = TRUE,
       r = 0.4,
       bin.size = 1e6,
       bin.breaks = seq(100, 500, 50),
       outward = FALSE,
       cir.chr.h = 1.3,
       chr.den.col = c("darkgreen", "yellow", "red"),
       file = "tiff",
       dpi = 600,
       threshold = c(NULL, 1e-5),
       threshold.col = c("white", "red"),
       file.output = TRUE, verbose = TRUE,
       width = 10, height = 10)
