# manhattan_qq.R
# Manhattan and QQ plots for PLINK association results (run locally in R).
#
# Input: a PLINK .assoc file. The script cleans chromosome codes (removes
# PAR1/PAR2 rows, maps X/Y/MT to 23/24/25) and drops invalid P values before
# plotting.

# Load the qqman package
library(qqman)

# --- User-configurable section ------------------------------------------------
assoc_file <- "merged-YAO.assoc"   # PLINK association result
out_dir <- "."                     # output directory for the PNG files
# -----------------------------------------------------------------------------

# Read GWAS results
gwas_data <- read.table(assoc_file, header = TRUE)

# Inspect chromosome codes
print(unique(gwas_data$CHR))

# Remove PAR1 / PAR2 rows
gwas_data <- gwas_data[!(gwas_data$CHR %in% c("PAR1", "PAR2")), ]

# Map non-numeric chromosome codes to numbers
gwas_data$CHR <- as.character(gwas_data$CHR)
gwas_data$CHR[gwas_data$CHR == "X"]  <- "23"
gwas_data$CHR[gwas_data$CHR == "Y"]  <- "24"
gwas_data$CHR[gwas_data$CHR == "MT"] <- "25"
print(unique(gwas_data$CHR))

gwas_data$CHR <- as.numeric(gwas_data$CHR)
print(unique(gwas_data$CHR))

# Drop invalid P values (NA / Inf / <= 0)
gwas_data <- gwas_data[!is.na(gwas_data$P) & !is.infinite(gwas_data$P) & gwas_data$P > 0, ]

# Manhattan plot
png(file.path(out_dir, "manhattan_plot.png"))
manhattan(gwas_data, chr = "CHR", bp = "BP", snp = "SNP", p = "P",
          main = "Manhattan Plot")
dev.off()

# QQ plot
png(file.path(out_dir, "qq_plot.png"))
qq(gwas_data$P)
dev.off()
