#!/usr/bin/env Rscript
# UpSet plots of significant variant / gene sets across the three reference
# genomes (T2T-CHM13 / GRCh38 / T2T-YAO).
#
# Input format (data/upsetr/UpSetR_example.txt):
#   - Tab-delimited, one column per set. The header names sets as
#     "<PANEL>-<REFERENCE>", e.g. "GCS-CHM13"; panels are separated by an
#     empty column. Each column below the header lists the element IDs
#     (variants or genes) belonging to that set, one per row (ragged columns
#     are allowed).
#   - PANEL naming: <CALLER><ANALYSIS><ANNOTATION>
#       1st letter - variant caller:   G = GATK, D = DeepVariant
#       2nd letter - differential set: C = chi-square test, P = PLINK
#       3rd letter - annotation:       S = SnpEff, V = VEP
#     e.g. GCS = GATK + chi-square + SnpEff, DPV = DeepVariant + PLINK + VEP.
#
# For every panel the three reference-specific sets (CHM13 / GRCh38 / YAO)
# are combined into one UpSet plot and written to <PANEL>.pdf in OUTPUT_DIR.
#
# Edit the user-configurable section before running.

# install.packages("UpSetR")   # uncomment on a fresh machine
library(UpSetR)

# ---------------------------------------------------------------------------
# User-configurable section
# ---------------------------------------------------------------------------

INPUT_FILE <- "data/upsetr/UpSetR_example.txt"
OUTPUT_DIR <- "output/upsetr"

# How many largest intersections to show per plot
N_INTERSECTIONS <- 12

# ---------------------------------------------------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

read_set_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  header <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  body <- strsplit(lines[-1], "\t", fixed = TRUE)
  sets <- list()
  for (j in seq_along(header)) {
    col_name <- trimws(header[j])
    if (!nzchar(col_name)) next                       # panel separator column
    values <- unique(trimws(unlist(
      lapply(body, function(row) if (j <= length(row)) row[j] else NA),
      use.names = FALSE)))
    values <- values[nzchar(values) & !is.na(values)]
    if (length(values)) sets[[col_name]] <- values
  }
  sets
}

panel_of <- function(set_name) sub("-.*$", "", set_name)
ref_of   <- function(set_name) sub("^[^-]*-", "", set_name)

plot_panel <- function(sets, panel, refs = c("CHM13", "GRCh38", "YAO")) {
  chosen <- paste0(panel, "-", refs)
  missing <- chosen[!chosen %in% names(sets)]
  if (length(missing)) {
    message("Panel ", panel, ": skipping, missing set(s): ",
            paste(missing, collapse = ", "))
    return(invisible(NULL))
  }
  set_list <- lapply(chosen, function(n) sets[[n]])
  names(set_list) <- refs
  mat <- fromList(set_list)

  pdf(file.path(OUTPUT_DIR, paste0(panel, ".pdf")), width = 8, height = 6)
  on.exit(dev.off(), add = TRUE)
  print(upset(mat,
              sets = refs,
              order.by = "freq",
              nintersects = N_INTERSECTIONS,
              keep.order = TRUE,
              point.size = 3,
              line.size = 1.2,
              main.bar.color = "#377EB8",
              sets.bar.color = "#E41A1C",
              text.scale = c(1.6, 1.4, 1.2, 1.0, 1.6, 1.2)))
  # Panel label drawn with grid (NOT via upset(plot.title=...): the internal
  # ggplot_build call of UpSetR 1.4.0 is incompatible with ggplot2 >= 3.5)
  grid::grid.text(panel, x = 0.5, y = 0.99, just = "top",
                  gp = grid::gpar(fontsize = 14, fontface = "bold"))
  message("Saved: ", file.path(OUTPUT_DIR, paste0(panel, ".pdf")))
}

sets <- read_set_file(INPUT_FILE)
message("Loaded ", length(sets), " set(s): ", paste(names(sets), collapse = ", "))

for (panel in unique(panel_of(names(sets)))) {
  plot_panel(sets, panel)
}
