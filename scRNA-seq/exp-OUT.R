# ============================== [User: modify only this section] ==============================
# 1. Specify the Seurat object to analyze (change only this line when switching objects)
seurat_obj <- filtered_cells_LY6E_self          # <-- Replace with your object name, e.g., PKU.combined.COVID, obj1, obj2, etc.
obj_name <- "self-LY6E"
data_dir <- "path/to/your/project/20251204/"
# 2. Genes of interest (modify as needed)
#genes_of_interest <- c("IL15RA", "CARD19", "IFNAR2", "LY6E", "FHL1")
genes_of_interest <- "LY6E"

# 3. Group column name (column in meta.data)
group_col <- "compub"                # <-- Your group column name (2+ groups supported)

# 4. (Optional) Output file prefix to distinguish results from different objects
file_prefix <- "result"              # Final filename: result_average_expression_YYYYMMDD.csv, etc.
# ============================================================================

# ---------------------- Code below: no modification needed, run directly ----------------------

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)

# Check whether the Seurat object exists
if (!exists(deparse(substitute(seurat_obj)))) {
  stop("Error: Seurat object '", deparse(substitute(seurat_obj)), "' does not exist! Please check the variable name.")
}

cat("=== Analyzing object:", obj_name, "===\n")
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of genes:", nrow(seurat_obj), "\n")
cat("Available assays:", paste(names(seurat_obj@assays), collapse = ", "), "\n")
cat("meta.data columns:", paste(colnames(seurat_obj@meta.data), collapse = ", "), "\n")

cat("Current cell types (Idents):", paste(levels(Idents(seurat_obj)), collapse = ", "), "\n")
cat("Cell type summary:\n")
print(table(Idents(seurat_obj)))

cat("First 10 gene names:", paste(head(rownames(seurat_obj), 10), collapse = ", "), "\n")

# Check whether requested genes exist in the object
available_genes <- rownames(seurat_obj)
missing_genes <- setdiff(genes_of_interest, available_genes)
if (length(missing_genes) > 0) {
  cat("Warning: the following genes were not found:", paste(missing_genes, collapse = ", "), "\n")
}
existing_genes <- intersect(genes_of_interest, available_genes)
if (length(existing_genes) == 0) {
  stop("Error: none of the specified genes were found! Please check gene names.")
}
cat("Genes to analyze:", paste(existing_genes, collapse = ", "), "\n")

# Check group column
if (!group_col %in% colnames(seurat_obj@meta.data)) {
  stop(paste("Error: group column", group_col, "does not exist in meta.data!"))
}
cat("Group summary:\n")
print(table(seurat_obj@meta.data[[group_col]]))

# ====================== Compute expression statistics (including quantiles) ======================
results_avg_expression <- data.frame()
results_expression_proportion <- data.frame()

for (gene in existing_genes) {
  cat("Processing gene:", gene, "\n")

  gene_data <- FetchData(seurat_obj, vars = c(gene, group_col))
  gene_data$cell_type <- as.character(Idents(seurat_obj))
  colnames(gene_data)[1:2] <- c("expression", "group")
  gene_data <- gene_data[complete.cases(gene_data), ]
  if (nrow(gene_data) == 0) next

  # Mean expression + median + quantiles
  avg_expression <- gene_data %>%
    group_by(cell_type, group) %>%
    summarise(
      mean_expression   = mean(expression, na.rm = TRUE),
      median_expression = median(expression, na.rm = TRUE),
      sd_expression     = sd(expression, na.rm = TRUE),
      q05 = quantile(expression, 0.05, na.rm = TRUE),
      q25 = quantile(expression, 0.25, na.rm = TRUE),
      q75 = quantile(expression, 0.75, na.rm = TRUE),
      q95 = quantile(expression, 0.95, na.rm = TRUE),
      n_cells = n(),
      .groups = 'drop'
    ) %>%
    mutate(gene = gene)

  # Expression proportion
  expression_proportion <- gene_data %>%
    group_by(cell_type, group) %>%
    summarise(
      total_cells = n(),
      expressing_cells = sum(expression > 0, na.rm = TRUE),
      proportion = expressing_cells / total_cells,
      .groups = 'drop'
    ) %>%
    mutate(gene = gene)

  results_avg_expression <- rbind(results_avg_expression, avg_expression)
  results_expression_proportion <- rbind(results_expression_proportion, expression_proportion)
}

# ====================== Multi-group pairwise significance testing ======================
perform_significance_tests_multi <- function(obj, genes, gcol) {
  sig_res <- data.frame()

  for (gene in genes) {
    cat("Significance test ->", gene, "\n")

    dat <- FetchData(obj, vars = c(gene, gcol))
    dat$cell_type <- as.character(Idents(obj))
    colnames(dat)[1:2] <- c("expression", "group")
    dat <- dat[complete.cases(dat), ]
    if (nrow(dat) == 0) next

    cell_types <- unique(dat$cell_type)
    groups <- unique(dat$group)
    if (length(groups) < 2) next

    group_pairs <- combn(groups, 2, simplify = FALSE)

    for (ct in cell_types) {
      ct_dat <- dat[dat$cell_type == ct, ]
      if (nrow(ct_dat) < 10) next

      for (pair in group_pairs) {
        g1 <- pair[1]; g2 <- pair[2]
        d1 <- ct_dat[ct_dat$group == g1, "expression"]
        d2 <- ct_dat[ct_dat$group == g2, "expression"]
        if (length(d1) < 5 || length(d2) < 5) next

        p_wilcox <- tryCatch(wilcox.test(d1, d2, exact = FALSE)$p.value, error = function(e) NA)

        tab <- matrix(c(sum(d1 > 0), sum(d1 <= 0), sum(d2 > 0), sum(d2 <= 0)), nrow = 2, byrow = TRUE)
        p_fisher <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA)

        sig_res <- rbind(sig_res, data.frame(
          gene = gene,
          cell_type = ct,
          group1 = g1,
          group2 = g2,
          n_group1 = length(d1),
          n_group2 = length(d2),
          expressing_group1 = sum(d1 > 0),
          expressing_group2 = sum(d2 > 0),
          wilcox_pval = p_wilcox,
          fisher_pval = p_fisher,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  if (nrow(sig_res) > 0) {
    sig_res$wilcox_pval_adj <- p.adjust(sig_res$wilcox_pval, method = "BH")
    sig_res$fisher_pval_adj <- p.adjust(sig_res$fisher_pval, method = "BH")
  }

  return(sig_res)
}

# Run significance tests
significance_results <- data.frame()
if (length(existing_genes) > 0) {
  significance_results <- perform_significance_tests_multi(seurat_obj, existing_genes, group_col)
}

# ====================== Output & save ======================
date_tag <- format(Sys.Date(), "%Y%m%d")

cat("\n=== Average expression results (with quantiles) ===\n")
print(head(results_avg_expression))

cat("\n=== Expression proportion results ===\n")
print(head(results_expression_proportion))

cat("\n=== Significance test results (pairwise) ===\n")
if (nrow(significance_results) > 0) {
  print(head(significance_results))
} else {
  cat("No significant results\n")
}

# Save files (prefix includes object name + date)
prefix <- paste0(file_prefix, "_", obj_name, "_")

write.csv(results_avg_expression,
          paste0(data_dir, prefix, "average_expression_", date_tag, ".csv"),
          row.names = FALSE)
cat("Average expression saved ->", paste0(data_dir, prefix, "average_expression_", date_tag, ".csv"), "\n")

write.csv(results_expression_proportion,
          paste0(data_dir, prefix, "expression_proportion_", date_tag, ".csv"),
          row.names = FALSE)
cat("Expression proportion saved ->", paste0(data_dir, prefix, "expression_proportion_", date_tag, ".csv"), "\n")

if (nrow(significance_results) > 0) {
  write.csv(significance_results,
            paste0(data_dir, prefix, "significance_test_", date_tag, ".csv"),
            row.names = FALSE)
  cat("Significance test saved ->", paste0(data_dir, prefix, "significance_test_", date_tag, ".csv"), "\n")
}

cat("\n", obj_name, "analysis complete!\n")
