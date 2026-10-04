# =============================================================================
# lasso_selection_glmnet.R
#
# LASSO (L1-penalised logistic regression) feature screening on the one-hot
# encoded genotype table, implemented with glmnet.
#
# Motivation: after one-hot encoding the feature count p far exceeds the
# sample count n (p >> n). Unpenalised logistic regression then suffers from
# multicollinearity (linkage disequilibrium between neighbouring loci) and
# overfitting. The L1 penalty produces a sparse solution, shrinking the
# coefficient of uninformative / redundant features exactly to zero, so
# variable selection and parameter estimation are done simultaneously.
#
# Model : glmnet, family = "binomial", alpha = 1 (LASSO)
# CV    : 10-fold cross-validation, type.measure = "auc"
# Rule  : lambda.1se — the most parsimonious model whose CV accuracy is
#         within 1 standard error of the best (lambda.min). The 1se rule is
#         preferred for feature screening in clinical studies.
#
# Outputs (in output_dir):
#   Fig1_LASSO_Coefficient_Path.pdf   coefficient shrinkage path vs log(lambda)
#   Fig2_LASSO_CV_Error.pdf           10-fold CV curve with lambda.min / 1se
#   Fig3_LASSO_Feature_Weights.pdf    non-zero coefficients (risk vs protective)
#   R_LASSO_Selected_Weights.csv      feature / coefficient table
#   Final_Data_For_Python_XGBoost.csv reduced dataset (ID, SCAP + selected
#                                      features) for downstream modelling
# =============================================================================

library(glmnet)
library(ggplot2)
library(dplyr)

# --- User-configurable section ------------------------------------------------
# One-hot encoded input (output of convert_to_onehot.py)
input_file <- "/path/to/scap_genotype/SCAP_GenoType_OneHot.csv"
output_dir <- "/path/to/scap_genotype/"

# Identifier column of the input table
id_col <- "ID"
# Label column
label_col <- "SCAP"
# Random seed for CV fold assignment
random_seed <- 42
# ------------------------------------------------------------------------------

set.seed(random_seed)

cat("1. Loading one-hot encoded data...\n")
data <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

# Data preparation: everything except ID and SCAP is a feature
X_df <- data[, -which(names(data) %in% c(id_col, label_col))]
X <- as.matrix(X_df)
y <- as.numeric(data[[label_col]])

cat("\n2. Running LASSO regression with 10-fold cross-validation...\n")
fit <- glmnet(X, y, family = "binomial", alpha = 1)
cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 1,
                    type.measure = "auc", nfolds = 10)
lambda_min <- cv_fit$lambda.min  # lambda with the best CV accuracy
lambda_1se <- cv_fit$lambda.1se  # most parsimonious lambda within 1 SE

# Extract non-zero coefficients at lambda.1se
coef_matrix <- as.matrix(coef(cv_fit, s = "lambda.1se"))
coef_df <- data.frame(
  Feature = rownames(coef_matrix),
  Coefficient = as.numeric(coef_matrix[, 1]),
  stringsAsFactors = FALSE
)
final_features <- coef_df %>%
  filter(Coefficient != 0 & Feature != "(Intercept)") %>%
  arrange(Coefficient)

cat(sprintf("\n3. LASSO (lambda.1se) retained %d core genotype features.\n",
            nrow(final_features)))

cat("\n4. Generating vector PDF figures...\n")
# Fig 1: coefficient shrinkage path — how features are eliminated as the
# penalty grows; lines surviving longest carry the strongest signal
pdf(file = paste0(output_dir, "Fig1_LASSO_Coefficient_Path.pdf"),
    width = 7, height = 5.5)
plot(fit, xvar = "lambda", label = TRUE, lwd = 2,
     xlab = "Log Lambda (Penalty parameter)", ylab = "Coefficients")
abline(v = log(c(lambda_min, lambda_1se)), lty = 2, col = "gray")
title("LASSO Coefficient Profile Plot", line = 2.5)
dev.off()

# Fig 2: 10-fold CV curve — the "certificate" against overfitting
pdf(file = paste0(output_dir, "Fig2_LASSO_CV_Error.pdf"),
    width = 7, height = 5.5)
plot(cv_fit)
title("10-Fold Cross-Validation Error Curve", line = 2.5)
dev.off()

# Fig 3: non-zero coefficients (risk factors vs protective factors)
if (nrow(final_features) > 0) {
  final_features$Effect <- ifelse(final_features$Coefficient > 0,
                                  "Risk Factor (>0)", "Protective Factor (<0)")
  final_features$Feature <- factor(final_features$Feature,
                                   levels = final_features$Feature)
  p3 <- ggplot(final_features, aes(x = Coefficient, y = Feature, fill = Effect)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = c("Risk Factor (>0)" = "#d62728",
                                 "Protective Factor (<0)" = "#1f77b4")) +
    theme_bw() +
    labs(title = "Non-zero Coefficients Selected by LASSO",
         x = "Standardized Coefficient (Log-Odds)",
         y = "Genotype Features") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 9),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed")
  # Dynamic height so long feature lists stay readable
  ggsave(filename = paste0(output_dir, "Fig3_LASSO_Feature_Weights.pdf"),
         plot = p3, device = "pdf", width = 8,
         height = max(4, nrow(final_features) * 0.15))
}

cat("\n5. Exporting results for downstream models...\n")
write.csv(final_features,
          paste0(output_dir, "R_LASSO_Selected_Weights.csv"),
          row.names = FALSE)
selected_cols <- c(id_col, label_col, as.character(final_features$Feature))
final_dataset <- data[, selected_cols]
write.csv(final_dataset,
          paste0(output_dir, "Final_Data_For_Python_XGBoost.csv"),
          row.names = FALSE)

cat(sprintf("\nAll done. PDF figures and CSV tables written to: %s\n", output_dir))
