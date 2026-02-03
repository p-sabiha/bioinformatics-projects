# Gene Expression Statistical Analysis
# Course: M.Sc. Medical Bioinformatics - Semester 1
# Skills: R, Statistics, Cell and Molecular Biology

# ============================================
# 1. SETUP AND DATA LOADING
# ============================================

# Set working directory (modify if needed)
# setwd("path/to/project")

# Load data
data <- read.csv("../data/gene_expression.csv", header = TRUE)

# View data structure
cat("=== DATA STRUCTURE ===\n")
str(data)
cat("\n=== FIRST FEW ROWS ===\n")
print(head(data))

# ============================================
# 2. DATA PREPROCESSING
# ============================================

# Extract expression values
normal_samples <- data[, c("normal_1", "normal_2", "normal_3")]
tumor_samples <- data[, c("tumor_1", "tumor_2", "tumor_3")]

# Calculate mean expression for each condition
data$normal_mean <- rowMeans(normal_samples)
data$tumor_mean <- rowMeans(tumor_samples)

# Calculate fold change (tumor vs normal)
data$fold_change <- data$tumor_mean / data$normal_mean
data$log2_fc <- log2(data$fold_change)

# ============================================
# 3. STATISTICAL ANALYSIS
# ============================================

# Perform t-test for each gene
p_values <- c()
for (i in 1:nrow(data)) {
  normal_vals <- as.numeric(normal_samples[i, ])
  tumor_vals <- as.numeric(tumor_samples[i, ])
  test_result <- t.test(normal_vals, tumor_vals)
  p_values <- c(p_values, test_result$p.value)
}

data$p_value <- p_values
data$neg_log10_p <- -log10(data$p_value)

# Determine significance (p < 0.05)
data$significant <- ifelse(data$p_value < 0.05, "Yes", "No")

# ============================================
# 4. RESULTS SUMMARY
# ============================================

cat("\n=== STATISTICAL RESULTS ===\n")
results <- data[, c("gene_id", "gene_name", "normal_mean", "tumor_mean",
                     "log2_fc", "p_value", "significant")]
print(results)

# Significant genes
cat("\n=== SIGNIFICANTLY DIFFERENTIALLY EXPRESSED GENES ===\n")
sig_genes <- data[data$significant == "Yes", c("gene_name", "log2_fc", "p_value")]
sig_genes <- sig_genes[order(sig_genes$p_value), ]
print(sig_genes)

# Upregulated in tumor (log2FC > 0)
cat("\n=== UPREGULATED IN TUMOR ===\n")
upregulated <- sig_genes[sig_genes$log2_fc > 0, ]
print(upregulated)

# Downregulated in tumor (log2FC < 0)
cat("\n=== DOWNREGULATED IN TUMOR ===\n")
downregulated <- sig_genes[sig_genes$log2_fc < 0, ]
print(downregulated)

# ============================================
# 5. VISUALIZATION
# ============================================

# Save plots to output folder
pdf("../output/gene_expression_plots.pdf", width = 10, height = 8)

# Plot 1: Boxplot comparing conditions for top genes
par(mfrow = c(2, 2))

# Select top 4 significant genes for boxplot
top_genes <- head(sig_genes$gene_name, 4)
for (gene in top_genes) {
  idx <- which(data$gene_name == gene)
  normal_vals <- as.numeric(normal_samples[idx, ])
  tumor_vals <- as.numeric(tumor_samples[idx, ])
  boxplot(list(Normal = normal_vals, Tumor = tumor_vals),
          main = paste("Expression of", gene),
          ylab = "Expression Level",
          col = c("lightblue", "salmon"))
}

# Plot 2: Volcano Plot
par(mfrow = c(1, 1))
plot(data$log2_fc, data$neg_log10_p,
     main = "Volcano Plot: Tumor vs Normal",
     xlab = "Log2 Fold Change",
     ylab = "-Log10(p-value)",
     pch = 19,
     col = ifelse(data$significant == "Yes", "red", "gray"))
abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = c(-1, 1), col = "green", lty = 2)
legend("topright",
       legend = c("Significant", "Not Significant"),
       col = c("red", "gray"), pch = 19)

# Plot 3: Bar plot of fold changes
par(mfrow = c(1, 1))
barplot(data$log2_fc,
        names.arg = data$gene_name,
        main = "Log2 Fold Change (Tumor vs Normal)",
        xlab = "Genes",
        ylab = "Log2 Fold Change",
        col = ifelse(data$log2_fc > 0, "salmon", "lightblue"),
        las = 2,
        cex.names = 0.7)
abline(h = 0, col = "black")

dev.off()

cat("\n=== PLOTS SAVED TO ../output/gene_expression_plots.pdf ===\n")

# ============================================
# 6. SAVE RESULTS
# ============================================

write.csv(results, "../output/analysis_results.csv", row.names = FALSE)
cat("=== RESULTS SAVED TO ../output/analysis_results.csv ===\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
