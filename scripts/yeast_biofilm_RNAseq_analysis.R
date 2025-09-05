# Yeast Biofilm RNA-Seq Analysis 


# Load required packages
library(edgeR)
library(ggplot2)
library(pheatmap)
library(org.Sc.sgd.db)
library(dplyr)
library(tidyr)

# 1. Data Preparation ----
count_data <- read.delim("./gene_counts_final.txt", 
                         header=TRUE, row.names=1, skip=1)
yeast_counts <- as.matrix(count_data[, 6:ncol(count_data)])  # Extract count columns
colnames(yeast_counts) <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(yeast_counts))

# Create metadata
sample_metadata <- data.frame(
  SampleID = colnames(yeast_counts),
  BiofilmStage = factor(rep(c("Early", "Thin", "Mature"), each = 3),
                        levels = c("Early", "Thin", "Mature")),
  row.names = colnames(yeast_counts)
)

# 2. Differential Expression Analysis ----
yeast_dge <- DGEList(counts = yeast_counts, 
                     samples = sample_metadata,
                     group = sample_metadata$BiofilmStage)

# Filter and normalize
keep <- filterByExpr(yeast_dge)
yeast_dge <- yeast_dge[keep, , keep.lib.sizes = FALSE]
yeast_dge <- calcNormFactors(yeast_dge, method = "TMM")

# Design matrix
design <- model.matrix(~0 + BiofilmStage, data = yeast_dge$samples)
colnames(design) <- levels(yeast_dge$samples$BiofilmStage)

# Dispersion and GLM
yeast_dge <- estimateDisp(yeast_dge, design, robust = TRUE)
fit <- glmQLFit(yeast_dge, design)

# Contrasts 
biofilm_contrasts <- list(
  Early_vs_Thin = c(-1, 1, 0),
  Early_vs_Mature = c(-1, 0, 1),
  Thin_vs_Mature = c(0, -1, 1),  
  Early_vs_LateStages = c(2, -1, -1)/1
)

# Run tests
deg_results <- lapply(biofilm_contrasts, function(con) {
  test <- glmQLFTest(fit, contrast = con)
  res <- topTags(test, n = Inf)$table
  res$GeneSymbol <- mapIds(org.Sc.sgd.db, 
                           keys = rownames(res),
                           column = "GENENAME",
                           keytype = "ORF",
                           multiVals = "first")
  return(res)
})

# Pre-calculate CPM for all visualizations
log_cpm <- cpm(yeast_dge, log = TRUE, prior.count = 2)

# 3. Create Output Directory ----
output_dir <- "Biofilm_RNAseq_Results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 4. Enhanced Volcano Plot ----
volcano_data <- deg_results$Early_vs_LateStages %>%
  mutate(Significance = case_when(
    FDR < 0.05 & logFC > 1 ~ "Up in Early",
    FDR < 0.05 & logFC < -1 ~ "Up in Late",
    TRUE ~ "NS"
  ))

ggplot(volcano_data, aes(x = logFC, y = -log10(FDR), color = Significance)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c("Up in Early" = "#E64B35", 
                                "Up in Late" = "#3182BD", 
                                "NS" = "grey80")) +
  labs(title = "Yeast Biofilm Developmental Transitions",
       subtitle = "Early vs Combined Thin+Mature Biofilm Stages",
       x = "log2 Fold Change (Early/Late)",
       y = "-log10 Adjusted p-value",
       color = "Expression Trend") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

ggsave(file.path(output_dir, "Biofilm_Volcano.png"), width = 8, height = 6, dpi = 300)

# 5. Heatmap of Top 30 DEGs (Early vs Late) ----
top_genes <- deg_results$Early_vs_LateStages %>%
  arrange(FDR) %>%
  head(30) %>%
  rownames()

heatmap_annot <- data.frame(
  BiofilmStage = sample_metadata$BiofilmStage,
  row.names = colnames(log_cpm)
)

pheatmap(log_cpm[top_genes, ],
         annotation_col = heatmap_annot,
         show_rownames = TRUE,
         scale = "row",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top 30 DEGs: Early vs Late Biofilm",
         filename = file.path(output_dir, "Biofilm_Heatmap.png"),
         width = 8, height = 6)



# 6.MA Plots ----
fast_maplot <- function(deg_data, contrast_name) {
  p <- deg_data %>%
    ggplot(aes(x = logCPM, y = logFC, 
               color = FDR < 0.05 & abs(logFC) > 1)) +
    geom_point(alpha = 0.3, size = 0.8) +
    scale_color_manual(values = c("grey70", "red")) +
    labs(title = paste("MA Plot:", contrast_name)) +
    theme_minimal()
  
  ggsave(file.path(output_dir, paste0("MA_", gsub(" ", "_", contrast_name), ".png")),
         plot = p, width = 6, height = 5, dpi = 150)
}

# Generate all MA plots
fast_maplot(deg_results$Early_vs_Thin, "Early vs Thin")
fast_maplot(deg_results$Early_vs_Mature, "Early vs Mature")
fast_maplot(deg_results$Thin_vs_Mature, "Thin vs Mature")

# 7. Save Required Results ----
# Objective 1: Pairwise comparisons
write.csv(deg_results$Early_vs_Thin, 
          file.path(output_dir, "DEGs_Early_vs_Thin.csv"),
          row.names = TRUE)
write.csv(deg_results$Early_vs_Mature, 
          file.path(output_dir, "DEGs_Early_vs_Mature.csv"),
          row.names = TRUE)
write.csv(deg_results$Thin_vs_Mature, 
          file.path(output_dir, "DEGs_Thin_vs_Mature.csv"),
          row.names = TRUE)

# Objective 2: Early vs Combined Late
write.csv(deg_results$Early_vs_LateStages, 
          file.path(output_dir, "Biofilm_DEG_Results.csv"),
          row.names = TRUE)

# 8. Generate Report ----
cat(paste0(
  "Yeast Biofilm RNA-Seq Analysis Complete\n",
  "======================================\n",
  "Key findings:\n",
  "- ", sum(deg_results$Early_vs_Thin$FDR < 0.05), " DEGs (Early vs Thin)\n",
  "- ", sum(deg_results$Early_vs_Mature$FDR < 0.05), " DEGs (Early vs Mature)\n",
  "- ", sum(deg_results$Thin_vs_Mature$FDR < 0.05), " DEGs (Thin vs Mature)\n",
  "- ", sum(volcano_data$Significance == "Up in Early"), " genes upregulated in early biofilm\n",
  "- ", sum(volcano_data$Significance == "Up in Late"), " genes upregulated in late stages\n",
  "- Top DEG: ", rownames(deg_results$Early_vs_LateStages)[1], 
  " (", deg_results$Early_vs_LateStages$GeneSymbol[1], ")\n",
  "\nOutput files saved to: ", normalizePath(output_dir)
))