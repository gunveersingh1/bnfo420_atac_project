# =========================================
# ATAC Peak Enrichment Analysis (FINAL CLEAN)
# =========================================

library(GenomicRanges)
library(ChIPseeker)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)
library(tibble)
library(stringr)
library(tidyr)
library(AnnotationDbi)
library(pheatmap)

setwd("/Users/gunveersingh/Documents/VCU_BNFO/atac_project")

# -----------------------------
# Load ATAC peaks
# -----------------------------
peaks_df <- read.table(
  "GSE296875_mergedConsensusPeaks_peaksCalledOnMergedBams_autosomalOnly_withCellTypeOverlap.bed",
  sep = "\t",
  header = FALSE
)

colnames(peaks_df)[1:3] <- c("chr", "start", "end")

peak_gr <- GRanges(
  seqnames = peaks_df$chr,
  ranges = IRanges(start = peaks_df$start,
                   end = peaks_df$end)
)

# -----------------------------
# Annotate peaks to genes
# -----------------------------
peak_annotation <- annotatePeak(
  peak_gr,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  annoDb = "org.Hs.eg.db"
)

annotated_df <- as.data.frame(peak_annotation)

# -----------------------------
# DEFINE CARDIOMETABOLIC GENES
# -----------------------------
go_terms <- c("GO:0006629","GO:0006006","GO:0042593","GO:0008286")

go_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = go_terms,
  keytype = "GOALL",
  columns = "SYMBOL"
)

cardiometabolic_genes <- unique(go_genes$SYMBOL)
cardio_in_dataset <- intersect(cardiometabolic_genes, annotated_df$SYMBOL)

# -----------------------------
# Build gene-level dataset
# -----------------------------
all_genes_background <- keys(org.Hs.eg.db, keytype = "SYMBOL")

gene_df <- data.frame(Gene = all_genes_background) %>%
  mutate(
    Group = ifelse(Gene %in% cardio_in_dataset, "Cardiometabolic", "Control"),
    HasPeak = ifelse(Gene %in% annotated_df$SYMBOL, 1, 0)
  )

# -----------------------------
# Fisher Test
# -----------------------------
table_2x2 <- table(gene_df$Group, gene_df$HasPeak)
print(table_2x2)

fisher_result <- fisher.test(table_2x2)
print(fisher_result)

# -----------------------------
# PRETTY BAR PLOT
# -----------------------------
plot_df <- gene_df %>%
  group_by(Group) %>%
  summarise(Proportion = mean(HasPeak))

barplot_fig <- ggplot(plot_df, aes(x = Group, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(Proportion * 100, 1)),
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Cardiometabolic" = "#1f78b4",
                               "Control" = "#33a02c")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Chromatin Accessibility Enrichment",
    x = "Gene Group",
    y = "Proportion with ATAC Peaks (%)"
  ) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("barplot_ATAC.png", barplot_fig, width = 6, height = 5, dpi = 300)
print(barplot_fig)

# -----------------------------
# HEATMAP (TOP 20 GENES)
# -----------------------------
cardio_peaks <- annotated_df %>%
  filter(SYMBOL %in% cardio_in_dataset)

cardio_peaks_clean <- cardio_peaks %>%
  mutate(annotation_group = case_when(
    str_detect(annotation, "Promoter") ~ "Promoter",
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron") ~ "Intron",
    TRUE ~ "Intergenic"
  ))

top_genes <- cardio_peaks_clean %>%
  group_by(SYMBOL) %>%
  summarise(total_peaks = n()) %>%
  arrange(desc(total_peaks)) %>%
  slice_head(n = 20) %>%
  pull(SYMBOL)

top_data <- cardio_peaks_clean %>%
  filter(SYMBOL %in% top_genes)

top_counts <- top_data %>%
  count(SYMBOL, annotation_group)

top_matrix <- top_counts %>%
  pivot_wider(
    names_from = annotation_group,
    values_from = n,
    values_fill = 0
  ) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

pheatmap(
  top_matrix,
  scale = "row",
  color = colorRampPalette(c("white", "#1f78b4"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  main = "Top 20 Cardiometabolic Genes (ATAC Peaks)"
)

# -----------------------------
# BOXPLOT + HISTOGRAM (SAFE)
# -----------------------------
peaks_per_gene <- annotated_df %>%
  group_by(SYMBOL) %>%
  summarise(peak_count = n())

plot_data <- gene_df %>%
  left_join(peaks_per_gene, by = c("Gene" = "SYMBOL")) %>%
  mutate(peak_count = ifelse(is.na(peak_count), 0, peak_count))

set.seed(1)
plot_data_sampled <- plot_data %>%
  group_by(Group) %>%
  slice_sample(n = 5000)

# BOXPLOT
boxplot_fig <- ggplot(plot_data_sampled, aes(x = Group, y = peak_count, fill = Group)) +
  geom_boxplot(width = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("Cardiometabolic" = "#1f78b4",
                               "Control" = "#33a02c")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Chromatin Accessibility per Gene",
    x = "Gene Group",
    y = "Number of ATAC Peaks"
  ) +
  theme(legend.position = "none")

ggsave("boxplot_ATAC.png", boxplot_fig, width = 6, height = 5, dpi = 300)
print(boxplot_fig)

# HISTOGRAM
hist_fig <- ggplot(plot_data_sampled, aes(x = peak_count, fill = Group)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  coord_cartesian(xlim = c(0, 100)) +
  scale_fill_manual(values = c("Cardiometabolic" = "#1f78b4",
                               "Control" = "#33a02c")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of ATAC Peaks per Gene",
    x = "Number of ATAC Peaks",
    y = "Number of Genes"
  )

ggsave("histogram_ATAC.png", hist_fig, width = 7, height = 5, dpi = 300)
print(hist_fig)

# -----------------------------
# TOP GENES BAR PLOT (EXTRA)
# -----------------------------
top_genes_plot <- annotated_df %>%
  filter(SYMBOL %in% cardio_in_dataset) %>%
  count(SYMBOL, sort = TRUE) %>%
  slice_head(n = 15)

ggplot(top_genes_plot, aes(x = reorder(SYMBOL, n), y = n)) +
  geom_bar(stat = "identity", fill = "#1f78b4") +
  coord_flip() +
  theme_minimal()

ggsave("top_genes_peaks.png")

# -----------------------------
# ANNOTATION DISTRIBUTION
# -----------------------------
annotation_summary <- cardio_peaks_clean %>%
  count(annotation_group)

ggplot(annotation_summary, aes(x = annotation_group, y = n, fill = annotation_group)) +
  geom_bar(stat = "identity") +
  theme_minimal()

ggsave("annotation_distribution.png")

# -----------------------------
# Save annotated data
# -----------------------------
write.csv(annotated_df,
          "consensus_peak_annotations.csv",
          row.names = FALSE)