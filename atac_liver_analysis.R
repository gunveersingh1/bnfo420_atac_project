# =========================================
# ATAC Consensus Peak Enrichment Analysis
# =========================================

library(GenomicRanges)
library(ChIPseeker)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)
library(pheatmap)
library(tibble)
library(stringr)


setwd("/Users/gunveersingh/Documents/VCU_BNFO/atac_project")

# -----------------------------
# Load consensus peak BED file
# -----------------------------

peaks_df <- read.table(
  "GSE296875_mergedConsensusPeaks_peaksCalledOnMergedBams_autosomalOnly_withCellTypeOverlap.bed.gz",
  sep = "\t",
  header = FALSE
)

colnames(peaks_df)[1:3] <- c("chr", "start", "end")

# Convert to GRanges
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
# Load cardiometabolic gene list
# -----------------------------

genes <- data.frame(
  GeneSymbol = c(
    "LDLR","APOB","PCSK9","HMGCR","APOE",
    "LPL","CETP","PPARG","PPARA","GCK",
    "GCKR","SREBF1","MTTP","ANGPTL3",
    "ANGPTL4","SORT1","ABCG5","ABCG8",
    "INSR","IRS1","TCF7L2","FTO"
  )
)


# -----------------------------
# Count overlaps
# -----------------------------

cardio_overlap <- sum(annotated_df$SYMBOL %in% genes$GeneSymbol)

total_peaks <- nrow(annotated_df)

# Background: number of unique genes in annotation
all_genes <- unique(annotated_df$SYMBOL)

background_overlap <- sum(all_genes %in% genes$GeneSymbol)

# Build enrichment table
enrichment_matrix <- matrix(
  c(cardio_overlap,
    total_peaks - cardio_overlap,
    background_overlap,
    length(all_genes) - background_overlap),
  nrow = 2
)

# -----------------------------
# Fisher’s Exact Test
# -----------------------------

fisher_result <- fisher.test(enrichment_matrix)

print(fisher_result)

# -----------------------------
# Save results
# -----------------------------

write.csv(annotated_df,
          "consensus_peak_annotations.csv",
          row.names = FALSE)

# -----------------------------
# Visualization
# -----------------------------

prop_cardio <- cardio_overlap / total_peaks
prop_other  <- background_overlap / length(all_genes)

plot_df <- data.frame(
  Category = c("Observed in Peaks", "Expected Background"),
  Proportion = c(prop_cardio, prop_other)
)

ggplot(plot_df, aes(x = Category, y = Proportion)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("Proportion of Cardiometabolic Genes in ATAC Peaks")



# Create simplified annotation category
cardio_peaks_clean <- cardio_peaks %>%
  mutate(annotation_group = case_when(
    str_detect(annotation, "Promoter") ~ "Promoter",
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron") ~ "Intron",
    str_detect(annotation, "Downstream") ~ "Downstream",
    TRUE ~ "Intergenic"
  ))

# Count per gene per simplified category
annotation_counts_clean <- cardio_peaks_clean %>%
  count(SYMBOL, annotation_group)

# Convert to matrix
annotation_matrix_clean <- annotation_counts_clean %>%
  tidyr::pivot_wider(
    names_from = annotation_group,
    values_from = n,
    values_fill = 0
  ) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

# Plot heatmap
pheatmap(
  annotation_matrix_clean,
  color = colorRampPalette(c("white", "blue"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "ATAC Peak Location per Cardiometabolic Gene"
)

# Save
png("clean_heatmap_annotation_per_gene.png", width = 800, height = 800)
pheatmap(
  annotation_matrix_clean,
  color = colorRampPalette(c("white", "blue"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "ATAC Peak Location per Cardiometabolic Gene"
)
dev.off()