# =========================================
# ATAC + STEATOSIS ANALYSIS (FINAL + DRAMATIC HISTOGRAM)
# =========================================

library(GenomicRanges)
library(ChIPseeker)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(tidyr)
library(AnnotationDbi)
library(pheatmap)
library(Seurat)
library(readxl)

setwd("/Users/gunveersingh/Documents/VCU_BNFO/atac_project")

# =========================================
# PART 1: ATAC ANALYSIS
# =========================================

peaks_df <- read.table(
  "GSE296875_mergedConsensusPeaks_peaksCalledOnMergedBams_autosomalOnly_withCellTypeOverlap.bed",
  sep = "\t",
  header = FALSE
)

colnames(peaks_df)[1:3] <- c("chr","start","end")

peak_gr <- GRanges(
  seqnames = peaks_df$chr,
  ranges = IRanges(start = peaks_df$start, end = peaks_df$end)
)

peak_annotation <- annotatePeak(
  peak_gr,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  annoDb = "org.Hs.eg.db"
)

annotated_df <- as.data.frame(peak_annotation)

# =========================================
# CREATE PEAK COUNT PER GENE
# =========================================

gene_peak_counts <- annotated_df %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  summarise(peak_count = n()) %>%
  ungroup()

# =========================================
# DEFINE CARDIOMETABOLIC GENES
# =========================================

go_terms <- c("GO:0006629","GO:0006006","GO:0042593","GO:0008286")

go_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = go_terms,
  keytype = "GOALL",
  columns = "SYMBOL"
)

cardio_genes <- unique(go_genes$SYMBOL)
cardio_in_data <- intersect(cardio_genes, gene_peak_counts$SYMBOL)

# Add group labels
gene_peak_counts <- gene_peak_counts %>%
  mutate(
    Group = ifelse(SYMBOL %in% cardio_in_data, "Cardiometabolic", "Control")
  )

# =========================================
# BOX PLOT + JITTER (BEST VISUAL)
# =========================================

box_plot <- ggplot(gene_peak_counts, aes(x = Group, y = peak_count, fill = Group)) +
  
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  
  geom_jitter(
    width = 0.2,
    alpha = 0.3,
    size = 1,
    color = "black"
  ) +
  
  coord_cartesian(ylim = c(0, 40)) +
  
  scale_fill_manual(values = c(
    "Cardiometabolic" = "#d73027",
    "Control" = "#4575b4"
  )) +
  
  theme_minimal(base_size = 16) +
  
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    legend.position = "none"
  ) +
  
  labs(
    title = "Chromatin Accessibility by Gene Group",
    x = "",
    y = "ATAC Peaks per Gene"
  )

print(box_plot)

ggsave("boxplot_jitter.png", box_plot, width = 7, height = 5)

ggsave("histogram_labeled.png", hist_labeled, width = 7, height = 5)
# =========================================
# FISHER TEST
# =========================================

all_genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")

gene_df <- data.frame(Gene = all_genes) %>%
  mutate(
    Group = ifelse(Gene %in% cardio_in_data, "Cardiometabolic", "Control"),
    HasPeak = ifelse(Gene %in% annotated_df$SYMBOL, 1, 0)
  )

print(fisher.test(table(gene_df$Group, gene_df$HasPeak)))

# =========================================
# PART 2: STEATOSIS
# =========================================

obj <- readRDS("GSE296875_allwells_RNA-ATAC-jointpass-UMI1000-harmonized_WNN_azimuth.rds")
meta <- obj@meta.data
rm(obj); gc()

meta <- meta %>%
  filter(!soc_sample_id %in% c("doublet","ambiguous"))

meta$soc_sample_id <- as.numeric(meta$soc_sample_id)

clinical <- read_excel("1-s2.0-S0002929725004343-mmc2.xlsx", skip = 1)
clinical$sample_id <- as.numeric(clinical$sample_id)

merged_df <- merge(meta, clinical,
                   by.x = "soc_sample_id",
                   by.y = "sample_id")

cat("Merged rows:", nrow(merged_df), "\n")

steatosis_df <- merged_df[, c("nCount_peaks","steatosis_categorical")]
colnames(steatosis_df) <- c("PeakCount","RawGroup")

steatosis_df$RawGroup <- tolower(as.character(steatosis_df$RawGroup))

steatosis_df$Group <- NA
steatosis_df$Group[grepl("none", steatosis_df$RawGroup)] <- "None"
steatosis_df$Group[grepl("severe", steatosis_df$RawGroup)] <- "Severe"
steatosis_df$Group[grepl("mild|moderate", steatosis_df$RawGroup)] <- "Mild/Moderate"

steatosis_df <- steatosis_df[!is.na(steatosis_df$Group), ]

set.seed(1)

split_list <- split(steatosis_df, steatosis_df$Group)

sampled_list <- lapply(split_list, function(df) {
  df[sample(nrow(df), size = min(3000, nrow(df)), replace = TRUE), ]
})

steatosis_df_sampled <- do.call(rbind, sampled_list)

steatosis_plot <- ggplot(steatosis_df_sampled,
                         aes(x = Group, y = PeakCount, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "None"="#33a02c",
    "Mild/Moderate"="#1f78b4",
    "Severe"="#e31a1c"
  )) +
  theme_minimal() +
  labs(title="Chromatin Accessibility vs Steatosis Severity")

print(steatosis_plot)

ggsave("steatosis_boxplot.png", steatosis_plot)

anova_result <- aov(PeakCount ~ Group, data = steatosis_df_sampled)
print(summary(anova_result))