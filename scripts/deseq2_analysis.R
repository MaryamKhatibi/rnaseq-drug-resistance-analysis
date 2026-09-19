# ============================================================
# Differential Expression Analysis with DESeq2
# ============================================================
#
# Project : RNA-seq analysis of gemcitabine resistance
# Dataset : GSE140077
# Organism: Homo sapiens
#
# Main objective:
# Identify transcriptional changes associated with gemcitabine
# resistance while accounting for cell-background effects.
#
# Main design:
#   ~ cell_background + resistance
#
# Author: Maryam Khatibi
# ============================================================

# ============================================================
# 1. Load packages
# ============================================================

library(DESeq2)
library(apeglm)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(msigdbr)

# ============================================================
# 2. Define project paths
# ============================================================

# Run this script from the project root directory.
project_dir <- normalizePath(getwd())

count_file <- file.path(
  project_dir,
  "results/counts/count_matrix.txt"
)

metadata_file <- file.path(
  project_dir,
  "data/metadata/sample_metadata.tsv"
)

results_dir <- file.path(
  project_dir,
  "results/deseq2"
)

figures_dir <- file.path(
  project_dir,
  "figures"
)

# Check that the script is being run from the project root.
if (!file.exists(count_file)) {
  stop(
    "Count matrix not found. Run this script from the project root."
  )
}

if (!file.exists(metadata_file)) {
  stop(
    "Sample metadata not found. Run this script from the project root."
  )
}

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 3. Import count matrix and sample metadata
# ============================================================

countData <- read.delim(
  count_file,
  row.names = 1,
  check.names = FALSE
)

colData <- read.delim(
  metadata_file,
  row.names = 1,
  check.names = FALSE
)

# ============================================================
# 4. Verify sample matching
# ============================================================

stopifnot(
  all(rownames(colData) %in% colnames(countData))
)

stopifnot(
  all(colnames(countData) %in% rownames(colData))
)

# Reorder count matrix columns to exactly match metadata rows.
countData <- countData[, rownames(colData)]

stopifnot(
  identical(colnames(countData), rownames(colData))
)

# ============================================================
# 5. Prepare experimental factors
# ============================================================

colData$cell_background <- factor(
  colData$cell_background,
  levels = c("BxPC-3", "CFPAC-1")
)

colData$resistance <- factor(
  colData$resistance,
  levels = c("Sensitive", "Resistant")
)

# ============================================================
# 6. Create DESeq2 dataset
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData,
  design    = ~ cell_background + resistance
)

# Confirm the experimental design.
print(design(dds))

# ============================================================
# 7. Filter low-count genes
# ============================================================

min_group_size <- min(
  table(colData$resistance)
)

keep <- rowSums(
  counts(dds) >= 10
) >= min_group_size

dds <- dds[keep, ]

cat(
  "Genes retained after filtering:",
  nrow(dds),
  "\n"
)

# ============================================================
# 8. Run DESeq2
# ============================================================

dds <- DESeq(dds)

# ============================================================
# 9. Differential expression: Resistant vs Sensitive
# ============================================================

res <- results(
  dds,
  contrast = c(
    "resistance",
    "Resistant",
    "Sensitive"
  ),
  alpha = 0.05
)

# ============================================================
# 10. Shrink log2 fold changes
# ============================================================

res_shrunk <- lfcShrink(
  dds,
  coef = "resistance_Resistant_vs_Sensitive",
  type = "apeglm"
)

# ============================================================
# 11. Prepare and annotate differential expression results
# ============================================================

res_df <- as.data.frame(res_shrunk)

# Preserve the original versioned Ensembl gene ID.
res_df$gene_id <- rownames(res_df)

# Remove Ensembl version suffix for database annotation.
rownames(res_df) <- sub(
  "\\..*$",
  "",
  rownames(res_df)
)

# Map Ensembl IDs to human gene annotation.
res_df$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys      = rownames(res_df),
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

res_df$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys      = rownames(res_df),
  column    = "ENTREZID",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

res_df$GENENAME <- mapIds(
  org.Hs.eg.db,
  keys      = rownames(res_df),
  column    = "GENENAME",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# ============================================================
# 12. Identify and classify differentially expressed genes
# ============================================================

significant_genes <- subset(
  res_df,
  !is.na(padj) &
    padj < 0.05
)

DEGs <- subset(
  res_df,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) > 1
)

DEGs_up <- subset(
  DEGs,
  log2FoldChange > 1
)

DEGs_down <- subset(
  DEGs,
  log2FoldChange < -1
)

cat(
  "Significant genes (padj < 0.05):",
  nrow(significant_genes),
  "\n"
)

cat(
  "DEGs (padj < 0.05 and |log2FC| > 1):",
  nrow(DEGs),
  "\n"
)

cat(
  "Upregulated DEGs:",
  nrow(DEGs_up),
  "\n"
)

cat(
  "Downregulated DEGs:",
  nrow(DEGs_down),
  "\n"
)

# ============================================================
# 13. Export differential expression results
# ============================================================

write.csv(
  res_df,
  file = file.path(
    results_dir,
    "DESeq2_all_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  DEGs,
  file = file.path(
    results_dir,
    "DEGs.csv"
  ),
  row.names = FALSE
)

write.csv(
  DEGs_up,
  file = file.path(
    results_dir,
    "DEGs_upregulated.csv"
  ),
  row.names = FALSE
)

write.csv(
  DEGs_down,
  file = file.path(
    results_dir,
    "DEGs_downregulated.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 14. Variance stabilizing transformation
# ============================================================

vsd <- vst(
  dds,
  blind = FALSE
)

# ============================================================
# 15. Principal component analysis
# ============================================================

pca_data <- plotPCA(
  vsd,
  intgroup = c(
    "cell_background",
    "resistance"
  ),
  returnData = TRUE
)

percentVar <- round(
  100 * attr(
    pca_data,
    "percentVar"
  )
)

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = resistance,
    shape = cell_background
  )
) +
  geom_point(size = 4) +
  xlab(
    paste0(
      "PC1: ",
      percentVar[1],
      "% variance"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      percentVar[2],
      "% variance"
    )
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    legend.position = "right"
  )

ggsave(
  filename = file.path(
    figures_dir,
    "PCA_plot.png"
  ),
  plot = pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# ============================================================
# 16. Sample distance heatmap
# ============================================================

sample_dist <- dist(
  t(assay(vsd))
)

sample_dist_matrix <- as.matrix(
  sample_dist
)

rownames(sample_dist_matrix) <- colnames(vsd)
colnames(sample_dist_matrix) <- colnames(vsd)

sample_annotation <- as.data.frame(
  colData[, c(
    "cell_background",
    "resistance"
  )]
)

pheatmap(
  sample_dist_matrix,
  annotation_col = sample_annotation,
  annotation_row = sample_annotation,
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  main = "Sample-to-sample distance",
  filename = file.path(
    figures_dir,
    "sample_distance_heatmap.png"
  ),
  width = 8,
  height = 7
)

# ============================================================
# 17. MA plot
# ============================================================

png(
  filename = file.path(
    figures_dir,
    "MA_plot.png"
  ),
  width = 2000,
  height = 1600,
  res = 300
)

plotMA(
  res_shrunk,
  alpha = 0.05,
  ylim = c(-5, 5),
  main = "MA plot: Resistant vs Sensitive"
)

dev.off()

# ============================================================
# 18. Volcano plot
# ============================================================

volcano_data <- res_df

volcano_data$significance <- "Not significant"

volcano_data$significance[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange > 1
] <- "Upregulated"

volcano_data$significance[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange < -1
] <- "Downregulated"

volcano_data$neg_log10_padj <- -log10(
  volcano_data$padj
)

volcano_data <- volcano_data[
  is.finite(volcano_data$neg_log10_padj),
]

volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj
  )
) +
  geom_point(
    aes(color = significance),
    alpha = 0.6,
    size = 1.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Volcano plot: Resistant vs Sensitive",
    x = "Log2 fold change",
    y = "-Log10 adjusted p-value",
    color = "Classification"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12)
  )

ggsave(
  filename = file.path(
    figures_dir,
    "volcano_plot.png"
  ),
  plot = volcano_plot,
  width = 7,
  height = 6,
  dpi = 300
)

# ============================================================
# 19. GO and KEGG over-representation analysis
# ============================================================

# Prepare DEG and background gene sets.
geneList <- unique(
  na.omit(DEGs$ENTREZID)
)

universe <- unique(
  na.omit(res_df$ENTREZID)
)

cat(
  "Genes used for enrichment:",
  length(geneList),
  "\n"
)

cat(
  "Background genes:",
  length(universe),
  "\n"
)

# ------------------------------------------------------------
# 19.1 GO Biological Process enrichment
# ------------------------------------------------------------

ego_BP <- enrichGO(
  gene          = geneList,
  universe      = universe,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(
  as.data.frame(ego_BP),
  file = file.path(
    results_dir,
    "GO_BP_enrichment.csv"
  ),
  row.names = FALSE
)

cat(
  "Significant GO Biological Process terms:",
  nrow(as.data.frame(ego_BP)),
  "\n"
)

if (nrow(as.data.frame(ego_BP)) > 0) {

  go_dotplot <- dotplot(
    ego_BP,
    showCategory = 20
  ) +
    ggtitle("GO Biological Process enrichment")

  ggsave(
    filename = file.path(
      figures_dir,
      "GO_BP_dotplot.png"
    ),
    plot = go_dotplot,
    width = 9,
    height = 7,
    dpi = 300
  )

  go_barplot <- barplot(
    ego_BP,
    showCategory = 20
  ) +
    ggtitle("GO Biological Process enrichment")

  ggsave(
    filename = file.path(
      figures_dir,
      "GO_BP_barplot.png"
    ),
    plot = go_barplot,
    width = 9,
    height = 7,
    dpi = 300
  )
}

# ------------------------------------------------------------
# 19.2 KEGG enrichment
# ------------------------------------------------------------

ekegg <- enrichKEGG(
  gene         = geneList,
  universe     = universe,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

write.csv(
  as.data.frame(ekegg),
  file = file.path(
    results_dir,
    "KEGG_Pathway.csv"
  ),
  row.names = FALSE
)

cat(
  "Significant KEGG pathways:",
  nrow(as.data.frame(ekegg)),
  "\n"
)

if (nrow(as.data.frame(ekegg)) > 0) {

  kegg_dotplot <- dotplot(
    ekegg,
    showCategory = 20
  ) +
    ggtitle("KEGG pathway enrichment")

  ggsave(
    filename = file.path(
      figures_dir,
      "KEGG_dotplot.png"
    ),
    plot = kegg_dotplot,
    width = 9,
    height = 7,
    dpi = 300
  )
}

# ------------------------------------------------------------
# 19.3 Prepare fold-change data for Pathview
# ------------------------------------------------------------

fc_data <- res_df[
  !is.na(res_df$ENTREZID) &
    !is.na(res_df$log2FoldChange),
  c("ENTREZID", "log2FoldChange")
]

# Aggregate duplicated Entrez IDs by mean log2 fold change.
foldChanges <- tapply(
  fc_data$log2FoldChange,
  fc_data$ENTREZID,
  mean,
  na.rm = TRUE
)

# ------------------------------------------------------------
# 19.4 Optional Pathview visualization
# ------------------------------------------------------------

# Pathview is generated only when significant KEGG ORA pathways
# are available. For the current dataset, KEGG ORA returned no
# significant pathways, so this section produces no figures.

if (nrow(as.data.frame(ekegg)) > 0) {

  top_pathways <- head(
    ekegg$ID,
    20
  )

  pathview_dir <- file.path(
    figures_dir,
    "pathview_output"
  )

  dir.create(
    pathview_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  pv.out <- lapply(
    top_pathways,
    function(pid) {

      pathview(
        gene.data   = foldChanges,
        gene.idtype = "ENTREZ",
        pathway.id  = pid,
        species     = "hsa",
        kegg.dir    = pathview_dir
      )
    }
  )
}

# ============================================================
# 20. Prepare ranked gene list for GSEA
# ============================================================

# Use the DESeq2 Wald statistic to rank the complete gene set.
# This preserves information across the full expression ranking.

gsea_data <- data.frame(
  ensembl_id = sub(
    "\\..*$",
    "",
    rownames(res)
  ),
  stat = res$stat,
  stringsAsFactors = FALSE
)

# Map Ensembl IDs to Entrez IDs.
gsea_data$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys      = gsea_data$ensembl_id,
  column    = "ENTREZID",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# Remove missing and non-finite values.
gsea_data <- gsea_data[
  !is.na(gsea_data$ENTREZID) &
    !is.na(gsea_data$stat) &
    is.finite(gsea_data$stat),
]

# Keep the strongest Wald statistic when multiple Ensembl IDs
# map to the same Entrez ID.
gsea_data <- gsea_data[
  order(
    abs(gsea_data$stat),
    decreasing = TRUE
  ),
]

gsea_data <- gsea_data[
  !duplicated(gsea_data$ENTREZID),
]

gene_rank <- gsea_data$stat
names(gene_rank) <- gsea_data$ENTREZID

gene_rank <- sort(
  gene_rank,
  decreasing = TRUE
)

cat(
  "Genes in GSEA ranking:",
  length(gene_rank),
  "\n"
)

# ============================================================
# 21. GO Biological Process GSEA
# ============================================================

gsea_GO <- gseGO(
  geneList      = gene_rank,
  ont           = "BP",
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = TRUE
)

gsea_GO_df <- as.data.frame(
  gsea_GO
)

write.csv(
  gsea_GO_df,
  file = file.path(
    results_dir,
    "GSEA_GO_BP.csv"
  ),
  row.names = FALSE
)

cat(
  "Significant GO GSEA terms:",
  nrow(gsea_GO_df),
  "\n"
)

if (nrow(gsea_GO_df) > 0) {

  gsea_GO_dotplot <- dotplot(
    gsea_GO,
    showCategory = 20
  ) +
    ggtitle("GSEA: GO Biological Process")

  ggsave(
    filename = file.path(
      figures_dir,
      "GSEA_GO_BP_dotplot.png"
    ),
    plot = gsea_GO_dotplot,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ============================================================
# 22. KEGG GSEA
# ============================================================

gsea_KEGG <- gseKEGG(
  geneList      = gene_rank,
  organism      = "hsa",
  keyType       = "ncbi-geneid",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = TRUE
)

gsea_KEGG_df <- as.data.frame(
  gsea_KEGG
)

write.csv(
  gsea_KEGG_df,
  file = file.path(
    results_dir,
    "GSEA_KEGG.csv"
  ),
  row.names = FALSE
)

cat(
  "Significant KEGG GSEA pathways:",
  nrow(gsea_KEGG_df),
  "\n"
)

if (nrow(gsea_KEGG_df) > 0) {

  gsea_KEGG_dotplot <- dotplot(
    gsea_KEGG,
    showCategory = 20
  ) +
    ggtitle("GSEA: KEGG pathways")

  ggsave(
    filename = file.path(
      figures_dir,
      "GSEA_KEGG_dotplot.png"
    ),
    plot = gsea_KEGG_dotplot,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ============================================================
# 23. Simplify GO GSEA results
# ============================================================

# Remove redundant GO terms based on semantic similarity.
gsea_GO_simplified <- simplify(
  gsea_GO,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

gsea_GO_simplified_df <- as.data.frame(
  gsea_GO_simplified
)

write.csv(
  gsea_GO_simplified_df,
  file = file.path(
    results_dir,
    "GSEA_GO_BP_simplified.csv"
  ),
  row.names = FALSE
)

cat(
  "GO terms before simplification:",
  nrow(gsea_GO_df),
  "\n"
)

cat(
  "GO terms after simplification:",
  nrow(gsea_GO_simplified_df),
  "\n"
)

if (nrow(gsea_GO_simplified_df) > 0) {

  gsea_GO_simplified_plot <- dotplot(
    gsea_GO_simplified,
    showCategory = 20
  ) +
    ggtitle("GSEA: Simplified GO Biological Processes")

  ggsave(
    filename = file.path(
      figures_dir,
      "GSEA_GO_BP_simplified_dotplot.png"
    ),
    plot = gsea_GO_simplified_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ============================================================
# 24. Separate significant GO GSEA terms by direction
# ============================================================

go_simple_df <- as.data.frame(
  gsea_GO_simplified
)

go_positive <- subset(
  go_simple_df,
  NES > 0 &
    p.adjust < 0.05
)

go_negative <- subset(
  go_simple_df,
  NES < 0 &
    p.adjust < 0.05
)

cat(
  "Significant positive GO terms:",
  nrow(go_positive),
  "\n"
)

cat(
  "Significant negative GO terms:",
  nrow(go_negative),
  "\n"
)

# ============================================================
# 25. Visualize top GO GSEA terms
# ============================================================

top_go_positive <- head(
  go_positive[
    order(go_positive$p.adjust),
    c(
      "ID",
      "Description",
      "NES",
      "p.adjust"
    )
  ],
  5
)

top_go_negative <- head(
  go_negative[
    order(go_negative$p.adjust),
    c(
      "ID",
      "Description",
      "NES",
      "p.adjust"
    )
  ],
  5
)

top_go <- rbind(
  top_go_positive,
  top_go_negative
)

top_go$Description <- factor(
  top_go$Description,
  levels = top_go$Description[
    order(top_go$NES)
  ]
)

top_go$neg_log10_padj <- -log10(
  top_go$p.adjust
)

if (nrow(top_go) > 0) {

  gsea_go_summary_plot <- ggplot(
    top_go,
    aes(
      x = NES,
      y = Description,
      size = neg_log10_padj,
      color = p.adjust
    )
  ) +
    geom_point() +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    labs(
      title = "Top GO Biological Processes from GSEA",
      x = "Normalized Enrichment Score (NES)",
      y = "Biological Process",
      size = "-Log10 adjusted p-value",
      color = "Adjusted p-value"
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      axis.text.y = element_text(size = 10)
    )

  ggsave(
    filename = file.path(
      figures_dir,
      "GSEA_GO_BP_top_terms.png"
    ),
    plot = gsea_go_summary_plot,
    width = 10,
    height = 8,
    dpi = 300
  )

  write.csv(
    top_go,
    file = file.path(
      results_dir,
      "GSEA_GO_BP_top_terms.csv"
    ),
    row.names = FALSE
  )
}

# ============================================================
# 26. Hallmark Gene Set Enrichment Analysis
# ============================================================

hallmark_sets <- msigdbr(
  db_species = "HS",
  species = "human",
  collection = "H"
)

cat(
  "Number of Hallmark gene sets:",
  length(unique(hallmark_sets$gs_name)),
  "\n"
)

hallmark_term2gene <- hallmark_sets[
  ,
  c(
    "gs_name",
    "ncbi_gene"
  )
]

hallmark_term2gene <- hallmark_term2gene[
  !is.na(hallmark_term2gene$ncbi_gene),
]

hallmark_term2gene$ncbi_gene <- as.character(
  hallmark_term2gene$ncbi_gene
)

hallmark_term2name <- unique(
  hallmark_sets[
    ,
    c(
      "gs_name",
      "gs_description"
    )
  ]
)

gsea_Hallmark <- GSEA(
  geneList      = gene_rank,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  TERM2GENE     = hallmark_term2gene,
  TERM2NAME     = hallmark_term2name,
  verbose       = TRUE
)

gsea_Hallmark_df <- as.data.frame(
  gsea_Hallmark
)

write.csv(
  gsea_Hallmark_df,
  file = file.path(
    results_dir,
    "GSEA_Hallmark.csv"
  ),
  row.names = FALSE
)

cat(
  "Significant Hallmark gene sets:",
  nrow(gsea_Hallmark_df),
  "\n"
)

if (nrow(gsea_Hallmark_df) > 0) {

  hallmark_dotplot <- dotplot(
    gsea_Hallmark,
    showCategory = 20
  ) +
    ggtitle("GSEA: Hallmark gene sets")

  ggsave(
    filename = file.path(
      figures_dir,
      "GSEA_Hallmark_dotplot.png"
    ),
    plot = hallmark_dotplot,
    width = 10,
    height = 8,
    dpi = 300
  )
}

hallmark_positive <- subset(
  gsea_Hallmark_df,
  NES > 0 &
    p.adjust < 0.05
)

hallmark_negative <- subset(
  gsea_Hallmark_df,
  NES < 0 &
    p.adjust < 0.05
)

cat(
  "Significant positive Hallmark sets:",
  nrow(hallmark_positive),
  "\n"
)

cat(
  "Significant negative Hallmark sets:",
  nrow(hallmark_negative),
  "\n"
)

# ============================================================
# 27. Visualize top Hallmark GSEA results
# ============================================================

top_hallmark_positive <- head(
  hallmark_positive[
    order(hallmark_positive$p.adjust),
    c(
      "ID",
      "Description",
      "NES",
      "p.adjust"
    )
  ],
  5
)

top_hallmark_negative <- head(
  hallmark_negative[
    order(hallmark_negative$p.adjust),
    c(
      "ID",
      "Description",
      "NES",
      "p.adjust"
    )
  ],
  5
)

top_hallmark <- rbind(
  top_hallmark_positive,
  top_hallmark_negative
)

top_hallmark$Hallmark <- sub(
  "^HALLMARK_",
  "",
  top_hallmark$ID
)

top_hallmark$Hallmark <- gsub(
  "_",
  " ",
  top_hallmark$Hallmark
)

top_hallmark$Hallmark <- factor(
  top_hallmark$Hallmark,
  levels = top_hallmark$Hallmark[
    order(top_hallmark$NES)
  ]
)

top_hallmark$neg_log10_padj <- -log10(
  top_hallmark$p.adjust
)

if (nrow(top_hallmark) > 0) {

  hallmark_summary_plot <- ggplot(
    top_hallmark,
    aes(
      x = NES,
      y = Hallmark,
      size = neg_log10_padj,
      color = p.adjust
    )
  ) +
    geom_point() +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    labs(
      title = "Top Hallmark Gene Sets from GSEA",
      x = "Normalized Enrichment Score (NES)",
      y = "Hallmark gene set",
      size = "-Log10 adjusted p-value",
      color = "Adjusted p-value"
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(
        size = 16,
        face = "bold"
      )
    )

  ggsave(
    filename = file.path(
      figures_dir,
      "GSEA_Hallmark_top_terms.png"
    ),
    plot = hallmark_summary_plot,
    width = 9,
    height = 7,
    dpi = 300
  )

  write.csv(
    top_hallmark,
    file = file.path(
      results_dir,
      "GSEA_Hallmark_top_terms.csv"
    ),
    row.names = FALSE
  )
}

# ============================================================
# 28. Detailed GSEA enrichment plots
# ============================================================

gsea_plots_dir <- file.path(
  figures_dir,
  "gsea_enrichment_plots"
)

dir.create(
  gsea_plots_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 28.1 Prepare plotting copies with concise labels
# ------------------------------------------------------------

gsea_Hallmark_plot <- gsea_Hallmark

gsea_Hallmark_plot@result$Description[
  gsea_Hallmark_plot@result$ID ==
    "HALLMARK_MYC_TARGETS_V2"
] <- "MYC TARGETS V2"

gsea_Hallmark_plot@result$Description[
  gsea_Hallmark_plot@result$ID ==
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
] <- "OXIDATIVE PHOSPHORYLATION"

gsea_KEGG_plot <- gsea_KEGG

gsea_KEGG_plot@result$Description[
  gsea_KEGG_plot@result$ID == "hsa03040"
] <- "Spliceosome"

gsea_KEGG_plot@result$Description[
  gsea_KEGG_plot@result$ID == "hsa00190"
] <- "Oxidative Phosphorylation"

# ------------------------------------------------------------
# 28.2 Plot Hallmark MYC Targets V2
# ------------------------------------------------------------

hallmark_myc_id <- which(
  gsea_Hallmark_plot@result$ID ==
    "HALLMARK_MYC_TARGETS_V2"
)

if (length(hallmark_myc_id) == 1) {

  hallmark_myc_plot <- gseaplot2(
    gsea_Hallmark_plot,
    geneSetID = hallmark_myc_id,
    title = "GSEA: Hallmark MYC Targets V2",
    pvalue_table = TRUE
  )

  ggsave(
    filename = file.path(
      gsea_plots_dir,
      "GSEA_Hallmark_MYC_TARGETS_V2.png"
    ),
    plot = hallmark_myc_plot,
    width = 9,
    height = 7,
    dpi = 300
  )
}

# ------------------------------------------------------------
# 28.3 Plot Hallmark oxidative phosphorylation
# ------------------------------------------------------------

hallmark_oxphos_id <- which(
  gsea_Hallmark_plot@result$ID ==
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)

if (length(hallmark_oxphos_id) == 1) {

  hallmark_oxphos_plot <- gseaplot2(
    gsea_Hallmark_plot,
    geneSetID = hallmark_oxphos_id,
    title = "GSEA: Hallmark Oxidative Phosphorylation",
    pvalue_table = TRUE
  )

  ggsave(
    filename = file.path(
      gsea_plots_dir,
      "GSEA_Hallmark_OXIDATIVE_PHOSPHORYLATION.png"
    ),
    plot = hallmark_oxphos_plot,
    width = 9,
    height = 7,
    dpi = 300
  )
}

# ------------------------------------------------------------
# 28.4 Plot KEGG Spliceosome
# ------------------------------------------------------------

kegg_spliceosome_id <- which(
  gsea_KEGG_plot@result$ID ==
    "hsa03040"
)

if (length(kegg_spliceosome_id) == 1) {

  kegg_spliceosome_plot <- gseaplot2(
    gsea_KEGG_plot,
    geneSetID = kegg_spliceosome_id,
    title = "GSEA: KEGG Spliceosome",
    pvalue_table = TRUE
  )

  ggsave(
    filename = file.path(
      gsea_plots_dir,
      "GSEA_KEGG_Spliceosome.png"
    ),
    plot = kegg_spliceosome_plot,
    width = 9,
    height = 7,
    dpi = 300
  )
}

# ------------------------------------------------------------
# 28.5 Plot KEGG oxidative phosphorylation
# ------------------------------------------------------------

kegg_oxphos_id <- which(
  gsea_KEGG_plot@result$ID ==
    "hsa00190"
)

if (length(kegg_oxphos_id) == 1) {

  kegg_oxphos_plot <- gseaplot2(
    gsea_KEGG_plot,
    geneSetID = kegg_oxphos_id,
    title = "GSEA: KEGG Oxidative Phosphorylation",
    pvalue_table = TRUE
  )

  ggsave(
    filename = file.path(
      gsea_plots_dir,
      "GSEA_KEGG_Oxidative_Phosphorylation.png"
    ),
    plot = kegg_oxphos_plot,
    width = 9,
    height = 7,
    dpi = 300
  )
}

cat(
  "Detailed GSEA plots generated.",
  "\n"
)

# ============================================================
# 29. Leading-edge gene analysis
# ============================================================

# Leading-edge genes are the core genes contributing to the
# enrichment signal of selected gene sets.

selected_pathways <- list(
  Hallmark_MYC_Targets_V2 = list(
    object = gsea_Hallmark,
    id = "HALLMARK_MYC_TARGETS_V2"
  ),
  Hallmark_Oxidative_Phosphorylation = list(
    object = gsea_Hallmark,
    id = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
  ),
  KEGG_Spliceosome = list(
    object = gsea_KEGG,
    id = "hsa03040"
  ),
  KEGG_Oxidative_Phosphorylation = list(
    object = gsea_KEGG,
    id = "hsa00190"
  )
)

leading_edge_list <- list()

for (pathway_name in names(selected_pathways)) {

  gsea_object <- selected_pathways[[pathway_name]]$object
  pathway_id <- selected_pathways[[pathway_name]]$id

  pathway_result <- as.data.frame(
    gsea_object
  )

  pathway_row <- pathway_result[
    pathway_result$ID == pathway_id,
    ,
    drop = FALSE
  ]

  if (nrow(pathway_row) != 1) {

    warning(
      "Pathway not found or duplicated: ",
      pathway_id
    )

    next
  }

  leading_genes <- unlist(
    strsplit(
      pathway_row$core_enrichment,
      split = "/"
    )
  )

  leading_edge_list[[pathway_name]] <- unique(
    leading_genes
  )
}

if (length(leading_edge_list) > 0) {

  # ----------------------------------------------------------
  # 29.1 Create leading-edge results table
  # ----------------------------------------------------------

  leading_edge_df <- do.call(
    rbind,
    lapply(
      names(leading_edge_list),
      function(pathway_name) {

        data.frame(
          pathway = pathway_name,
          ENTREZID = leading_edge_list[[pathway_name]],
          stringsAsFactors = FALSE
        )
      }
    )
  )

  rownames(leading_edge_df) <- NULL

  # ----------------------------------------------------------
  # 29.2 Annotate leading-edge genes
  # ----------------------------------------------------------

  leading_edge_df$SYMBOL <- mapIds(
    org.Hs.eg.db,
    keys = leading_edge_df$ENTREZID,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )

  leading_edge_df$GENENAME <- mapIds(
    org.Hs.eg.db,
    keys = leading_edge_df$ENTREZID,
    column = "GENENAME",
    keytype = "ENTREZID",
    multiVals = "first"
  )

  # ----------------------------------------------------------
  # 29.3 Export leading-edge genes
  # ----------------------------------------------------------

  write.csv(
    leading_edge_df,
    file = file.path(
      results_dir,
      "GSEA_leading_edge_genes.csv"
    ),
    row.names = FALSE
  )

  leading_edge_counts <- as.data.frame(
    table(leading_edge_df$pathway)
  )

  colnames(leading_edge_counts) <- c(
    "Pathway",
    "LeadingEdgeGeneCount"
  )

  write.csv(
    leading_edge_counts,
    file = file.path(
      results_dir,
      "GSEA_leading_edge_counts.csv"
    ),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # 29.4 Identify shared leading-edge genes
  # ----------------------------------------------------------

  gene_occurrence <- table(
    leading_edge_df$ENTREZID
  )

  shared_genes <- names(
    gene_occurrence[
      gene_occurrence > 1
    ]
  )

  shared_leading_edge_df <- leading_edge_df[
    leading_edge_df$ENTREZID %in% shared_genes,
    ,
    drop = FALSE
  ]

  write.csv(
    shared_leading_edge_df,
    file = file.path(
      results_dir,
      "GSEA_shared_leading_edge_genes.csv"
    ),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # 29.5 Summarize shared leading-edge genes
  # ----------------------------------------------------------

  if (nrow(shared_leading_edge_df) > 0) {

    shared_gene_summary <- aggregate(
      pathway ~ ENTREZID + SYMBOL + GENENAME,
      data = shared_leading_edge_df,
      FUN = paste,
      collapse = "; "
    )

    shared_gene_summary$PathwayCount <- sapply(
      strsplit(
        shared_gene_summary$pathway,
        split = "; "
      ),
      length
    )

    shared_gene_summary <- shared_gene_summary[
      order(
        shared_gene_summary$PathwayCount,
        decreasing = TRUE
      ),
    ]

    write.csv(
      shared_gene_summary,
      file = file.path(
        results_dir,
        "GSEA_shared_leading_edge_summary.csv"
      ),
      row.names = FALSE
    )

    cat(
      "Genes shared by multiple selected gene sets:",
      length(shared_genes),
      "\n"
    )
  }
}

# ============================================================
# End of analysis
# ============================================================
