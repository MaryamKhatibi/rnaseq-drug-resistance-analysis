# RNA-seq Analysis of Gemcitabine Resistance in Pancreatic Cancer

Independent bulk RNA-seq reanalysis of human pancreatic cancer cell lines to identify transcriptional and pathway-level changes associated with acquired gemcitabine resistance.

## Project overview

This project reanalyzes the public GEO dataset **GSE140077**, which contains transcriptome profiles from two human pancreatic cancer cell lines and their gemcitabine-resistant derivatives. The study includes four experimental conditions with three biological replicates per condition (12 samples total).

The main analysis asks:

> Which transcriptional changes are associated with gemcitabine resistance after accounting for the underlying cell-line background?

The primary statistical model is:

`~ cell_background + resistance`

This separates the effect of resistance status from the strong baseline expression differences between the two cell backgrounds.

## Dataset

- **GEO:** [GSE140077](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140077)
- **BioProject:** PRJNA588180
- **Organism:** *Homo sapiens*
- **Experiment:** bulk RNA-seq
- **Samples:** 12
- **Conditions:**
  - BxPC-3 Sensitive (n=3)
  - BxPC-3-GR Resistant (n=3)
  - CFPAC-1 Sensitive (n=3)
  - CFPAC-1-GR Resistant (n=3)
- **Original publication:** Zhou et al., *Cancer Medicine* (2020), DOI: [10.1002/cam4.2764](https://doi.org/10.1002/cam4.2764)

The original study established gemcitabine-resistant derivatives of BxPC-3 and CFPAC-1 and characterized their transcriptomic profiles. This repository contains an **independent reanalysis** of the public sequencing data rather than a reproduction of the original computational pipeline.

## Analysis workflow

```text
Public SRA data
      |
      v
FastQC / MultiQC
      |
      v
Adapter and quality trimming
(Trimmomatic)
      |
      v
FastQC / MultiQC
      |
      v
Splice-aware alignment
(HISAT2)
      |
      v
BAM sorting and indexing
(Samtools)
      |
      v
Gene-level quantification
(featureCounts)
      |
      v
DESeq2 differential expression
      |
      +---- PCA
      +---- Sample-distance heatmap
      +---- MA plot
      +---- Volcano plot
      |
      v
Pathway analysis
      +---- GO / KEGG over-representation analysis
      +---- GO GSEA
      +---- KEGG GSEA
      +---- Hallmark GSEA
      +---- Leading-edge analysis
```

## Key methods

### Quality control and preprocessing

Raw paired-end reads were evaluated with FastQC and summarized with MultiQC. Reads were adapter-trimmed with Trimmomatic and the paired surviving reads were used for downstream alignment.

### Alignment

Reads were aligned to the human GRCh38 reference genome using HISAT2 with paired-end reverse-stranded RNA-seq settings. BAM files were sorted and indexed with Samtools.

### Quantification

Gene-level counts were generated with featureCounts using the GENCODE v50 annotation and reverse-stranded paired-end settings. Read pairs were counted as fragments/templates.

### Differential expression

DESeq2 was used with the design:

```text
~ cell_background + resistance
```

The primary contrast is:

```text
Resistant vs Sensitive
```

Low-count genes were filtered by requiring at least 10 counts in at least the smallest resistance group. Log2 fold changes were shrunk using `apeglm`.

### Functional analysis

Gene-level results were analyzed using:

- Gene Ontology (GO) Biological Process enrichment
- KEGG over-representation analysis
- GO GSEA
- KEGG GSEA
- MSigDB Hallmark GSEA
- Leading-edge gene extraction for selected enriched gene sets

GSEA rankings were based on the DESeq2 Wald statistic and therefore used the full ranked gene list rather than only statistically significant DEGs.

## Main results

The DESeq2 analysis retained **15,423 genes** after low-count filtering.

Using `padj < 0.05` and `|log2FC| > 1`:

- **2,523 DEGs**
- **1,379 upregulated** in Resistant
- **1,144 downregulated** in Resistant

An additional **8,115 genes** had `padj < 0.05` without imposing the fold-change threshold.

### Over-representation analysis

GO and KEGG over-representation analysis did not identify terms/pathways passing the `FDR < 0.05` threshold with the selected DEG set and background.

### GSEA

GSEA revealed much stronger pathway-level structure when the full ranked gene list was used.

Strong negative enrichment was observed for themes including:

- ribosome biogenesis and rRNA processing
- RNA splicing / spliceosome-related processes
- mitochondrial gene expression and translation
- oxidative phosphorylation and ATP production
- nucleotide biosynthesis

Positive enrichment included signaling, cytoskeletal, and cell-adhesion-related processes.

### Hallmark GSEA

The Hallmark collection provided a compact complementary view of pathway-level changes. Strong negative enrichment included:

- MYC target programs
- oxidative phosphorylation
- DNA repair
- E2F target programs
- unfolded protein response
- xenobiotic metabolism
- epithelial-mesenchymal transition
- p53 pathway
- G2/M checkpoint

Positive enrichment included Hallmark gene sets related to heme metabolism, mitotic spindle, protein secretion, and other regulatory programs.

These findings describe **expression-level associations with resistance** and should not be interpreted as evidence of direct causality.

## Selected figures

The main figures generated by the analysis include:

- `figures/PCA_plot.png`
- `figures/sample_distance_heatmap.png`
- `figures/MA_plot.png`
- `figures/volcano_plot.png`
- `figures/GSEA_GO_BP_top_terms.png`
- `figures/GSEA_Hallmark_top_terms.png`
- `figures/gsea_enrichment_plots/`

## Repository structure

```text
rnaseq-drug-resistance-analysis/
|
├── README.md
├── scripts/
│   ├── prepare_metadata.sh
│   ├── prepare_count_matrix.sh
│   ├── qc_raw.sh
│   ├── multiqc_raw.sh
│   ├── trim_reads.sh
│   ├── qc_trimmed.sh
│   ├── multiqc_trimmed.sh
│   ├── align_hisat2.sh
│   ├── multiqc_mapping.sh
│   ├── featurecounts.sh
│   └── deseq2_analysis.R
│
├── data/
│   ├── metadata/
│   ├── raw/
│   ├── fastq/
│   ├── sra/
│   └── reference/
│
├── results/
│   ├── qc/
│   ├── mapping/
│   ├── counts/
│   └── deseq2/
│
├── figures/
├── docs/
└── environment/
```

Large raw sequencing files, BAM files, genome FASTA files, and HISAT2 index files should not be committed to GitHub. The repository is intended to store code, metadata, selected derived results, figures, and documentation needed to reproduce or understand the analysis.

## Reproducibility notes

The main DESeq2 script expects to be run from the project root:

```r
rm(list = ls())
source("scripts/deseq2_analysis.R")
```

Reference resources used in the analysis include GRCh38 and GENCODE v50. The prebuilt HISAT2 index is also GRCh38, but its underlying Ensembl release differs from the GENCODE annotation release; this provenance difference is documented rather than hidden.

## Interpretation and limitations

This analysis uses established pancreatic cancer cell lines and resistant derivatives, so the results describe transcriptional differences between these experimental cell models. They should not be assumed to represent all pancreatic tumors or all mechanisms of gemcitabine resistance.

The project focuses on association rather than causation. Pathway enrichment indicates coordinated expression patterns and does not by itself demonstrate that a pathway causes resistance.

## Reference

Zhou J, Zhang L, Zheng H, et al. Identification of chemoresistance-related mRNAs based on gemcitabine-resistant pancreatic cancer cell lines. *Cancer Medicine*. 2020;9(3):1115-1130. DOI: https://doi.org/10.1002/cam4.2764

Dataset source: NCBI GEO GSE140077.
