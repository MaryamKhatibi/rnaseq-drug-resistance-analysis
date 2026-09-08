# Bulk RNA-seq Analysis of Drug Resistance in Cancer

An end-to-end bioinformatics project investigating transcriptional changes associated with drug resistance in cancer cell models using bulk RNA-seq data.

## Project Overview

Drug resistance is a major challenge in cancer treatment. Changes in gene expression can provide insights into the molecular mechanisms associated with resistance and may reveal candidate genes and biological pathways for further investigation.

This project performs an independent reanalysis of publicly available bulk RNA-seq data from drug-sensitive and drug-resistant cancer cell models.

The analysis is designed as a reproducible RNA-seq workflow, starting from raw sequencing reads and progressing through quality control, read preprocessing, alignment, gene-level quantification, differential expression analysis, visualization, and pathway enrichment.

> **Project status:** In progress

---

## Research Question

**Which genes and biological pathways show altered expression between drug-sensitive and drug-resistant cancer cell models?**

The analysis aims to identify transcriptional patterns associated with drug resistance and characterize the biological processes and pathways represented by differentially expressed genes.

---

## Dataset

The final dataset and study accession will be documented here once the dataset selection is finalized.

* **Repository:** NCBI Gene Expression Omnibus (GEO) / Sequence Read Archive (SRA)
* **Organism:** Human
* **Data type:** Bulk RNA-seq
* **Experimental design:** Drug-sensitive vs. drug-resistant cancer cell models
* **Biological replicates:** To be documented after final dataset selection
* **Raw sequencing data:** FASTQ files
* **Reference genome:** To be specified according to the selected dataset

### Data Availability

Raw sequencing files are not included in this repository because of their large size.

Instructions for obtaining the original data will be provided in:

```text
data/README.md
```

---

## Analysis Workflow

```text
Raw FASTQ files
       │
       ▼
Quality Control
FastQC + MultiQC
       │
       ▼
Adapter / Quality Trimming
Trimmomatic
       │
       ▼
Post-trimming QC
FastQC + MultiQC
       │
       ▼
Read Alignment
HISAT2
       │
       ▼
Alignment Quality Assessment
       │
       ▼
Genome Visualization
IGV
       │
       ▼
Gene-level Quantification
featureCounts
       │
       ▼
Count Matrix
       │
       ▼
Filtering & Normalization
       │
       ▼
Differential Expression Analysis
DESeq2
       │
       ├──────────────┐
       ▼              ▼
PCA            MA / Volcano Plot
       │
       ▼
Heatmap of Differentially
Expressed Genes
       │
       ▼
Functional Enrichment
GO / KEGG / Pathway Analysis
       │
       ▼
Biological Interpretation
```

---

## Tools & Technologies

### Command-line / Linux

* Linux / WSL2
* Bash
* SRA Toolkit
* FastQC
* MultiQC
* Trimmomatic
* HISAT2
* SAMtools
* IGV
* featureCounts

### R / Bioconductor

* R
* RStudio
* Bioconductor
* DESeq2
* ggplot2
* pheatmap
* clusterProfiler
* Additional annotation and pathway-analysis packages as required

### Version Control

* Git
* GitHub

---

## Repository Structure

```text
rnaseq-drug-resistance-analysis/
│
├── data/
│   ├── metadata/
│   └── README.md
│
├── scripts/
│   ├── 01_download.sh
│   ├── 02_fastqc.sh
│   ├── 03_multiqc.sh
│   ├── 04_trimming.sh
│   ├── 05_mapping.sh
│   ├── 06_featurecounts.sh
│   ├── 07_dge.R
│   └── 08_pathway_analysis.R
│
├── results/
│   ├── qc/
│   ├── mapping/
│   ├── counts/
│   ├── differential_expression/
│   └── pathways/
│
├── figures/
│   ├── PCA.png
│   ├── volcano.png
│   ├── MA_plot.png
│   └── heatmap.png
│
├── docs/
│   └── methods.md
│
├── environment/
│   └── software_versions.txt
│
├── README.md
└── LICENSE
```

---

## Quality Control

Raw sequencing reads will first be assessed using FastQC and MultiQC.

Quality metrics will be evaluated for:

* Per-base sequence quality
* Adapter contamination
* Sequence duplication
* GC content
* Sequence length
* Overrepresented sequences

If adapter contamination or low-quality bases are identified, reads will be trimmed using Trimmomatic.

Post-trimming FastQC and MultiQC reports will be generated to evaluate the effect of preprocessing.

---

## Read Alignment

Preprocessed reads will be aligned to the appropriate reference genome using HISAT2.

Alignment quality will be assessed using alignment statistics and relevant QC metrics.

Selected alignments will also be inspected using IGV to visually evaluate read coverage and genomic alignment.

---

## Gene Quantification

Aligned reads will be assigned to genomic features using featureCounts.

The resulting gene-level count matrix will be used as input for downstream differential expression analysis.

---

## Differential Expression Analysis

Differential expression will be performed using DESeq2.

The analysis will include:

* Count filtering
* Library-size normalization
* Exploratory data analysis
* Principal Component Analysis (PCA)
* Differential expression testing
* Multiple-testing correction
* MA plots
* Volcano plots
* Heatmaps

Genes will be evaluated based on statistical significance and effect size rather than statistical significance alone.

---

## Functional and Pathway Analysis

Differentially expressed genes will be investigated using functional enrichment and pathway analysis.

Depending on the final dataset and annotation availability, analyses may include:

* Gene Ontology (GO) enrichment
* KEGG pathway enrichment
* Gene set enrichment analysis
* Pathway visualization

The goal is to identify biological processes and signaling pathways potentially associated with drug resistance.

---

## Reproducibility

All major analysis steps will be implemented using scripts rather than manual execution.

The repository will document:

* Software versions
* Analysis parameters
* Dataset accession
* Reference genome and annotation
* Experimental design
* R package versions
* Processing steps

Large raw sequencing files and intermediate alignment files will not be stored in the GitHub repository.

---

## Results

This section will be updated after completion of the analysis.

### Quality Control

*To be added.*

### Principal Component Analysis

*To be added.*

### Differential Expression

*To be added.*

### Pathway Analysis

*To be added.*

---

## Biological Interpretation

The final results will be interpreted in the context of drug resistance biology.

Because the analysis is based on a specific set of cancer cell models, identified genes and pathways will be considered **candidate transcriptional changes associated with drug resistance**, rather than definitive causal mechanisms.

---

## References

The primary study, GEO accession, SRA BioProject, reference genome, annotation, and major software/tools used in the analysis will be cited here after the final dataset has been selected.

---

## Author

**Maryam Khatibi**

GitHub: [MaryamKhatibi](https://github.com/MaryamKhatibi)
