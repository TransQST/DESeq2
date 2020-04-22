# Differential gene expression with RNA-seq data (DESEq2 R package)
> Originally developed for use in the TransQST IMI project (transqst.org). All the scripts are written in R. 
> It is suited for toxicogenomics designs with multiple doses and time points.


## Prerequisites
Datasets: inputs are **RSEM output tables (ending in ".genes.results")** and a **metadata file** containing multiple columns (filenames that match table names), dose and time attributes.

Necessary packages:
```{r}
DESeq2
RColorBrewer
pheatmap
ggplot2
```
## Outline

### Input files
RSEM outputs (*.genes.results) containing expected number of read counts and a metadata file.
Example files (example.genes.results and example_metadata.txt) are provided.

### Usage
There are two sets of scripts, ```DESeq2_call.R``` and ```DESeq2_run.R```.
```DESeq2_call.R``` is customized by the user with input parameters (files, paths, type of normalization and outputs).
```DESeq2_run.R``` is sourced by ```DESeq2_call.R``` to run the analysis and produce outputs. No customization is needed unles the user wants to change specific parameters (e.g., use TPM instead of expected read counts).

### Outputs
1. Table with normalized read counts.
2. Tables containing differentially expressed genes for each contrast.
3. A PCA and a hierarchical clustering of normalized data.
