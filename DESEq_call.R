##### File to set parameters for DESEq2 #####

## Set paths and files
basedir <- "/Users/terezinhasouza/Downloads/"
gene.path <- "/Users/terezinhasouza/Downloads/genes/genes/"
metadata <- "/Users/terezinhasouza/Downloads/metadata.txt" #file name or path to file


## Params for DESeq2
control.name <- "Ctrl" #name of control, or to which contrasts should be made
min.reads <- 10
norm.method <- "vsd"
p.adj.method <- "bonferroni"
save.files <- TRUE #write .csv files for each contrast treated vs control
clustering <- TRUE #clustering of samples
pca <- TRUE #pca of normalized read counts

source("DESeq2_run.R")