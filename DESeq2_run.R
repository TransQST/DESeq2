###### Script to run th parameters passed to DESeq_call.R ####

# Check to see if packages are installed. 
#Install them if they are not, then load them into the R session.
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("RColorBrewer", "pheatmap", "ggplot2", "dplyr")
check.packages(packages)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")

setwd(basedir)
## Read and create a copy of metadata file
meta <- read.delim(metadata, header = T, check.names = F, 
                       sep = "\t", stringsAsFactors = T)
meta2 <- meta
meta2$condition <- paste(meta$dose, meta$time, sep = "_")
row.names(meta2) <- paste(meta$dose, meta$time, meta$replicate, sep = "_")
##Set input for contrasts and order factors
times <- unique(meta2$time) %>% factor(., levels = sort(.))
doses <- setdiff(unique(meta2$dose), control.name) %>% factor(., levels = c(control.name, sort(.)))

# Check if files are in place
if(!all(paste0(meta2$filename, ".genes.results") %in% 
  list.files(gene.path))){
  stop("Metadata and input data files do not match")
}
cat("Metadata and input data files match, moving on...")

## Merge genes.results files into a single file
setwd(gene.path)
file.names <- meta2$filename
files.list <- list()
cat("Merging 'genes.results' files")
for(file in meta2$filename){
  read.file <- read.delim(paste0(file, ".genes.results"), header = T, sep = "\t")
  cond <- row.names(meta2[meta2$filename == file,])
  read.file2 <- read.file[c("gene_id", "expected_count")] #change here if want to use TPM or other quantification method
  colnames(read.file2)[2] <- cond
  files.list[[cond]] <- read.file2
}
gene.counts <- Reduce(function(...) merge(..., by = "gene_id"), files.list)
row.names(gene.counts) <- gene.counts$gene_id
gene.counts <- gene.counts[,-1]

## Transform gene counts file into a matrix
setwd(basedir)
matrix.count <- as.matrix(gene.counts)
storage.mode(matrix.count) = "integer"

##Run DESEq2
cat("Start run with DESeq2")
dds <- DESeqDataSetFromMatrix(countData = matrix.count,
                              colData = meta2,
                              design = ~ condition)
keep <- rowMeans(counts(dds)) >= min.reads
dds <- dds[keep,]
dds3 <-DESeq(dds)
resultsNames(dds3) ##design = ~condition

## Extract tables from contrast analyses
log2FC <- data.frame(row.names(dds3))
names <- vector()
for(i in doses){
  for(j in times){
    comb.c <- paste0(control.name, "_", j)
    comb.t <- paste0(i, "_", j)
    names <- c(names, comb.t)
    res <- results(dds3, contrast = c("condition", paste0(comb.t), paste0(comb.c)), 
                   pAdjustMethod = p.adj.method)
    res.full <- as.data.frame(res)
    if(save.files){
    write.csv(res.full, paste0(basedir, comb.t, ".csv"), row.names = T, quote = F)
    }
    #print(paste0("Analysis of ", comb.t, " done, here's a summary:"))
    summary(res, alpha = 0.05)
    log2FC <- cbind(log2FC, res.full[, "log2FoldChange", drop = T])
  }
}

##Do some plots on the whole dataset
message("Writing output files")
#Write normalized table
if(norm.method == "vsd"){
  norm.table <- vst(dds3, blind=FALSE)
} else if(norm.method == "rlog"){
  norm.table <- rlog(dds3, blind=FALSE)
}
write.csv(assay(norm.table), paste0("norm_table_", norm.method, ".csv"), row.names = T, quote = F)

## Create summary plots
pdf("summary.pdf", onefile = T)
if(clustering){
  sampleDists <- dist(t(assay(norm.table)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(norm.table$dose, norm.table$time, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  
}
### ATTENTION: by default, DESeq filters some genes with low variance prior to PCA.
### keep that in mind when performing PCA on the normalized set
if(pca){
pca.plot <- plotPCA(norm.table, intgroup=c("time", "dose"),
                    returnData = T)
ggplot(pca.plot, aes(PC1, PC2)) + 
  geom_point(aes(shape = factor(time), color = factor(dose), size = 10)) + 
  ggtitle("PCA normalized data") +
  xlab("PC1") + ylab("PC2") +
  labs(color = "Dose", shape = "Time") +
  guides(size = FALSE)
}
dev.off()

###### the end ######
