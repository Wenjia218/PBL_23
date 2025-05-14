library( "DESeq2" )
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggpubr)
library("org.Sc.sgd.db")

# Introduce( what do we do, why is it important) - methods( DESeq ) - results
# Run DESeq( input, output, scripts, graphics)
# ATAC foldchanges, find matches --> 100 base - 500 base, bereich offen (more active) oder geschlossen( less TFs binds on the promoter)
# for each gene several foldchanges
# DESeq for ATAC and DESeq for RNA (heatshock) Seq of foldchanges, and that they are combined with each other
# one looks into one gene along and the other looks into the co-regulation

# verarbeitet werden -> normal bei 30 degree 30 minutes -> not very hot

heatshock <- read.delim("complex_yeast_heatshock.tsv")
df_ATAC <- read.delim("ATACcounts_promotor_us500_ds100.tsv")
TFLink <- read.delim('TFLink_Saccharomyces_cerevisiae_interactions_SS_simpleFormat_v1.0.tsv')


# prepare data from heat shock
names <- grep("Wildtype", names(heatshock), value = TRUE)
df_WT <- heatshock[, names]
rownames(df_WT) <- heatshock$gene_id
# df_ATAC and df_WT will be used for following process of DEseq2

# filter only genes that both in ATAC data and RNA data are
gene_ids <- intersect(rownames(df_WT), rownames(df_ATAC))
ATAC_data <- df_ATAC[, c(1:4, 9:16, 21:28)]
RNA_data <- df_WT[gene_ids,]
order <- order(colnames(RNA_data))
RNA_data <- RNA_data[,order]

## DESeq 

# condition: 1 = 25/30_10, 2 = 25/30_30, 3 = 25/42_10, 4 = 25/42_30

calculateDESeq <- function(ATAC_data, RNA_data, condition){
  # fold changes of ATAC_data
  compare1 <- ATAC_data[, c(1:4, condition*4+1, condition*4+2, condition*4+3, condition*4+4)]
  coldata <- data.frame(condition = factor(c('control', 'control', 'control', 'control', 'treat', 'treat', 'treat','treat'), levels = c('control', 'treat')))
  dds <- DESeqDataSetFromMatrix(countData = compare1, colData = coldata, design= ~condition)
  dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
  res <- results(dds1, contrast = c('condition', 'treat', 'control'))

  # fold changes of RNA_data
  compare2 <- RNA_data[, c(1:3, condition*3+1, condition*3+2, condition*3+3)]
  coldata2 <- data.frame(condition = factor(c('control','control','control','treat', 'treat','treat'), levels = c('control', 'treat')))
  dds2 <- DESeqDataSetFromMatrix(countData = compare2, colData = coldata2, design= ~condition)
  dds3 <- DESeq(dds2, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
  res2 <- results(dds3, contrast = c('condition', 'treat', 'control'))

  # output of the DESeq data in files
  res_df <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  write.table(res_df, 'ATACData.DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)
  res_df2 <- data.frame(res2, stringsAsFactors = FALSE, check.names = FALSE)
  write.table(res_df2, 'RNAData.DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)
  res_ATAC <- res_df
  res_ATAC$gene_id <- rownames(res_ATAC)
  res_RNA <- res_df2
  res_RNA$gene_id <- rownames(res_RNA)

  join_ATAC_RNA <- inner_join(res_ATAC, res_RNA, by= "gene_id")
  return(join_ATAC_RNA)
}

makeplot <- function(ATAC_data, RNA_data, condition, highlights, TF){
    join_ATAC_RNA <- calculateDESeq(ATAC_data, RNA_data, condition)
    highlights <- join_ATAC_RNA %>% filter(gene_id %in% highlights)
    title <- switch(condition, "25/30_10", "25/30_30", "25/42_10", "25/42_30")
    title <- paste("foldchanges of the same gene ", title, sep=" ")
    g <- ggplot(join_ATAC_RNA, aes(x = log2FoldChange.x, y = log2FoldChange.y)) +
      geom_point() +
      geom_point(data = highlights, aes(color = TF)) + 
      geom_hline(yintercept=0, linetype="dashed", color = "blue") +
      geom_vline(xintercept=0, linetype="dashed", color = "blue") +
      labs(title = title) +
      ylab("RNA") + 
      xlab("ATAC") 
    return(g)
}

# tried to show four plots of four different condition in one figure
mergedplot <- function(ATAC_data, RNA_data, highlights, TF){
  g1 <- makeplot(ATAC_data, RNA_data, 1, highlights, TF)
  g2 <- makeplot(ATAC_data, RNA_data, 2, highlights, TF)
  g3 <- makeplot(ATAC_data, RNA_data, 3, highlights, TF)
  g4 <- makeplot(ATAC_data, RNA_data, 4, highlights, TF)
  figure <- ggarrange(g1,g2,g3,g4,
                    labels = c("25-30_10", "25-30_30", "25-42_10", "25-42_30"),
                    ncol = 2, nrow = 2)
  figure
}


TFs = c("YGL073W","YMR037C", "YKL062W" ) # HSF1, MSN2 and MSN4


# find targets from TFLink Database
TFtargets_HSF1 <- TFLink %>% filter(`Name.TF` == "HSF1")
TF_HSF1<-mapIds(org.Sc.sgd.db, keys=TFtargets_HSF1$UniprotID.Target, column="ORF", keytype="UNIPROT", multiVals="first")
TFtargets_MSN2 <- TFLink %>% filter(`Name.TF` == "MSN2")
TF_MSN2<-mapIds(org.Sc.sgd.db, keys=TFtargets_MSN2$UniprotID.Target, column="ORF", keytype="UNIPROT", multiVals="first")
TFtargets_MSN4 <- TFLink %>% filter(`Name.TF` == "MSN4")
TF_MSN4<-mapIds(org.Sc.sgd.db, keys=TFtargets_MSN4$UniprotID.Target, column="ORF", keytype="UNIPROT", multiVals="first")
TF_MSN24 <- intersect(TF_MSN2, TF_MSN4)


mergedplot(ATAC_data, RNA_data, TF_HSF1, "targets of HSF1")
mergedplot(ATAC_data, RNA_data, TF_MSN24, "targets of MSN24")


# Violin plots
heatshockDESeq <- function(RNA_data, condition1){
 # fold changes of RNA_data
  compare2 <- RNA_data[, c(1:3, condition*3+1, condition*3+2, condition*3+3)]
  coldata2 <- data.frame(condition = factor(c('control','control','control','treat', 'treat','treat'), levels = c('control', 'treat')))
  dds2 <- DESeqDataSetFromMatrix(countData = compare2, colData = coldata2, design= ~condition)
  dds3 <- DESeq(dds2, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
  res2 <- results(dds3, contrast = c('condition', 'treat', 'control'))
  # output of the DESeq data in files
  res_df2 <- data.frame(res2, stringsAsFactors = FALSE, check.names = FALSE)
  write.table(res_df2, 'RNAData.DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)
  res_RNA <- res_df2
  res_RNA$gene_id <- rownames(res_RNA)
  return(res_RNA)
}

# condition: Wildtype, MSN24.KO, HSF1.KD
violinplot <- function(data, condition){
  names <- grep(condition, names(data), value = TRUE)
  df <- heatshock[, names]
  rownames(df) <- heatshock$gene_id
  RNA_data <- df
  order <- order(colnames(RNA_data))
  RNA_data <- RNA_data[,order]
  
  foldchanges1 <- heatshockDESeq(RNA_data, 1)
  foldchanges2 <- heatshockDESeq(RNA_data, 2)
  foldchanges3 <- heatshockDESeq(RNA_data, 3)
  foldchanges4 <- heatshockDESeq(RNA_data, 4)
  newdata <- data.frame(gene_ids = foldchanges1$gene_id,
                        X37_10 = foldchanges1$log2FoldChange,
                        X37_30 = foldchanges2$log2FoldChange,
                        X42_10 = foldchanges3$log2FoldChange,
                        X42_30 = foldchanges4$log2FoldChange)
  longdata <- newdata %>%
    pivot_longer(
      cols = `X37_10`:`X42_30`, 
      names_to = "condition",
      values_to = "foldchanges",
   )
  longdata$foldchanges <- as.numeric(longdata$foldchanges)
  title <- paste("foldchanges between 25° and condition in ", condition)
  g <- ggplot(longdata, aes(x = condition, y = foldchanges, fill = condition)) +
    labs(title = title) +
    geom_violin() +
    ylab("foldchanges")
  return(g)
}

g1 <- violinplot(heatshock,"Wildtype")
g2 <- violinplot(heatshock,"MSN24.KO")
g3 <- violinplot(heatshock,"HSF1.KD")

figure <- ggarrange(g1,g2,g3,
                    ncol = 2, nrow = 2)
figure


# genes with the biggest foldchanges

DESeqRNA <- heatshockDESeq(RNA_data, 2)
data <- filter(DESeqRNA, row_number(desc(log2FoldChange)) <= 10)
genes <- rownames(data)
# YBR072W: Small heat shock protein (sHSP) with chaperone activity
# YDR171W: Small heat shock protein (sHSP) with chaperone activity, forms barrel-shaped oligomers that suppress unfolded protein aggregation
# YDR258C: Oligomeric mitochondrial matrix chaperone
# YER103W: Heat shock protein that is highly induced upon stress
# YFL014W: Plasma membrane protein involved in maintaining membrane organization
# YGR142W: v-SNARE binding protein
# LL026W: heat shock protein that cooperates with Ydj1p (Hsp40) and Ssa1p (Hsp70) to refold and reactivate denatured, protein aggregates



