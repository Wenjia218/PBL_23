condition <- c("X25_0","X25_0","X25_0","X37_10","X37_10","X37_10","X37_30","X37_30","X37_30","X42_10","X42_10","X42_10","X42_30","X42_30","X42_30")
condition2 <- c("X25","X25","X25","X25","X37_5","X37_5","X37_5","X37_5","X37_10","X37_10","X37_10","X37_10","X37_30","X37_30","X37_30","X37_30","X42_5","X42_5","X42_5","X42_5","X42_10","X42_10","X42_10","X42_10","X42_30","X42_30","X42_30","X42_30")


mergeData <- function(data){
  data <- rbind(data, colnames(data))
  data <- data.frame(t(data), stringsAsFactors = F)
  colnames(data) <- c( "values","conditions")
  rownames(data) <- c(1:15)
  data$values <- as.numeric(data$values)
  return(data)
}

mergeData2 <- function(data){
  data <- rbind(data, colnames(data))
  data <- data[,-29]
  data <- data.frame(t(data), stringsAsFactors = F)
  colnames(data) <- c( "values","conditions")
  rownames(data) <- c(1:28)
  data$values <- as.numeric(data$values)
  return(data)
}

violin_plots <- function(gene){
  df_gene <- 
    df %>% 
    filter(`gene_id` == gene)
  names <- grep("HSF1", names(df_gene), value = TRUE)
  HSF1_KD <- df_gene[, names]
  HSF1_KD <- mergeData(HSF1_KD)
  
  names <- grep("MSN24", names(df_gene), value = TRUE)
  MSN24_KO <- df_gene[, names]
  MSN24_KO <- mergeData(MSN24_KO)
  
  names <- grep("Wildtype", names(df_gene), value = TRUE)
  WT <- df_gene[, names]
  WT<- mergeData(WT)
  
  myTitle1 = paste("Overview HSF1 KD", gene, sep=" ")
  g1 <- ggplot(HSF1_KD, aes(x = condition, y = values, fill = condition)) +
    geom_violin() +
    labs(title = myTitle1) +
    ylab("Value") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  myTitle2 = paste("Overview MSN24 KO", gene, sep=" ")
  g2 <- ggplot(MSN24_KO, aes(x = condition, y = values, fill = condition)) +
    geom_violin() +
    labs(title = myTitle2) +
    ylab("Value") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  myTitle3 = paste("Overview Wildtype", gene, sep=" ")
  g3 <- ggplot(WT, aes(x = condition, y = values, fill = condition)) +
    geom_violin() +
    labs(title = myTitle3) +
    ylab("Value") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  
  # add a plot of counts binding sites
  df_ATAC$gene_id = rownames(df_ATAC)
  df_gene2 <- 
    df_ATAC %>% 
    filter(`gene_id` == gene)
  df_gene2 <- mergeData2(df_gene2)
  
  myTitle4 = "Counts of binding sites"
  g4 <- ggplot(df_gene2, aes(x = condition2, y = values, fill = condition2)) +
    geom_violin() +
    labs(title = myTitle4) +
    ylab("Value") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  library(gridExtra)
  combined_plots <- grid.arrange(g1, g2, g3, g4, ncol = 2)
  
}

violin_plots("YGL073W")
violin_plots("YMR037C")
violin_plots("YKL062W")

