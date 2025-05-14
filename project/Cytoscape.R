
library(BiocManager)
library(RCy3)
library(igraph)
library(dplyr)
library(tidyverse)
library("org.Sc.sgd.db")


#TF Link
#get dataframe from TFLink
TFLink <- read.delim('TFLink_Saccharomyces_cerevisiae_interactions_SS_simpleFormat_v1.0.tsv')
TFLink$Name.Target<-mapIds(org.Sc.sgd.db, keys=TFLink$UniprotID.Target, column="ORF", keytype="UNIPROT", multiVals="first")


TF_HSF1 <-
  TFLink %>% filter(`Name.TF` == "HSF1")
TF_HSF1 <- TF_HSF1[,5:6]
TF_MSN2 <-
  TFLink %>% filter(`Name.TF` == "MSN2")
TF_MSN2 <- TF_MSN2[,5:6]
TF_MSN4 <-
  TFLink %>% filter(`Name.TF` == "MSN4")
TF_MSN4 <- TF_MSN4[,5:6]

#prepare data for the network
nodesHSF1 <- levels(factor(as.vector(TF_HSF1[,2])))
nodesHSF1 <- append(nodesHSF1, "HSF1")

nodesMSN2 <- levels(factor(as.vector(TF_MSN2[,2])))
nodesMSN2 <- append(nodesMSN2, "MSN2")

nodesMSN4 <- levels(factor(as.vector(TF_MSN4[,2])))
nodesMSN4 <- append(nodesMSN4, "MSN4")

TF_HSF1_MSN24 <- rbind(TF_HSF1,TF_MSN2, TF_MSN4)
nodes <- c(nodesHSF1, nodesMSN2, nodesMSN4)
nodes <- unique(nodes)
ig1 <- graph_from_data_frame(TF_HSF1_MSN24, directed=FALSE, vertices=nodes)
createNetworkFromIgraph(ig1,"myIgraphHSF1_MSN24")



