library(clusterProfiler)
library(org.Hs.eg.db)
enrichGO<-function( data,aspect="BP"){
  geneList <- data$logFC
  names(geneList) <- data$Row.names

  # MF
  # BP
  # CC
ego <- clusterProfiler::enrichGO(
  gene          = genes,
  universe      = names(geneList),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = aspect,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
  ego2 <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
ego2 <- clusterProfiler::pairwise_termsim(ego2)
ego2@result$logFC <- geneList[ego2@result$geneID]
    return(ego2)
  }
