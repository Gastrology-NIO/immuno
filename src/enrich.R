library(clusterProfiler)
library(org.Hs.eg.db)
enrichGO<-function( data,aspect="BP"){
  geneList <- data$Row.names
  #names(geneList) <- data$Row.names

  # MF
  # BP
  # CC
ego <- clusterProfiler::enrichGO(
  gene          = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = aspect,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

    return(ego)
  }
