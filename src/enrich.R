library(clusterProfiler)
library(org.Hs.eg.db)
enrichGO<-function( data,aspect="BP"){
  geneList <- data$logFC
  names(geneList) <- data$Row.names

  # MF
  # BP
  # CC
ego <- clusterProfiler::enrichGO(
  gene          = data$Row.names,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = aspect,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

    return(ego)
  }
