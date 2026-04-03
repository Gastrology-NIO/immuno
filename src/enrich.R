library(clusterProfiler)
library(org.Hs.eg.db)
enrichGO<-function( genes,aspect="BP"){
  # MF
  # BP
  # CC
ego2 <- clusterProfiler::enrichGO(
  gene         = genes,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENSEMBL",
  ont          = aspect,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)
    return(ego2)
  }
