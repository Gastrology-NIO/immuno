
BiocManager::install('enrichplot', character.only = TRUE)
install.packages('ggtree')
library(clusterProfiler)
library(enrichplot)
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)



runGSEA<-function(differetial){
  original_gene_list <- top.table$logFC
  names(original_gene_list) <- rownames(top.table)
  gene_list<-na.omit(original_gene_list)
  c = sort(gene_list, decreasing = TRUE)
  gse <- gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "ENSEMBL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.1, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "fdr")
  gse[gse@result$p.adjust<0.1]->gse_tmp
  df<-as.data.frame(gse)
  write.csv2(df, "gseGO-01.csv")
  return(gse)  
}

plotDotPlot<-function(gse, save_path){
    require(DOSE)
    D<-dotplot(gse, showCategory=10, split=c(".sign")) + facet_grid(~.sign) +
      theme(
        axis.text.y = element_text(size = 8)  # Change y-axis label size
      )
    ggsave(
      save_path,
      plot = D,
    )
  }


plotEmaPlot<-function(gse, save_path="./emapplot.svg"){
  edo <- pairwise_termsim(gse)
  E<-enrichplot::emapplot(edo, showCategory = 20, layout = "kk")
  ggsave(
    save_path,
    plot = E,
  )
}


plotEmaPlot<-function(gse,gene_list, save_path="./cnetplot2.png"){
  cnet<-cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
  ggsave(
    save_path,
    plot = cnet,
    width=10, height=8, dpi=300
  )
}
