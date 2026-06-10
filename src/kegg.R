BiocManager::install(pathview)
library(pathview)

see_pathview <- function(..., save_image = FALSE)
{
  msg <- capture.output(pathview::pathview(...), type = "message")
  msg <- grep("image file", msg, value = T)
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  img <- png::readPNG(filename)
  grid::grid.raster(img)
  if(!save_image) invisible(file.remove(filename))
}
                     

rungseKegg<-function(res, original_gene_list){
  df<-as.data.frame(res)
  ids<-bitr(original_gene_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
  dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  df2 = df[df$X %in% dedup_ids$ENSEMBL,]
  
  # Create a new column in df2 with the corresponding ENTREZ IDs
  df2$Y = dedup_ids$ENTREZID
  
  # Create a vector of the gene unuiverse
  kegg_gene_list <- df2$log2FoldChange
  
  # Name vector with ENTREZ ids
  names(kegg_gene_list) <- df2$Y
  
  # omit any NA values 
  kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  kegg_gene_list <- kegg_gene_list[!duplicated(names(kegg_gene_list))]
  kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)
  kk2 <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = 'hsa',
                 nPerm        = 10000,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr",
                 keyType       = "ncbi-geneid")
  df<-as.data.frame(kk2)
  write.csv2(df, "gseKEGG.csv")
  return(kk2)
}
        

runenrichKegg<-function(res, original_gene_list, prefix){
  df<-as.data.frame(res)
  organism<-	org.Hs.eg.db
  ids<-bitr(original_gene_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
  dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  df2 = df[df$gene_id %in% dedup_ids$ENSEMBL,]
  
  # Create a new column in df2 with the corresponding ENTREZ IDs
  df2$Y = dedup_ids$ENTREZID
  
  # Create a vector of the gene unuiverse
  kegg_gene_list <- df2$log2FoldChange
  
  # Name vector with ENTREZ ids
  names(kegg_gene_list) <- df2$Y
  
  # omit any NA values 
  kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  kegg_gene_list <- kegg_gene_list[!duplicated(names(kegg_gene_list))]
  kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)
  kk2 <- clusterProfiler::enrichKEGG(names(kegg_gene_list),
                 organism     = 'hsa',
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr",
                  # keyType = "kegg"
                 keyType       = "ncbi-geneid")

  return(kk2)
}
plotDotPlot<-function(kk2, save_path){
  d<-dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
  ggsave(
    save_path,
    plot = d,
  )
}

plotKeggPath<-function(kk2,save_path="./hsa04814_kegg.svg", path_id="hsa04814"){

  # Produce the native KEGG plot (PNG)
  dme <- see_pathview(gene.data=kegg_gene_list, pathway.id=path_id, 
                      species = 'hsa', kegg.native = F)
  ggsave(
    save_path,
    plot = dme,
  )
  # Produce a different plot (PDF) (not displayed here)
  # dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04814", species = kegg_organism, kegg.native = F)

}


