
run_analyse <- function(folder, metadata, name, analyse_pairs=T, analyse_sex=T) {
    # limma
    x<-load_DGE(metadata,  "./htseq/") 
    output_file<-paste0(folder, "limma", name,".csv")
    run<-run_limma(x, metadata, output_file)
    
    length(which(run$adj.P.Val < 0.05))
    tmp<-run[run$adj.P.Val < 0.05,]
    tmp<-tmp[c(tmp$logFC < -1 | tmp$logFC >1),]
    nrow(tmp)
    
    #deseq2

    
    # deseq2 with genes
    if (analyse_pairs==T & analyse_sex==T){
        
    } else if (analyse_sex==T){

    run_deseq<-run_deseq2(x, metadata)
    run_deseq <- run_deseq[!is.na(run_deseq$padj),]

    } else if (analyse_pairs==T){

    # deseq2 with pairs
    metadata <- metadata %>%
      group_by(patient_id) %>%
      filter(n_distinct(type) == 2)
    x<-load_DGE(metadata,  "./htseq2/") 
    run_deseq<-run_deseq2_paired(x, metadata)
    }

    
    lncRNA<-result_sign[result_sign$gene_biotype=="lncRNA",]
    nrow(lncRNA)
    nrow(lncRNA[lncRNA$log2FoldChange< -1,])
    nrow(lncRNA[lncRNA$log2FoldChange>1,])
    
    run_deseq <- run_deseq[!is.na(run_deseq$padj),]
    length(which(run_deseq$padj < 0.05))
    output_file<-paste0(folder, "deseq_", name,".csv")
    write.csv2(run_deseq, output_file)
    tmp<-run_deseq[run_deseq$padj < 0.05,]
    tmp<-tmp[c(tmp$log2FoldChange < -1 | tmp$log2FoldChange >1),]
    nrow(tmp)

        
    result<-add_genes_info(run_deseq, ah)
    result[result$padj<0.05,] -> result_sign
    
    result_sign2 <- data.frame(lapply(result_sign, function(x) {
      if (is.list(x)) sapply(x, paste, collapse = ",") else x
    }))
    output_file<-paste0(folder, "deseq_genes_info_", name,".csv")
    write.table(result_sign2, output_file, sep = ";", row.names = FALSE)

    
    #plot PCA
    output_file<-paste0(folder, "plotPCA_", name,".svg")
    plotPCA(x, metadata, output_file)
    
    # enrichment GO
    enrich<-enrichGO(result_sign)
    enrich<- as.data.frame(enrich)
    output_file<-paste0(folder, "enrichmentGO_padj_0_05_", name,".csv")
    
    write.csv2(enrich, output_file)
    result_sign[c(result_sign$log2FoldChange<-1 | result_sign$log2FoldChange>1),] -> result_sign
    enrich<-enrichGO(result_sign)
    enrich_df<- as.data.frame(enrich)
    output_file<-paste0(folder, "enrichmentGO_padj_0_05_lfc_1_", name,".csv")
    
    write.csv2(enrich_df, output_file)
    
    
    pdf(paste0(folder,"emmaplot_", name, "_goEnrich.pdf"), width = 7, height = 7)
    ego <- pairwise_termsim(enrich)
    emapplot(ego)
    dev.off()

    original_gene_list <- result_sign$log2FoldChange
    names(original_gene_list) <- result_sign$Row.names
    enrich_tmp<-setReadable(enrich, 'org.Hs.eg.db', 'ENSEMBL')
    cnet<-cnetplot(enrich_tmp, foldChange=original_gene_list, showCategory=5)
      ggsave(
          paste0(folder,"cnet_", name, "_goEnrich.pdf"),
        plot = cnet,
      )


    # enrichment KEGG
    
    kegg<-runenrichKegg(result_sign, result_sign$gene_id, name)
    
    pdf(paste0(folder,"keggEnrich_",name,".pdf"), width = 7, height = 7)
    dotplot(kegg, showCategory=30, label_format=NULL) + ggtitle("dotplot for kegg enrichment")
    dev.off()
    
    
    kegg_tmp<-setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')
    cnet<-cnetplot(kegg_tmp, foldChange=original_gene_list, showCategory=5)
      ggsave(
          paste0(folder,"cnet_", name, "_KEGGEnrich.pdf"),
        plot = cnet,
      )

}
