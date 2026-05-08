
run_analyse <- function(folder, metadata_0_analyse, name, analyse_pairs=T, analyse_sex=T) {
    # limma
    x<-load_DGE(metadata_0_analyse,  "./htseq/") 
    output_file<-paste0(folder, "limma", name,".csv")
    run_0<-run_limma(x, metadata_0_analyse, output_file)
    
    length(which(run_0$adj.P.Val < 0.05))
    tmp<-run_0[run_0$adj.P.Val < 0.05,]
    tmp<-tmp[c(tmp$logFC < -1 | tmp$logFC >1),]
    nrow(tmp)
    
    #deseq2
    run_0_deseq<-run_deseq2(x, metadata_0_analyse)
    run_0_deseq <- run_0_deseq[!is.na(run_0_deseq$padj),]
    length(which(run_0_deseq$padj < 0.05))
    
    output_file<-paste0(folder, "deseq_", name,".csv")
    write.csv2(run_0_deseq, output_file)
    tmp<-run_0_deseq[run_0_deseq$padj < 0.05,]
    tmp<-tmp[c(tmp$log2FoldChange < -1 | tmp$log2FoldChange >1),]
    nrow(tmp)
    
    # deseq2 with genes
    
    result<-add_genes_info(run_0_deseq,ah)
    result[result$padj<0.05,] -> result_sign
    
    result_sign2 <- data.frame(lapply(result_sign, function(x) {
      if (is.list(x)) sapply(x, paste, collapse = ",") else x
    }))
    output_file<-paste0(folder, "deseq_genes_info_", name,".csv")
    write.table(result_sign2, output_file, sep = ";", row.names = FALSE)
    
    
    
    lncRNA<-result_sign[result_sign$gene_biotype=="lncRNA",]
    nrow(lncRNA)
    nrow(lncRNA[lncRNA$log2FoldChange< -1,])
    nrow(lncRNA[lncRNA$log2FoldChange>1,])
    
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
    
    
    
    pdf(paste0(folder,"emmaplot_", name, "_goEnrich_0.pdf"), width = 7, height = 7)
    ego <- pairwise_termsim(enrich)
    emapplot(ego)
    dev.off()

    original_gene_list <- result_sign$log2FoldChange
    names(original_gene_list) <- result_sign$Row.names
    enrich_tmp<-setReadable(enrich, 'org.Hs.eg.db', 'ENSEMBL')
    cnet<-cnetplot(enrich_tmp, foldChange=original_gene_list, showCategory=5)
      ggsave(
          paste0(folder,"cnet_", name, "_goEnrich_0.pdf"),
        plot = cnet,
      )
    
}
