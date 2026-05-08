
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

        deseq_res<-run_deseq2(x, metadata)
        run_deseq<-deseq_res$res
        dds<-deseq_res$dds
        model <- dist_mat ~ s + research

    } else if (analyse_pairs==T){

        # deseq2 with pairs
        metadata <- metadata %>%
          group_by(patient_id) %>%
          filter(n_distinct(type) == 2)
        x<-load_DGE(metadata,  "./htseq2/") 
        deseq_res<-run_deseq2_paired(x, metadata)
        run_deseq<-deseq_res$res
        dds<-deseq_res$dds
                model <- dist_mat ~ patient_id + research

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
    plotPCA(dds, metadata, output_file, model)
    
    # enrichment GO
    enrich_BP<-enrichGO(result_sign, "BP")
    enrich_BP<- as.data.frame(enrich_BP)
    output_file<-paste0(folder, "enrichmentGO_BH_padj_0_05_", name,".csv")
    write.csv2(enrich_BP, output_file)
    
    result_sign[c(result_sign$log2FoldChange<-1 | result_sign$log2FoldChange>1),] -> result_sign
    enrich_BP<-enrichGO(result_sign, "BP")
    enrich_BP_df<- as.data.frame(enrich)
    output_file<-paste0(folder, "enrichmentGO_BP_padj_0_05_lfc_1_", name,".csv")
    write.csv2(enrich_BP_df, output_file)

        enrich_MF<-enrichGO(result_sign, "MF")
    enrich_MF_df<- as.data.frame(enrich_MF)
    output_file<-paste0(folder, "enrichmentGO_MF_padj_0_05_lfc_1_", name,".csv")
    write.csv2(enrich_MF_df, output_file)

        enrich_CC<-enrichGO(result_sign, "CC")
    enrich_CC_df<- as.data.frame(enrich)
    output_file<-paste0(folder, "enrichmentGO_CC_padj_0_05_lfc_1_", name,".csv")
    write.csv2(enrich_CC_df, output_file)
    
    pdf(paste0(folder,"emmaplot_", name, "_goEnrich.pdf"), width = 7, height = 7)
    ego <- pairwise_termsim(enrich_BP)
    e1<-emapplot(ego)

    ego <- pairwise_termsim(enrich_MF)
    e2<-emapplot(ego)

    ego <- pairwise_termsim(enrich_CC)
    e3<-emapplot(ego)

    e<-patchwork::wrap_plots(
        e1,e2,e3)
    print(e)
    dev.off()

    original_gene_list <- result_sign$log2FoldChange
    names(original_gene_list) <- result_sign$Row.names
    enrich_tmp<-setReadable(enrich_BP, 'org.Hs.eg.db', 'ENSEMBL')
    cnet<-cnetplot(enrich_tmp, foldChange=original_gene_list, showCategory=5)
    ggsave(
          paste0(folder,"cnet_", name, "_goEnrich.pdf"),
        plot = cnet,
      )


    # enrichment KEGG
    
    kegg<-runenrichKegg(result_sign, result_sign$gene_id, name)
    df<-as.data.frame(kegg)
        output_file<-paste0(folder, "enrichmentKEGG_padj_0_05_", name,".csv")
    write.csv2(df, output_file)
    
    pdf(paste0(folder,"dotplot_keggEnrich_",name,".pdf"), width = 7, height = 7)
       d<- dotplot(kegg, showCategory=30, label_format=NULL) + ggtitle("dotplot for kegg enrichment")
    print(d)
    dev.off()
    
    
    kegg_tmp<-setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')
    
    original_gene_list <- result_sign$log2FoldChange
    names(original_gene_list) <- result_sign$entrezid
    cnet<-cnetplot(kegg_tmp, foldChange=original_gene_list, showCategory=5)
      ggsave(
          paste0(folder,"cnet_", name, "_KEGGEnrich.pdf"),
        plot = cnet,
      )

}
