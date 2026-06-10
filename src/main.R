
run_analyse <- function(folder, metadata, name, analyse_pairs=T, analyse_sex=T, max_genes=0) {
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

    

    
    run_deseq <- run_deseq[!is.na(run_deseq$padj),]
    length(which(run_deseq$padj < 0.05))
    output_file<-paste0(folder, "deseq_", name,".csv")
    write.csv2(run_deseq, output_file)
    tmp<-run_deseq[run_deseq$padj < 0.05,]
    tmp<-tmp[c(tmp$log2FoldChange < -1 | tmp$log2FoldChange >1),]
    nrow(tmp)

    result<-add_genes_info(run_deseq, ah)
    result[result$padj<0.05,] -> result_sign
    lncRNA<-result_sign[result_sign$gene_biotype=="lncRNA",]
    nrow(lncRNA)
    nrow(lncRNA[lncRNA$log2FoldChange< -1,])
    nrow(lncRNA[lncRNA$log2FoldChange>1,])
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
    enrich_BP_df<- as.data.frame(enrich_BP)
    output_file<-paste0(folder, "enrichmentGO_BP_padj_0_05_lfc_1_", name,".csv")
    write.csv2(enrich_BP_df, output_file)

    enrich_MF<-enrichGO(result_sign, "MF")
    enrich_MF_df<- as.data.frame(enrich_MF)
    output_file<-paste0(folder, "enrichmentGO_MF_padj_0_05_lfc_1_", name,".csv")
    write.csv2(enrich_MF_df, output_file)

    enrich_CC<-enrichGO(result_sign, "CC")
    enrich_CC_df<- as.data.frame(enrich_CC)
    output_file<-paste0(folder, "enrichmentGO_CC_padj_0_05_lfc_1_", name,".csv")
    write.csv2(enrich_CC_df, output_file)

    library(dplyr)
    if (max_genes>0){
    # wybór 5 kategorii, które będą pokazane
    
    # dla każdej kategorii zostaw top20 genów wg |logFC|
    enrich_CC@result$geneID <- sapply(
      enrich_CC@result$geneID,
      function(x) {
        genes <- unlist(strsplit(x, "/"))
    
        genes_fc <- original_gene_list[genes]
    
        top_genes <- names(sort(abs(genes_fc),
                                decreasing = TRUE))[1:min(20, length(genes_fc))]
    
        paste(top_genes, collapse = "/")
      }
    )
    
    # dla każdej kategorii zostaw top20 genów wg |logFC|
    enrich_MF@result$geneID <- sapply(
      enrich_MF@result$geneID,
      function(x) {
        genes <- unlist(strsplit(x, "/"))
    
        genes_fc <- original_gene_list[genes]
    
        top_genes <- names(sort(abs(genes_fc),
                                decreasing = TRUE))[1:min(20, length(genes_fc))]
    
        paste(top_genes, collapse = "/")
      }
    )
    
    # dla każdej kategorii zostaw top20 genów wg |logFC|
    enrich_BP@result$geneID <- sapply(
      enrich_BP@result$geneID,
      function(x) {
        genes <- unlist(strsplit(x, "/"))
    
        genes_fc <- original_gene_list[genes]
    
        top_genes <- names(sort(abs(genes_fc),
                                decreasing = TRUE))[1:min(20, length(genes_fc))]
    
        paste(top_genes, collapse = "/")
      }
    )
    }
    
    ego <- pairwise_termsim(enrich_BP)
    e1<-emapplot(ego)

    ego <- pairwise_termsim(enrich_MF)
    e2<-emapplot(ego)
    if (nrow(enrich_CC_df)> 1){

    ego <- pairwise_termsim(enrich_CC)
    e3<-emapplot(ego)
    }
    pdf(paste0(folder,"emmaplot_BH_", name, "_goEnrich.pdf"), width = 7, height = 7)
    print(e1)
    dev.off()
        if (nrow(enrich_CC_df)> 1){

    pdf(paste0(folder,"emmaplot_CC_", name, "_goEnrich.pdf"), width = 7, height = 7)
    print(e3)
    dev.off()
}
    pdf(paste0(folder,"emmaplot_MF_", name, "_goEnrich.pdf"), width = 7, height = 7)
    print(e2)
    dev.off()
    
    pdf(paste0(folder,"dotplot_GOEnrich_BP_",name,".pdf"), width = 7, height = 7)
       d_bp<- dotplot(enrich_BP, showCategory=30, label_format=NULL) + ggtitle("Biological Process enrichment")
    print(d_bp)
    dev.off()

    pdf(paste0(folder,"dotplot_GOEnrich_MF_",name,".pdf"), width = 7, height = 7)
       d_mf<- dotplot(enrich_MF, showCategory=30, label_format=NULL) + ggtitle("Molecular Function enrichment")
    print(d_mf)
    dev.off()
    if (nrow(enrich_CC_df)> 1){
        pdf(paste0(folder,"dotplot_GOEnrich_CC_",name,".pdf"), width = 7, height = 7)
       d_CC<- dotplot(enrich_CC, showCategory=30, label_format=NULL) + ggtitle("Cellular Component enrichment")
        print(d_CC)
        dev.off()
    }
    
    original_gene_list <- result_sign$log2FoldChange
    names(original_gene_list) <- result_sign$Row.names
    if (nrow(enrich_BP)>1){
            enrich_tmp<-setReadable(enrich_BP, 'org.Hs.eg.db', 'ENSEMBL')
            cnet_BP<-cnetplot(enrich_tmp, foldChange=original_gene_list, showCategory=5)
            ggsave(
                  paste0(folder,"cnet_", name, "_goEnrich.pdf"),
                plot = cnet_BP,
              )
    }
    
    if (nrow(enrich_MF)>1){
        enrich_tmp<-setReadable(enrich_MF, 'org.Hs.eg.db', 'ENSEMBL')

        cnet_MF<-cnetplot(enrich_tmp, foldChange=original_gene_list, showCategory=5)
        ggsave(
              paste0(folder,"cnet_", name, "_MF_goEnrich.pdf"),
            plot = cnet_MF,
          )
    }
    
    if (nrow(enrich_CC)>1){
        enrich_tmp<-setReadable(enrich_CC, 'org.Hs.eg.db', 'ENSEMBL')
    cnet_CC<-cnetplot(enrich_tmp, foldChange=original_gene_list, showCategory=5)
    ggsave(
          paste0(folder,"cnet_", name, "_cc_goEnrich.pdf"),
        plot = cnet_CC,
      )
        }
    # enrichment KEGG
    
    kegg<-runenrichKegg(result_sign, result_sign$gene_id, name)
    df<-as.data.frame(kegg)
        output_file<-paste0(folder, "enrichmentKEGG_padj_0_05_", name,".csv")
    write.csv2(df, output_file)
    
    pdf(paste0(folder,"dotplot_keggEnrich_",name,".pdf"), width = 7, height = 7)
       d_kegg<- dotplot(kegg, showCategory=30, label_format=NULL) + ggtitle("KEGG enrichment")
    print(d_kegg)
    dev.off()
    
    
    kegg_tmp<-setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')
    
    original_gene_list <- result_sign$log2FoldChange
    names(original_gene_list) <- result_sign$entrezid
    cnet_kegg<-cnetplot(kegg_tmp, foldChange=original_gene_list, showCategory=5)
    ggsave(
          paste0(folder,"cnet_", name, "_KEGGEnrich.pdf"),
        plot = cnet_kegg
      )
    
      plots <- list()
    id=1
    if (nrow(enrich_BP)>1){
      plots[[id]]<-d_bp
        id<-id+1
    }
    if (nrow(enrich_MF)>1){
          plots[[id]]<-d_mf
            id<-id+1
    }
    if (nrow(enrich_CC)>1){
        plots[[id]]<-d_CC
        id<-id+1
    }
      plots[[id]]<-d_kegg
        p<-wrap_plots(plots, ncol = 2) +
          plot_annotation(tag_levels = "A")
       ggsave(
             paste0("enrichgokegg_genes_dotplot_",name,".svg"),
            plot = p,
              width = 18, height = 12,
        )
    print(paste0("enrichgokegg_genes_dotplot_",name,".svg"))
           
    plots <-list()
    id=1
    if (nrow(enrich_BP)>1){
        plots[[id]]<-cnet_BP
        id<-id+1
    }
    if (nrow(enrich_MF)>1){
      plots[[id]]<-cnet_MF
    id<-id+1
    }
    if (nrow(enrich_CC)>1){
          plots[[id]]<-cnet_CC
        id<-id+1

        }
      plots[[id]]<-cnet_kegg
        p<-wrap_plots(plots, ncol = 2) +
          plot_annotation(tag_levels = "A")
       ggsave(
             paste0("enrichgokegg_genes_Cnet_",name,".svg"),
            plot = p,
              width = 12, height = 12,
        )
    print(paste0("enrichgokegg_genes_Cnet_",name,".svg"))
               
}
