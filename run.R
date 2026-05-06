
add_genes_info <- function(data){
  gene_symbols<-data$Row.names
  ah <- AnnotationHub()
  query(ah, c("Homo sapiens", "EnsDb"))
  edb <- ah[["AH119325"]]
  genes <- genes(edb, filter = GeneIdFilter(gene_symbols))
  genes<-as.data.frame(genes)
  data<-merge(genes, data, by=0)
  return(data)
}




source("./src/load_data.R")


metadata<-load_data_meta("./data/metadata_jag.csv")
gene_counts<-load_data_counts("./htseq/", "./data/gene_counts.csv")




# todo: dodać nazwy genów itd

# 1 pobranie VS kontrola
conditions <- list(type=c('Pobranie 1', 'kontrola'))
metadata_0_analyse <- cut_metadata(metadata, conditions)
metadata_0_analyse$research<-metadata_0_analyse$type
metadata_0_analyse$research[metadata_0_analyse$type=='kontrola']<-"Control"
metadata_0_analyse$research[metadata_0_analyse$type=='Pobranie 1']<-"Searched"
metadata_0_analyse<- metadata_0_analyse[metadata_0_analyse$probe_name != "",]

x<-load_DGE(metadata_0_analyse,  "./htseq/") 
output_file<-"./result/pobranie1_vs_kontrola.csv"
run_0<-run_limma(x, metadata_0_analyse, output_file)
length(which(run_0$adj.P.Val < 0.05))
tmp<-run_0[run_0$adj.P.Val < 0.05,]
tmp<-tmp[c(tmp$logFC < -1 | tmp$logFC >1),]
nrow(tmp)


run_0_deseq<-run_deseq2(x, metadata_0_analyse)
run_0_deseq <- run_0_deseq[!is.na(run_0_deseq$padj),]
length(which(run_0_deseq$padj < 0.05))

write.csv2(run_0_deseq, "deseq_0.csv")
tmp<-run_0_deseq[run_0_deseq$padj < 0.05,]
tmp<-tmp[c(tmp$log2FoldChange < -1 | tmp$log2FoldChange >1),]
nrow(tmp)


result<-add_genes_info(run_0_deseq,ah)
result[result$padj<0.05,] -> result_sign
enrich<-enrichGO(result_sign)
enrich<- as.data.frame(enrich)
write.csv2(enrich, "./result/pobranie1_vs_kontrola_enrich_padj_0_05.csv")
result_sign[c(result_sign$log2FoldChange<-1 | result_sign$log2FoldChange>1),] -> result_sign
enrich_0<-enrichGO(result_sign)
enrich_0_df<- as.data.frame(enrich_0)
write.csv2(enrich_0_df, "./result/pobranie1_vs_kontrola_enrich_padj_0_05_logfc.csv")



pdf("goEnrich_0.pdf", width = 7, height = 7)
ego <- pairwise_termsim(enrich_0)
emapplot(ego)
dev.off()



           

# 1 pobranie (3 miesiące vs powyżej 2 lat)
conditions <- list(type='Pobranie 1', 'time.of.OS'=c('poni\xbfej 3 miesi\xeacy', '> 2 lata'))
metadata_1_analyse <- cut_metadata(metadata, conditions)
metadata_1_analyse$research<-metadata_1_analyse$'time.of.OS'
metadata_1_analyse$research[metadata_1_analyse$research=='poni\xbfej 3 miesi\xeacy']<-"Control"
metadata_1_analyse$research[metadata_1_analyse$research=='> 2 lata']<-"Searched"
x<-load_DGE(metadata_1_analyse,  "./htseq/") 
output_file<-"./result/1pobranie_3miesiace_vs_powyżej_2_lata.csv"
run_1<-run_limma(x, metadata_1_analyse, output_file)
length(which(run_1$adj.P.Val < 0.05))

tmp<-run_1[run_1$adj.P.Val < 0.05,]
tmp<-tmp[c(tmp$logFC < -1 | tmp$logFC >1),]
nrow(tmp)


run_1_deseq<-run_deseq2(x, metadata_1_analyse)
run_1_deseq <- run_1_deseq[!is.na(run_1_deseq$padj),]
length(which(run_1_deseq$padj < 0.05))


write.csv2(run_1_deseq, "deseq_1.csv")
tmp<-run_1_deseq[run_1_deseq$padj < 0.05,]
tmp<-tmp[c(tmp$log2FoldChange < -1 | tmp$log2FoldChange >1),]
nrow(tmp)


result<-add_genes_info(run_1_deseq, ah)
result[result$padj<0.05,] -> result_sign
enrich<-enrichGO(result_sign)
enrich<- as.data.frame(enrich)
write.csv2(enrich, "./result/1pobranie_3miesiace_vs_powyżej_2_lata_enrich_padj_0_05.csv")
result_sign[c(result_sign$log2FoldChange<-1 | result_sign$log2FoldChange>1),] -> result_sign
enrich_1<-enrichGO(result_sign)
enrich_1_df<- as.data.frame(enrich)
write.csv2(enrich, "./result//1pobranie_3miesiace_vs_powyżej_2_latae_enrich_padj_0_05_logfc.csv")


pdf("goEnrich_1.pdf", width = 7, height = 7)
ego <- pairwise_termsim(enrich_1)
emapplot(ego)
dev.off()



# 1 pobranie (I linia vs kolejna linia)
conditions <- list(type='Pobranie 1', 'linia'=c('I linia', ''))
metadata_2_analyse <- cut_metadata(metadata, conditions)
metadata_2_analyse$research<-metadata_2_analyse$'linia'
metadata_2_analyse$research[metadata_2_analyse$research=='I linia']<-"Control"
metadata_2_analyse$research[metadata_2_analyse$research=='']<-"Searched"
x<-load_DGE(metadata_2_analyse,  "./htseq/") 
output_file<-"./result/1pobranie_1_linia_vs_kolejne.csv"
run_2<-run_limma(x, metadata_2_analyse, output_file)
length(which(run_2$adj.P.Val < 0.05))
run_2_deseq<-run_deseq2(x, metadata_2_analyse)
run_2_deseq <- run_3_deseq[!is.na(run_2_deseq$padj),]
length(which(run_2_deseq$padj < 0.05))

write.csv2(run_2_deseq, "deseq_2.csv")
tmp<-run_2_deseq[run_2_deseq$padj < 0.05,]
tmp<-tmp[c(tmp$log2FoldChange < -1 | tmp$log2FoldChange >1),]
nrow(tmp)

# NtproBNP==0 (1 pobranie vs 2 pobranie)
conditions <- list(type=c('Pobranie 1', 'Pobranie 2'), 'wzrost.NtproBNP'=0)
metadata_3_analyse <- cut_metadata(metadata, conditions)
metadata_3_analyse$research<-metadata_3_analyse$type
metadata_3_analyse$research[metadata_3_analyse$research=='Pobranie 1']<-"Control"
metadata_3_analyse$research[metadata_3_analyse$research=='Pobranie 2']<-"Searched"
metadata_3_analyse<- metadata_3_analyse[!c(metadata_3_analyse$probe_name %in% c("41IM", "31IM","117IM","45IM")),]
metadata_3_analyse[metadata_3_analyse$probe_name != "",] -> metadata_3_analyse
# metadata_3_analyse<- metadata_3_analyse[!c(metadata_3_analyse$probe_name %in% c("123IM", "96IM",  "95IM", "107IM", "92IM" , "116IM", "117IM", "110IM", "130IM","27IM",
# "34IM","36IM","33IM","39IM","43IM","35IM","64IM","45IM","72IM","47IM","49IM","66IM","61IM","78IM","70IM","87IM","94IM")),]
# samtools sort -n -o ./sorted/41IM.nameSorted.bam ./star/41IMAligned.sortedByCoord.out.bam
# htseq-count --stranded=reverse -f bam -r name ./sorted/117IM.nameSorted.bam ./reference/Homo_sapiens.GRCh38.99.gtf > ./htseq/117IM.txt &

x<-load_DGE(metadata_3_analyse,  "./htseq/") 
output_file<-"./result/NtproBNP_0_pobranie_vs_2_pobranie.csv"
run_3<-run_limma(x, metadata_3_analyse, output_file)
length(which(run_3$adj.P.Val < 0.05))
run_3_deseq<-run_deseq2(x, metadata_3_analyse)
run_3_deseq <- run_3_deseq[!is.na(run_3_deseq$padj),]
length(which(run_3_deseq$padj < 0.05))

write.csv2(run_3_deseq, "deseq_3.csv")
tmp<-run_3_deseq[run_3_deseq$padj < 0.05,]
tmp<-tmp[c(tmp$log2FoldChange < -1 | tmp$log2FoldChange >1),]
nrow(tmp)

# NtproBNP==1 (1 pobranie vs 2 pobranie)
conditions <- list(type=c('Pobranie 1', 'Pobranie 2'), 'wzrost.NtproBNP'=1)
metadata_4_analyse <- cut_metadata(metadata, conditions)
metadata_4_analyse$research<-metadata_4_analyse$type
metadata_4_analyse$research[metadata_4_analyse$research=='Pobranie 1']<-"Control"
metadata_4_analyse$research[metadata_4_analyse$research=='Pobranie 2']<-"Searched"
metadata_4_analyse<- metadata_4_analyse[metadata_4_analyse$probe_name != "",]
metadata_4_analyse<- metadata_4_analyse[!c(metadata_4_analyse$probe_name %in% c("63IM")),]
# "46IM", "54IM", "63IM", "68IM", "84IM", "102IM", "113IM", "111IM", "125IM"



x<-load_DGE(metadata_4_analyse,  "./htseq/") 
output_file<-"./result/NtproBNP_1_pobranie_vs_2_pobranie.csv"
run_4<-run_limma(x, metadata_4_analyse, output_file)
length(which(run_4$adj.P.Val < 0.05))
run_4_deseq<-run_deseq2(x, metadata_4_analyse)
run_4_deseq <- run_4_deseq[!is.na(run_4_deseq$padj),]
length(which(run_4_deseq$padj < 0.05))

write.csv2(run_4_deseq, "deseq_4.csv")
tmp<-run_4_deseq[run_4_deseq$padj < 0.05,]
tmp<-tmp[c(tmp$log2FoldChange < -1 | tmp$log2FoldChange >1),]
nrow(tmp)


result<-add_genes_info(run_4_deseq, ah)
result[result$padj<0.05,] -> result_sign
enrich<-enrichGO(result_sign)
enrich<- as.data.frame(enrich)
write.csv2(enrich, "./result/NtproBNP_1_pobranie_vs_2_pobranie_enrich_padj_0_05.csv")
result_sign[c(result_sign$log2FoldChange<-1 | result_sign$log2FoldChange>1),] -> result_sign
enrich_4<-enrichGO(result_sign)
enrich_4<- as.data.frame(enrich)
write.csv2(enrich, "./result/NtproBNP_1_pobranie_vs_2_pobranie_enrich_padj_0_05_logfc.csv")

# 2 pobranie (NtproBNP==1, 0)
conditions <- list(type='Pobranie 2', 'wzrost.NtproBNP'=c(0,1))
metadata_5_analyse <- cut_metadata(metadata, conditions)
metadata_5_analyse$research<-metadata_5_analyse$'wzrost.NtproBNP'
metadata_5_analyse$research[metadata_5_analyse$research=='0']<-"Control"
metadata_5_analyse$research[metadata_5_analyse$research=='1']<-"Searched"
metadata_5_analyse<- metadata_5_analyse[metadata_5_analyse$probe_name != "",]
metadata_5_analyse<- metadata_5_analyse[!c(metadata_5_analyse$probe_name %in% c("63IM")),]
x<-load_DGE(metadata_5_analyse,  "./htseq2/") 
y<-get_voom(x, metadata_5_analyse)
#res<-test(metadata_5_analyse, y)

x<-load_DGE(metadata_5_analyse,  "./htseq2/") 
# x <- calcNormFactors(x)
run_5<-run_limma(x, metadata_5_analyse, output_file)
length(which(run_5$adj.P.Val < 0.05))
tmp<-run_5[run_5$adj.P.Val < 0.05,]
tmp<-tmp[c(tmp$logFC < -1 | tmp$logFC >1),]
nrow(tmp)

run_5_deseq<-run_deseq2(x, metadata_5_analyse)
run_5_deseq <- run_5_deseq[!is.na(run_5_deseq$padj),]
length(which(run_5_deseq$padj < 0.05))

write.csv2(run_5_deseq, "deseq_5.csv")
tmp<-run_5_deseq[run_5_deseq$padj < 0.05,]
tmp<-tmp[c(tmp$log2FoldChange < -1 | tmp$log2FoldChange >1),]
nrow(tmp)


result<-add_genes_info(run_5_deseq, ah)
result[result$padj<0.05,] -> result_sign
# enrich<-enrichGO(result_sign)
# enrich<- as.data.frame(enrich)
write.csv2(enrich, "./result/pobranie2_NtproBNP_1_vs_0_enrich_padj_0_05.csv")
result_sign[c(result_sign$log2FoldChange<-1 | result_sign$log2FoldChange>1),] -> result_sign
  gene_symbols_5<-rownames(result_sign)
enrich_5<-enrichGO(result_sign)
enrich_5<- as.data.frame(enrich_5)
write.csv2(enrich, "./result/pobranie2_NtproBNP_1_vs_0_enrich_padj_0_05_logfc.csv")

enrich_5<-enrichGO(result_sign)
enrich_5_df<- as.data.frame(enrich_5)
pdf("goEnrich_5.pdf", width = 7, height = 7)
ego <- pairwise_termsim(enrich_5)
emapplot(ego)
dev.off()



result<-add_genes_info(run_4_deseq, ah)
result[result$padj<0.05,] -> result_sign
# enrich<-enrichGO(result_sign)
# enrich<- as.data.frame(enrich)
write.csv2(enrich, "./result/NtproBNP_1_pobranie_vs_2_pobranie_enrich_padj_0_05_logfc.csv")
result_sign[c(result_sign$log2FoldChange<-1 | result_sign$log2FoldChange>1),] -> result_sign
enrich_4<-enrichGO(result_sign)
enrich_4_df<- as.data.frame(enrich_4)
write.csv2(enrich_4_df, "./result/NtproBNP_1_pobranie_vs_2_pobranie_enrich_padj_0_05_logfc.csv")

pdf("goEnrich_4.pdf", width = 7, height = 7)
ego <- pairwise_termsim(enrich_4)
emapplot(ego)
dev.off()


  gene_symbols_4<-rownames(result_sign)


library(clusterProfiler)
library(enrichplot)

geneClusters <- list(
  pobranie2_NtproBNP_1_vs_0 = gene_symbols_5,
  NtproBNP_1_pobranie_vs_2_pobranie = gene_symbols_4
)

cc <- compareCluster(
  geneCluster = geneClusters,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,   # zmień jeśli trzeba
  ont = "BP"
)
cc <- pairwise_termsim(cc)
pdf("goEnrich.pdf", width = 7, height = 7)

emapplot(cc, pie = "count")
dev.off()



enrich_5<-enrichGO(result_sign)
enrich_5<- as.data.frame(enrich_5)
write.csv2(enrich, "./result/pobranie2_NtproBNP_1_vs_0_enrich_padj_0_05_logfc.csv")


library(clusterProfiler)
library(enrichplot)



pdf("sensitivity_analysis_68IM.pdf", width = 7, height = 7)

plot(run_5[common, "logFC"],
     run_5_no_68[common, "logFC"],
     xlab = "logFC full",
     ylab = "logFC no 68IM")

abline(0,1,col="red")
dev.off()

cor(run_5[common, "logFC"],
    run_5_no_68[common, "logFC"])



sig_full <- rownames(run_5)[run_5$adj.P.Val < 0.05]
sig_no68 <- rownames(run_5_no_68)[run_5_no_68$adj.P.Val < 0.05]

length(setdiff(sig_full, sig_no68))
length(setdiff(sig_no68, sig_full))



tmp<-merge(run_5[which(run_5$adj.P.Val < 0.05),], count_5, by=0)
write.csv2(tmp, "run_5.csv")

result<-add_genes_info(run_5)
result_fixed <- result |>
  mutate(across(where(is.list), ~ sapply(., paste, collapse = ";")))
write.csv2(result_fixed, "./result/pobranie2_NtproBNP_1_vs_0_with_genes.csv")


# 1 pobranie vs 2 pobranie

conditions <- list(type=c('Pobranie 1', 'Pobranie 2'),'wzrost.NtproBNP'=0)
metadata_5_2_analyse <- cut_metadata(metadata, conditions)
metadata_5_2_analyse$research<-metadata_5_2_analyse$type
metadata_5_2_analyse$research[metadata_5_2_analyse$research=='Pobranie 1']<-"Control"
metadata_5_2_analyse$research[metadata_5_2_analyse$research=='Pobranie 2']<-"Searched"
metadata_5_2_analyse<- metadata_5_2_analyse[metadata_5_2_analyse$probe_name != "",]
metadata_5_2_analyse<- metadata_5_2_analyse[!c(metadata_5_2_analyse$probe_name %in% c("63IM")),]


x<-load_DGE(metadata_5_2_analyse,  "./htseq/") 
output_file<-"./result/pobranie1_pobranie2_with_genes.csv"
run_5_2<-run_limma(x, metadata_5_2_analyse, output_file)
length(which(run_5_2$adj.P.Val < 0.05))
tmp<-run_5_2[run_5_2$adj.P.Val < 0.05,]
tmp<-tmp[c(tmp$logFC < -1 | tmp$logFC >1),]
nrow(tmp)

run_5_2_deseq<-run_deseq2(x, metadata_5_2_analyse)
run_5_2_deseq <- run_5_2_deseq[!is.na(run_5_2_deseq$padj),]
length(which(run_5_2_deseq$padj < 0.05))

write.csv2(run_5_2_deseq, "deseq_5_2.csv")
tmp<-run_5_2_deseq[run_5_2_deseq$padj < 0.05,]
tmp<-tmp[c(tmp$log2FoldChange < -1 | tmp$log2FoldChange >1),]
nrow(tmp)


result<-add_genes_info(run_5_2_deseq)
result[result$padj<0.05,] -> result_sign
enrich<-enrichGO(result_sign)
enrich<- as.data.frame(enrich)
write.csv2(enrich, "./result/pobranie1_pobranie2_enrich_padj_0_05.csv"
result_sign[c(result_sign$log2FoldChange<-1 | result_sign$log2FoldChange>1),] -> result_sign
enrich<-enrichGO(result_sign)
enrich<- as.data.frame(enrich)
write.csv2(enrich, "./result/pobranie1_pobranie2_enrich_padj_0_05_logfc.csv")
# 1 pobranie (NtproBNP==1, 0)


conditions <- list(type='Pobranie 1', 'wzrost.NtproBNP'=c(0,1))
metadata_6_analyse <- cut_metadata(metadata, conditions)
metadata_6_analyse$research<-metadata_6_analyse$'wzrost.NtproBNP'
metadata_6_analyse$research[metadata_6_analyse$research=='0']<-"Control"
metadata_6_analyse$research[metadata_6_analyse$research=='1']<-"Searched"
metadata_6_analyse<- metadata_6_analyse[metadata_6_analyse$probe_name != "",]
metadata_6_analyse<- metadata_6_analyse[!c(metadata_6_analyse$probe_name %in% c("63IM")),]


x<-load_DGE(metadata_6_analyse,  "./htseq/") 
output_file<-"./result/pobranie1_NtproBNP_1_vs_0_with_genes.csv"
run_6<-run_limma(x, metadata_6_analyse, output_file)
length(which(run_6$adj.P.Val < 0.05))
tmp<-run_6[run_6$adj.P.Val < 0.05,]
tmp<-tmp[c(tmp$logFC < -1 | tmp$logFC >1),]
nrow(tmp)

run_6_deseq<-run_deseq2(x, metadata_6_analyse)
run_6_deseq <- run_6_deseq[!is.na(run_6_deseq$padj),]
length(which(run_6_deseq$padj < 0.05))

write.csv2(run_6_deseq, "deseq_6.csv")
tmp<-run_6_deseq[run_6_deseq$padj < 0.05,]
tmp<-tmp[c(tmp$log2FoldChange < -1 | tmp$log2FoldChange >1),]
nrow(tmp)

           
result<-add_genes_info(run_6_deseq,ah)
result[result$padj<0.05,] -> result_sign
enrich<-enrichGO(result_sign)
enrich<- as.data.frame(enrich)
write.csv2(enrich, "./result/pobranie1_NtproBNP_1_vs_0_enrich_padj_0_05.csv")
result_sign[c(result_sign$log2FoldChange<-1 | result_sign$log2FoldChange>1),] -> result_sign
enrich_6<-enrichGO(result_sign)
enrich_6<- as.data.frame(enrich_6)
write.csv2(enrich, "./result/pobranie1_NtproBNP_1_vs_0_enrich_padj_0_05_logfc.csv")
           
# 3,4 pobranie (NtproBNP==1, 0)

conditions <- list(type=c('Pobranie 3', 'Pobranie 4'), 'wzrost.NtproBNP'=c(0,1))
metadata_7_analyse <- cut_metadata(metadata, conditions)
metadata_7_analyse$research<-metadata_7_analyse$'wzrost.NtproBNP'
metadata_7_analyse$research[metadata_7_analyse$research=='0']<-"Control"
metadata_7_analyse$research[metadata_7_analyse$research=='1']<-"Searched"
metadata_7_analyse<- metadata_7_analyse[metadata_7_analyse$probe_name != "",]
metadata_7_analyse<- metadata_7_analyse[!c(metadata_7_analyse$probe_name %in% c("63IM")),]

x<-load_DGE(metadata_7_analyse,  "./htseq2/") 
output_file<-"./result/pobranie1_NtproBNP_1_vs_0_with_genes.csv"
run_7<-run_limma(x, metadata_7_analyse, output_file)
length(which(run_7$adj.P.Val < 0.05))
run_7_deseq<-run_deseq2(x, metadata_7_analyse)
run_7_deseq <- run_7_deseq[!is.na(run_7_deseq$padj),]
length(which(run_7_deseq$padj < 0.05))

write.csv2(run_7_deseq, "deseq_7.csv")
tmp<-run_7_deseq[run_7_deseq$padj < 0.05,]
tmp<-tmp[c(tmp$log2FoldChange < -1 | tmp$log2FoldChange >1),]
nrow(tmp)


conditions <- list(type=c('Pobranie 3', 'Pobranie 4'), 'wzrost.NtproBNP'=c(0,1))
metadata_8_analyse <- cut_metadata(metadata, conditions)
metadata_8_analyse$research<-metadata_8_analyse$'wzrost.NtproBNP'
metadata_8_analyse$research[metadata_8_analyse$research=='0']<-"Control"
metadata_8_analyse$research[metadata_8_analyse$research=='1']<-"Searched"
metadata_8_analyse<- metadata_8_analyse[metadata_8_analyse$probe_name != "",]
metadata_8_analyse<- metadata_8_analyse[!c(metadata_8_analyse$probe_name %in% c("63IM")),]

x<-load_DGE(metadata_8_analyse,  "./htseq2/") 

  files<-sapply(metadata_8_analyse$probe_name,function(x) paste0(folder,x,".txt"))

                
output_file<-"./result/pobranie_3_4_NtproBNP_1_vs_0_with_genes.csv"
run_8<-run_limma(x, metadata_8_analyse, output_file)
length(which(run_8$adj.P.Val < 0.05))
run_8_deseq<-run_deseq2(x, metadata_8_analyse)
run_8_deseq <- run_8_deseq[!is.na(run_8_deseq$padj),]
length(which(run_8_deseq$padj < 0.05))

write.csv2(run_8_deseq, "deseq_8.csv")



x<-load_DGE(metadata_5_analyse,  "./htseq/") 
output_file<-"./result/pobranie2_NtproBNP_1_vs_0.csv"
run_5<-run_limma(x, metadata_5_analyse, output_file)
length(which(run_5$adj.P.Val < 0.05))
result<-add_genes_info(run_5)
result_fixed <- result |>
  mutate(across(where(is.list), ~ sapply(., paste, collapse = ";")))
write.csv2(result_fixed, "./result/pobranie2_NtproBNP_1_vs_0_with_genes.csv")

enrich<-enrichGO(result$Row.names)
enrich<- as.data.frame(enrich)
write.csv2(enrich, "./result/pobranie2_NtproBNP_1_vs_0_enrich.csv")

result_significant<-result[result$adj.P.Val<0.05,]
enrich<-enrichGO(result_significant)
enrich<- as.data.frame(enrich)
write.csv2(enrich, "./result/pobranie2_NtproBNP_1_vs_0_enrich_result_significant.csv")




result_significant_up<-result_significant[result_significant$logFC>0,]
enrich<-enrichGO(result_significant_up)
enrich<- as.data.frame(enrich)
write.csv2(enrich, "./result/pobranie2_NtproBNP_1_vs_0_enrich_result_result_significant_up.csv")



result_significant_down<-result_significant[result_significant$logFC < 0,]
enrich<-enrichGO(result_significant_down, "BP")
enrich<- as.data.frame(enrich)
write.csv2(enrich, "./result/pobranie2_NtproBNP_1_vs_0_enrich_result_result_significant_down_BP.csv")

enrich<-enrichGO(result_significant_down, "MF")
enrich<- as.data.frame(enrich)
write.csv2(enrich, "./result/pobranie2_NtproBNP_1_vs_0_enrich_result_result_significant_down_MF.csv")

tmp<-merge(run_5[which(run_5$adj.P.Val < 0.05),],run_6, by=0)
colnames(tmp)<-c("Row.names_0",paste0(colnames(run_5), "_pobr2"),paste0(colnames(run_6), "_pobr1"))
rownames(tmp)<-tmp$Row.names_0
tmp<-merge(tmp,run_8, by=0)
rownames(tmp)<-tmp$Row.names
colnames(tmp)<-c("Row.names_0", "Row.names_1",paste0(colnames(run_5), "_pobr2"),paste0(colnames(run_6), "_pobr1"),paste0(colnames(run_6), "_pobr3_4"))
tmp<-merge(tmp,run_7, by=0)
rownames(tmp)<-tmp$Row.names
colnames(tmp)<-c("Row.names_0", "Row.names_1", "Row.names_2",paste0(colnames(run_5), "_pobr2"),paste0(colnames(run_6), "_pobr1"),paste0(colnames(run_6), "_pobr3_4"),paste0(colnames(run_6), "_pobr2_3_4"))


samtools sort -n -o ./sorted/46IM.nameSorted.bam ./star/46IMAligned.sortedByCoord.out.bam &
samtools sort -n -o ./sorted/54IM.nameSorted.bam ./star/54IMAligned.sortedByCoord.out.bam &
samtools sort -n -o ./sorted/68IM.nameSorted.bam ./star/68IMAligned.sortedByCoord.out.bam &
samtools sort -n -o ./sorted/84IM.nameSorted.bam ./star/84IMAligned.sortedByCoord.out.bam &
samtools sort -n -o ./sorted/102IM.nameSorted.bam ./star/102IMAligned.sortedByCoord.out.bam &
samtools sort -n -o ./sorted/113IM.nameSorted.bam ./star/113IMAligned.sortedByCoord.out.bam &
samtools sort -n -o ./sorted/111IM.nameSorted.bam ./star/111IMAligned.sortedByCoord.out.bam &
samtools sort -n -o ./sorted/125IM.nameSorted.bam ./star/125IMAligned.sortedByCoord.out.bam &
111 <- error

htseq-count --stranded=reverse -f bam -r name ./sorted/46IM.nameSorted.bam ./reference/Homo_sapiens.GRCh38.99.gtf > ./htseq/46IM.txt &
htseq-count --stranded=reverse -f bam -r name ./sorted/54IM.nameSorted.bam ./reference/Homo_sapiens.GRCh38.99.gtf > ./htseq/54IM.txt &
htseq-count --stranded=reverse -f bam -r name ./sorted/84IM.nameSorted.bam ./reference/Homo_sapiens.GRCh38.99.gtf > ./htseq/84IM.txt &
htseq-count --stranded=reverse -f bam -r name ./sorted/102IM.nameSorted.bam ./reference/Homo_sapiens.GRCh38.99.gtf > ./htseq/102IM.txt &
htseq-count --stranded=reverse -f bam -r name ./sorted/113IM.nameSorted.bam ./reference/Homo_sapiens.GRCh38.99.gtf > ./htseq/113IM.txt &

htseq-count --stranded=reverse -f bam -r name ./sorted/68IM.nameSorted.bam ./reference/Homo_sapiens.GRCh38.99.gtf > ./htseq/68IM.txt &
htseq-count --stranded=reverse -f bam -r name ./sorted/111IM.nameSorted.bam ./reference/Homo_sapiens.GRCh38.99.gtf > ./htseq/111IM.txt &
htseq-count --stranded=reverse -f bam -r name ./sorted/125IM.nameSorted.bam ./reference/Homo_sapiens.GRCh38.99.gtf > ./htseq/125IM.txt &



samtools sort -n -o ./sorted/91IM.nameSorted.bam ./star/91IMAligned.sortedByCoord.out.bam &
htseq-count --stranded=reverse -f bam -r name ./sorted/91IM.nameSorted.bam ./reference/Homo_sapiens.GRCh38.99.gtf > ./htseq/91IM.txt &




                
output_file<-"./result/pobranie2_NtproBNP_1_vs_0.csv"
y<-get_voom(x, metadata_5_analyse)
plotPCA(y,metadata_5_analyse, "tmp.svg")
Metadata_5_analyse[metadata_5_analyse$probe_name %in% c("130IM", "78IM", "123IM", "96IM", "111IM", "95IM", "87IM", "95IM", "33IM", "110IM", "102IM", "46IM"),]
control<-metadata_5_analyse$probe_name[!c(metadata_5_analyse$probe_name %in% c("130IM", "78IM", "123IM", "96IM", "111IM", "95IM", "87IM", "95IM", "33IM", "110IM", "102IM", "46IM"))] 
searched<-c("130IM", "78IM", "123IM", "96IM", "111IM", "95IM", "87IM", "95IM", "33IM", "110IM", "102IM", "46IM")

metadata_5_analyse$research[metadata_5_analyse$probe_name %in% control]<-"Control"
metadata_5_analyse$research[metadata_5_analyse$probe_name %in% searched]<-"Searched"
pdf("MDS_sex.pdf", width = 7, height = 6)
plotMDS(y, col = as.numeric(as.factor(metadata_5_analyse$s)))
dev.off()

pdf("MDS_research.png", width = 7, height = 6)
plotMDS(y, col = as.numeric(as.factor(metadata_5_analyse$research)))
dev.off()


png("MDS_research.png", width = 800, height = 600)
plotMDS(y, col = as.factor(metadata_5_analyse$research))
dev.off()

mds <- plotMDS(y, plot = FALSE)
mds$var.explained
df <- data.frame(
  Dim1 = mds$x,
  Dim2 = mds$y,
  sex = metadata_5_analyse$s,
  research = metadata_5_analyse$research
)
pca <- prcomp(t(y$E))

p1<-plot(pca$x[,1], pca$x[,2], col = as.factor(metadata_5_analyse$s))
p2<-plot(pca$x[,1], pca$x[,2], col = as.factor(metadata_5_analyse$research))

ggsave("p1.png",
       plot = p1,
       width = 8,
       height = 6,
       dpi = 100)

ggsave("p2.png",
       plot = p2,
       width = 8,
       height = 6,
       dpi = 100)


v_adj <- removeBatchEffect(y, covariates = model.matrix(~ s, data = metadata_5_analyse))
pdf("MDS_research.png", width = 7, height = 6)

plotMDS(v_adj, col = as.numeric(as.factor(metadata_5_analyse$research)))
dev.off()





                metadata_5_analyse<- metadata_5_analyse[!c(metadata_5_analyse$probe_name %in% c("63IM")),]

x<-load_DGE(metadata_5_analyse,  "./htseq2/") 
x$samples$lib.size
y<-get_voom(x, metadata_5_analyse)
count_5<-as.data.frame(x$counts)

                

                
