
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

# 2 pobranie (NtproBNP==1, 0)
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

# NtproBNP==0 (1 pobranie vs 2 pobranie)
conditions <- list(type=c('Pobranie 1', 'Pobranie 2'), 'wzrost.NtproBNP'=0)
metadata_3_analyse <- cut_metadata(metadata, conditions)
metadata_3_analyse$research<-metadata_3_analyse$type
metadata_3_analyse$research[metadata_3_analyse$research=='Pobranie 1']<-"Control"
metadata_3_analyse$research[metadata_3_analyse$research=='Pobranie 2']<-"Searched"
metadata_3_analyse<- metadata_3_analyse[!c(metadata_3_analyse$probe_name %in% c("117IM","45IM")),]
metadata_3_analyse[metadata_3_analyse$probe_name != "",] -> metadata_3_analyse
# metadata_3_analyse<- metadata_3_analyse[!c(metadata_3_analyse$probe_name %in% c("123IM", "96IM",  "95IM", "107IM", "92IM" , "116IM", "117IM", "110IM", "130IM","27IM",
# "34IM","36IM","33IM","39IM","43IM","35IM","64IM","45IM","72IM","47IM","49IM","66IM","61IM","78IM","70IM","87IM","94IM")),]
# samtools sort -n -o ./sorted/41IM.nameSorted.bam ./star/41IMAligned.sortedByCoord.out.bam
# htseq-count --stranded=reverse -f bam -r name ./sorted/117IM.nameSorted.bam ./reference/Homo_sapiens.GRCh38.99.gtf > ./htseq/117IM.txt &

x<-load_DGE(metadata_3_analyse,  "./htseq/") 
output_file<-"./result/NtproBNP_0_pobranie_vs_2_pobranie.csv"
run_3<-run_limma(x, metadata_3_analyse, output_file)
length(which(run_3$adj.P.Val < 0.05))



# NtproBNP==1 (1 pobranie vs 2 pobranie)
conditions <- list(type=c('Pobranie 1', 'Pobranie 2'), 'wzrost.NtproBNP'=1)
metadata_4_analyse <- cut_metadata(metadata, conditions)
metadata_4_analyse$research<-metadata_4_analyse$type
metadata_4_analyse$research[metadata_4_analyse$research=='Pobranie 1']<-"Control"
metadata_4_analyse$research[metadata_4_analyse$research=='Pobranie 2']<-"Searched"
metadata_4_analyse<- metadata_4_analyse[metadata_4_analyse$probe_name != "",]
metadata_4_analyse<- metadata_4_analyse[!c(metadata_4_analyse$probe_name %in% c("117IM","68IM", "63IM")),]
# "46IM", "54IM", "63IM", "68IM", "84IM", "102IM", "113IM", "111IM", "125IM"



x<-load_DGE(metadata_4_analyse,  "./htseq/") 
output_file<-"./result/NtproBNP_1_pobranie_vs_2_pobranie.csv"
run_4<-run_limma(x, metadata_4_analyse, output_file)
length(which(run_4$adj.P.Val < 0.05))


GeneRatio	BgRatio	RichFactor	FoldEnrichment	zScore


# 2 pobranie (NtproBNP==1, 0)
conditions <- list(type='Pobranie 2', 'wzrost.NtproBNP'=c(0,1))
metadata_5_analyse <- cut_metadata(metadata, conditions)
metadata_5_analyse$research<-metadata_5_analyse$'wzrost.NtproBNP'
metadata_5_analyse$research[metadata_5_analyse$research=='0']<-"Control"
metadata_5_analyse$research[metadata_5_analyse$research=='1']<-"Searched"
metadata_5_analyse<- metadata_5_analyse[metadata_5_analyse$probe_name != "",]
metadata_5_analyse<- metadata_5_analyse[!c(metadata_5_analyse$probe_name %in% c("117IM","68IM", "63IM")),]


x<-load_DGE(metadata_5_analyse,  "./htseq/") 
output_file<-"./result/pobranie2_NtproBNP_1_vs_0.csv"
run_5<-run_limma(x, metadata_5_analyse, output_file)
length(which(run_5$adj.P.Val < 0.05))
result<-add_genes_info(run_5)
result_fixed <- result |>
  mutate(across(where(is.list), ~ sapply(., paste, collapse = ";")))
write.csv2(result_fixed, "./result/pobranie2_NtproBNP_1_vs_0_with_genes.csv")

# 1 pobranie (NtproBNP==1, 0)

conditions <- list(type='Pobranie 1', 'wzrost.NtproBNP'=c(0,1))
metadata_6_analyse <- cut_metadata(metadata, conditions)
metadata_6_analyse$research<-metadata_6_analyse$'wzrost.NtproBNP'
metadata_6_analyse$research[metadata_6_analyse$research=='0']<-"Control"
metadata_6_analyse$research[metadata_6_analyse$research=='1']<-"Searched"
metadata_6_analyse<- metadata_6_analyse[metadata_6_analyse$probe_name != "",]
metadata_6_analyse<- metadata_6_analyse[!c(metadata_6_analyse$probe_name %in% c("117IM","68IM", "63IM")),]


x<-load_DGE(metadata_6_analyse,  "./htseq/") 
output_file<-"./result/pobranie1_NtproBNP_1_vs_0_with_genes.csv"
run_6<-run_limma(x, metadata_6_analyse, output_file)
length(which(run_6$adj.P.Val < 0.05))

# 3,4 pobranie (NtproBNP==1, 0)

conditions <- list(type=c('Pobranie 3', 'Pobranie 4'), 'wzrost.NtproBNP'=c(0,1))
metadata_7_analyse <- cut_metadata(metadata, conditions)
metadata_7_analyse$research<-metadata_7_analyse$'wzrost.NtproBNP'
metadata_7_analyse$research[metadata_7_analyse$research=='0']<-"Control"
metadata_7_analyse$research[metadata_7_analyse$research=='1']<-"Searched"
metadata_7_analyse<- metadata_7_analyse[metadata_7_analyse$probe_name != "",]



x<-load_DGE(metadata_7_analyse,  "./htseq/") 
output_file<-"./result/pobranie1_NtproBNP_1_vs_0_with_genes.csv"
run_7<-run_limma(x, metadata_7_analyse, output_file)
length(which(run_7$adj.P.Val < 0.05))

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
