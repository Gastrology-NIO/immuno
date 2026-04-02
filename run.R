
source("./src/load_data.R")


metadata<-load_data_meta("./data/metadata_jag.csv")
gene_counts<-load_data_counts("./htseq/", "./data/gene_counts.csv")


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
metadata_3_analyse<- metadata_3_analyse[!c(metadata_3_analyse$probe_name %in% c("123IM", "96IM",  "95IM", "107IM", "92IM" , "116IM", "117IM", "110IM", "130IM","27IM",
"34IM","36IM","33IM","39IM","43IM","35IM","64IM","45IM","72IM","47IM","49IM","66IM","61IM","78IM","70IM","87IM","94IM")),]
# samtools sort -n -o ./sorted/41IM.nameSorted.bam ./star/41IMAligned.sortedByCoord.out.bam
# htseq-count --stranded=reverse -f bam -r name ./sorted/RNA_41IM.nameSorted.bam ./reference/omo_sapiens.GRCh38.99.gtf > ./htseq/41IM.txt

htseq-count --stranded=reverse -f bam -r name ./sorted/117IM.nameSorted.bam ./reference/Homo_sapiens.GRCh38.99.gtf > ./htseq/117IM.txt &
htseq-count --stranded=reverse -f bam -r name ./sorted/45IM.nameSorted.bam ./reference/Homo_sapiens.GRCh38.99.gtf > ./htseq/45IM.txt &

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
x<-load_DGE(metadata_4_analyse,  "./htseq/") 
output_file<-"./result/NtproBNP_1_pobranie_vs_2_pobranie.csv"
run_4<-run_limma(x, metadata_4_analyse, output_file)
length(which(run_4$adj.P.Val < 0.05))


# 2 pobranie (NtproBNP==1, 0)
conditions <- list(type='Pobranie 2', 'wzrost.NtproBNP'=c(0,1))
metadata_5_analyse <- cut_metadata(metadata, conditions)
metadata_5_analyse$research<-metadata_5_analyse$'wzrost.NtproBNP'
metadata_5_analyse$research[metadata_5_analyse$research=='0']<-"Control"
metadata_5_analyse$research[metadata_5_analyse$research=='1']<-"Searched"
metadata_5_analyse<- metadata_5_analyse[metadata_5_analyse$probe_name != "",]
x<-load_DGE(metadata_5_analyse,  "./htseq/") 
output_file<-"./result/pobranie2_NtproBNP_1_vs_0.csv"
run_4<-run_limma(x, metadata_5_analyse, output_file)
length(which(run_4$adj.P.Val < 0.05))

