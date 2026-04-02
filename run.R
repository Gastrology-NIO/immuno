
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

# 1 pobranie (I linia vs kolejna linia)
conditions <- list(type='Pobranie 1', 'linia'=c('I linia', ''))

# NtproBNP==0 (1 pobranie vs 2 pobranie)
conditions <- list(type=c('Pobranie 1', 'Pobranie 2'), 'wzrost NtproBNP'=0)

# NtproBNP==1 (1 pobranie vs 2 pobranie)
conditions <- list(type=c('Pobranie 1', 'Pobranie 2'), 'wzrost NtproBNP'=1)

# 2 pobranie (NtproBNP==0 vs NtproBNP==1)
conditions <- list(type=c('Pobranie 2'), 'wzrost NtproBNP'=c(0, 1))
