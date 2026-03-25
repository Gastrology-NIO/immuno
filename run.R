
source("./src/load_data.R")


metadata<-load_data_meta("./data/meta_jag.csv")
gene_counts<-load_data_counts("./htseq/")


# 1 pobranie (3 miesiące vs powyżej 2 lat)
conditions <- list(type='Pobranie 1', 'time of OS'=c('poniżej 3 miesięcy', '> 2 lata'))
metadata_1_analyse <- cut_metadata(metadata, conditions)
metadata_1_analyse$research<-metadata_1_analyse$'time of OS'

# 1 pobranie (I linia vs kolejna linia)
conditions <- list(type='Pobranie 1', 'linia'=c('I linia', ''))

# NtproBNP==0 (1 pobranie vs 2 pobranie)
conditions <- list(type=c('Pobranie 1', 'Pobranie 2'), 'wzrost NtproBNP'=0)

# NtproBNP==1 (1 pobranie vs 2 pobranie)
conditions <- list(type=c('Pobranie 1', 'Pobranie 2'), 'wzrost NtproBNP'=1)

# 2 pobranie (NtproBNP==0 vs NtproBNP==1)
conditions <- list(type=c('Pobranie 2'), 'wzrost NtproBNP'=c(0, 1))
