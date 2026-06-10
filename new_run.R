
# 1 pobranie VS kontrola
name<-"pobranie1_vs_kontrola"
folder<-paste0("./result/",name,"/")

conditions <- list(type=c('Pobranie 1', 'kontrola'))
metadata_0_analyse <- cut_metadata(metadata, conditions)
metadata_0_analyse$research<-metadata_0_analyse$type
metadata_0_analyse$research[metadata_0_analyse$type=='kontrola']<-"Control"
metadata_0_analyse$research[metadata_0_analyse$type=='Pobranie 1']<-"Searched"
metadata_0_analyse<- metadata_0_analyse[metadata_0_analyse$probe_name != "",]
run_analyse(folder, metadata_0_analyse, name, analyse_pairs=F, analyse_sex=T, max_genes=20)


folder<-"./result/pobranie1_vs_kontrola/"

# 1 pobranie (3 miesiące vs powyżej 2 lat)

name<-"1pobranie_3miesiace_vs_powyżej_2_lata"
folder<-paste0("./result/",name,"/")

conditions <- list(type='Pobranie 1', 'time.of.OS'=c('poni\xbfej 3 miesi\xeacy', '> 2 lata'))
metadata_1_analyse <- cut_metadata(metadata, conditions)
metadata_1_analyse$research<-metadata_1_analyse$'time.of.OS'
metadata_1_analyse$research[metadata_1_analyse$research=='poni\xbfej 3 miesi\xeacy']<-"Control"
metadata_1_analyse$research[metadata_1_analyse$research=='> 2 lata']<-"Searched"

run_analyse(folder, metadata_1_analyse, name, analyse_pairs=F, analyse_sex=T)


# NtproBNP==1 (1 pobranie vs 2 pobranie)
name<-"NtproBNP_1_1_pobranie_vs_2_pobranie"
folder<-paste0("./result/",name,"/")

conditions <- list(type=c('Pobranie 1', 'Pobranie 2'), 'wzrost.NtproBNP'=1)
metadata_4_analyse <- cut_metadata(metadata, conditions)
metadata_4_analyse$research<-metadata_4_analyse$type
metadata_4_analyse$research[metadata_4_analyse$research=='Pobranie 1']<-"Control"
metadata_4_analyse$research[metadata_4_analyse$research=='Pobranie 2']<-"Searched"
metadata_4_analyse<- metadata_4_analyse[metadata_4_analyse$probe_name != "",]
metadata_4_analyse<- metadata_4_analyse[!c(metadata_4_analyse$probe_name %in% c("63IM")),]

run_analyse(folder, metadata_4_analyse, name, analyse_pairs=T, analyse_sex=F)



# 2 pobranie (NtproBNP==1, 0)
name<-"2_pobranie_NtproBNP_1_vs_0"
folder<-paste0("./result/",name,"/")

conditions <- list(type='Pobranie 2', 'wzrost.NtproBNP'=c(0,1))
metadata_5_analyse <- cut_metadata(metadata, conditions)
metadata_5_analyse$research<-metadata_5_analyse$'wzrost.NtproBNP'
metadata_5_analyse$research[metadata_5_analyse$research=='0']<-"Control"
metadata_5_analyse$research[metadata_5_analyse$research=='1']<-"Searched"
metadata_5_analyse<- metadata_5_analyse[metadata_5_analyse$probe_name != "",]
metadata_5_analyse<- metadata_5_analyse[!c(metadata_5_analyse$probe_name %in% c("63IM")),]

run_analyse(folder, metadata_5_analyse, name, analyse_pairs=F, analyse_sex=T)


