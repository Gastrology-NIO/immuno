
metadata<-load_data_meta("./data/metadata_jag.csv")

# checkpoint1 VS control
name<-"checkpoint1_vs_control"
folder<-paste0("./result/",name,"/")

conditions <- list(type=c('checkpoint 1', 'control'))
metadata_0_analyse <- cut_metadata(metadata, conditions)
metadata_0_analyse$research<-metadata_0_analyse$type
metadata_0_analyse$research[metadata_0_analyse$type=='control']<-"Control"
metadata_0_analyse$research[metadata_0_analyse$type=='checkpoint 1']<-"Searched"
metadata_0_analyse<- metadata_0_analyse[metadata_0_analyse$probe_name != "",]
run_analyse(folder, metadata_0_analyse, name, analyse_pairs=F, analyse_sex=T, max_genes=20, max_categories=5)
run_analyse(folder, metadata_0_analyse, paste0(name, "_2"), analyse_pairs=F, analyse_sex=T)
run_analyse(folder, metadata_0_analyse, paste0(name, "_3"), analyse_pairs=F, analyse_sex=T, max_categories=5)
deconvolution_difference(name, metadata_0_analyse)

# additional info: difference.days
# Checkpoint 1 (under 3 months vs over 2 years)
name<-"checkpoint1_3months_vs_Over2Years"
folder<-paste0("./result/",name,"/")
conditions <- list(type='checkpoint 1', 'time.of.OS'=c('under 3 months', 'over 2 years'))
metadata_1_analyse <- cut_metadata(metadata, conditions)
metadata_1_analyse$research<-metadata_1_analyse$'time.of.OS'
metadata_1_analyse$research[metadata_1_analyse$research=='under 3 months']<-"Control"
metadata_1_analyse$research[metadata_1_analyse$research=='over 2 years']<-"Searched"

run_analyse(folder, metadata_1_analyse, name, analyse_pairs=F, analyse_sex=T)
deconvolution_difference(name, metadata_1_analyse)

name<-"checkpoint1_3months_vs_Over2Years"
folder<-paste0("./result/",name,"/")
conditions <- list(type='checkpoint 1')
metadata_checkpoint1 <- cut_metadata(metadata, conditions)
x<-load_DGE(metadata_checkpoint1,  "./htseq/") 
dge <- calcNormFactors(x)

# log-CPM
expr <- cpm(
    dge,
    log = TRUE,
    prior.count = 1
)
rownames(expr) <- sub("\\..*$", "", rownames(expr))
colnames(expr) <- basename(colnames(expr))

# NtproBNP==1 (checkpoint 1 vs 2 checkpoint 2)
name<-"NtproBNP_1_checkpoint1_vs_checkpoint2"
folder<-paste0("./result/",name,"/")
conditions <- list(type=c('checkpoint 1', 'checkpoint 2'), 'wzrost.NtproBNP'=1)
metadata_4_analyse <- cut_metadata(metadata, conditions)
metadata_4_analyse$research<-metadata_4_analyse$type
metadata_4_analyse$research[metadata_4_analyse$research=='checkpoint1']<-"Control"
metadata_4_analyse$research[metadata_4_analyse$research=='checkpoint2']<-"Searched"
metadata_4_analyse<- metadata_4_analyse[metadata_4_analyse$probe_name != "",]
run_analyse(folder, metadata_4_analyse, name, analyse_pairs=T, analyse_sex=F)
deconvolution_difference(name, metadata_4_analyse)



# checkpoint 2 (NtproBNP==1, 0)
name<-"checkpoint2_NtproBNP_1_vs_0"
folder<-paste0("./result/",name,"/")
conditions <- list(type='checkpoint 2', 'wzrost.NtproBNP'=c(0,1))
metadata_5_analyse <- cut_metadata(metadata, conditions)
metadata_5_analyse$research<-metadata_5_analyse$'wzrost.NtproBNP'
metadata_5_analyse$research[metadata_5_analyse$research=='0']<-"Control"
metadata_5_analyse$research[metadata_5_analyse$research=='1']<-"Searched"
metadata_5_analyse<- metadata_5_analyse[metadata_5_analyse$probe_name != "",]
run_analyse(folder, metadata_5_analyse, name, analyse_pairs=F, analyse_sex=T)
deconvolution_difference(name, metadata_5_analyse)


# checkpoint1 (NtproBNP==1, 0)
name<-"checkpoint1_NtproBNP_1_vs_0"
folder<-paste0("./result/",name,"/")
conditions <- list(type='checkpoint 1', 'wzrost.NtproBNP'=c(0,1))
metadata_2_analyse <- cut_metadata(metadata, conditions)
metadata_2_analyse$research<-metadata_2_analyse$'wzrost.NtproBNP'
metadata_2_analyse$research[metadata_2_analyse$research=='0']<-"Control"
metadata_2_analyse$research[metadata_2_analyse$research=='1']<-"Searched"
metadata_2_analyse<- metadata_2_analyse[metadata_2_analyse$probe_name != "",]
run_analyse(folder, metadata_2_analyse, name, analyse_pairs=F, analyse_sex=T)
deconvolution_difference(name, metadata_2_analyse)
