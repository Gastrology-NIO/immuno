
load_data_meta<- function(path_metadata){
  metadata<-read.csv(path_metadata, sep=";")
  #metadata_tmp<-metadata[c(metadata$type %in% c("Pobranie1", "kontrola")),]
  metadata$probe_name<-gsub("S15156", "RNA_S15156", metadata$probe_name)
  return(metadata)
}

load_data_counts<- function(path_htseq, gene_no_count_save_file){
  experiments_all<-data.frame()
  experiments_all<-read.csv(paste(path_htseq,metadata_tmp$probe_name[1], ".txt", sep=""),header=FALSE, sep="\t")
  names(experiments_all)<-c('gene', metadata_tmp$probe_name[1])
  for (i in metadata_tmp$probe_name[2:length(metadata_tmp$probe_name)]){
    exp_id<-i
    experiments<-read.csv(paste(path_htseq,i, ".txt", sep=""),header=FALSE, sep="\t")
    names(experiments)<-c('gene', exp_id)
    experiments_all<-merge(experiments_all, experiments)
  }
  write.csv2(experiments_all, gene_no_count_save_file)
  return(experiments_all)
}
  
cut_metadata<-function(metadata, conditions){
  for (i in 1:length(conditions)){
    print(i)
    df[df[,c("sample_type")]  %in% c("skora"),]
    name=names(conditions)[i]
    metadata<-metadata[metadata[,name] %in% conditions[name]]
  }
  return(metadata)
}
