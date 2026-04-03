library(edgeR)
load_data_meta<- function(path_metadata){
  metadata<-read.csv(path_metadata, sep=";")
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
  
load_DGE <- function(metadata, folder){
  files<-sapply(metadata$probe_name,function(x) paste0(folder,x,".txt"))
  x <- readDGE(files, columns=c(1,2))
  samplenames<-sapply(colnames(x), function(x) strsplit(x, "/")[[1]][3])
  colnames(x) <- samplenames
  x$samples$research <- metadata$research
  x$samples$s <- metadata$s
  return(x)
}

                      
cut_metadata<-function(metadata, conditions){
  for (i in 1:length(conditions)){
    name=names(conditions)[i]
    metadata<-metadata[metadata[,name] %in% conditions[name][[1]],]
  }
  return(metadata)
}
library(AnnotationHub)
library(ensembldb)

add_genes_info<-function(data){
  gene_symbols<-rownames(data)
  ah <- AnnotationHub()
  query(ah, c("Homo sapiens", "EnsDb"))
  edb <- ah[["AH119325"]]
  genes <- genes(edb, filter = GeneIdFilter(gene_symbols))
  genes_info<-as.data.frame(genes)
  data<-merge(genes_info, data, by=0)
  return(data)
}
