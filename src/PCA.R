


plotPCA<-function(dds,metadata, save_path){

  vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
pca <- prcomp(t(mat), scale. = FALSE)
pcaData <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    research = colData(vsd)$research,
    patient_id = colData(vsd)$patient_id
)
percentVar <- round(
    100 * (pca$sdev^2 / sum(pca$sdev^2)),
    1
)
p<-ggplot(pcaData,
       aes(PC1, PC2,
           color = research)) +
    geom_point(size = 4) +
    theme_minimal() +
    xlab(paste0("PC1: ", percentVar[1], "%")) +
    ylab(paste0("PC2: ", percentVar[2], "%"))
    ggsave(
    save_path,
    plot = p,
  )
  
#   mat_corr <- limma::removeBatchEffect(
#     assay(vsd),
#     batch = vsd$patient_id
# )
  # pca <- prcomp(t(mat_corr))
#   mds <- plotMDS(x, plot = FALSE)

#   df <- data.frame(
#     Dim1 = mds$x,
#     Dim2 = mds$y,
#     group = metadata$research, 
#     label = metadata$probe_name,
#         research = metadata$research

#   )
  
#   library(ggplot2)
  
# p <- ggplot(df, aes(Dim1, Dim2, color = research)) +
#   geom_point(size = 4) +
#   geom_text(aes(label = label), vjust = -0.5) +
#   theme_minimal()

#     ggsave(
#     save_path,
#     plot = p,
#   )
}
