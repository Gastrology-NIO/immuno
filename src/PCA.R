
library(vegan)

plotPCA<-function(dds,metadata, save_path, model){
  vsd <- vst(dds, blind = FALSE)
  mat <- counts(dds, normalized = TRUE)
  
  mat <- mat[rowMeans(mat) > 1, ]
  
  mat_t <- t(mat)
  dist_bc <- vegdist(mat_t, method = "euclidean")
  
  meta <- as.data.frame(colData(dds))
  
  # 6. PERMANOVA
  adonis_res <- adonis2(
    dist_bc ~ research,
    data = meta,
    permutations = 999
  )
  pval <- adonis_res$`Pr(>F)`[1]
  
  r2 <- adonis_res$R2[1]
  label_text <- paste0(
    "PERMANOVA\nR² = ",
    round(r2, 3),
    "\np = ",
    signif(pval, 3)
)
  

  
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
p <- ggplot(pcaData,
       aes(PC1, PC2,
           color = research)) +

    geom_point(size = 4) +

    stat_ellipse(
      aes(fill = research),
      geom = "polygon",
      alpha = 0.2,
      color = NA
    ) +



    annotate(
        "text",
        x = Inf,
        y = Inf,
      
      label = label_text,
        hjust = 1.1,
        vjust = 1.5,
        size = 5
    ) +

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


}
