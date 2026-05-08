


plotPCA<-function(x,metadata, save_path){
  mds <- plotMDS(x, plot = FALSE)

  df <- data.frame(
    Dim1 = mds$x,
    Dim2 = mds$y,
    group = metadata$research, 
    label = metadata$probe_name,
        research = metadata$research

  )
  
  library(ggplot2)
  
p <- ggplot(df, aes(Dim1, Dim2, color = s)) +
  geom_point(size = 4) +
  geom_text(aes(label = label), vjust = -0.5) +
  theme_minimal()

    ggsave(
    save_path,
    plot = p,
  )
}
