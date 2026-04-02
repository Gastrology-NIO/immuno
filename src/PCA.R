


plotPCA<function(x, save_path){
  mds <- plotMDS(x, plot = FALSE)

  df <- data.frame(
    Dim1 = mds$x,
    Dim2 = mds$y,
    group = research
  )
  
  library(ggplot2)
  
  p<-ggplot(df, aes(Dim1, Dim2, color = group)) +
    geom_point(size = 4) +
    theme_minimal()

    ggsave(
    save_path,
    plot = p,
  )
}
