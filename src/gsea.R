
BiocManager::install('enrichplot', character.only = TRUE)
install.packages('ggtree')
library(clusterProfiler)
library(enrichplot)
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
original_gene_list <- top.table$logFC

# name the vector
names(original_gene_list) <- rownames(top.table)

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
c = sort(gene_list, decreasing = TRUE)


gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.1, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "fdr")
gse[gse@result$p.adjust<0.1]->gse_tmp
df<-as.data.frame(gse)
write.csv2(df, "gseGO-01.csv")
# !!simplify

df<-as.data.frame(gse)
write.csv2(df, "gseGO.csv")
require(DOSE)
D<-dotplot(gse, showCategory=10, split=c(".sign")) + facet_grid(~.sign) +
  theme(
    axis.text.y = element_text(size = 8)  # Change y-axis label size
  )

ggsave(
  "./dotplot_GeneRatio.svg",
  plot = D,
)


D<-dotplot(gse, showCategory=10, split=c(".sign")) + facet_grid(~.sign) +
  theme(
    axis.text.y = element_text(size = 8)  # Change y-axis label size
  )

ggsave(
  "./dotplot_GeneRatio.svg",
  plot = D,
)


edo <- pairwise_termsim(gse)
E<-enrichplot::emapplot(edo, showCategory = 20, layout = "kk")
ggsave(
  "./emapplot.svg",
  plot = E,
)

cnet<-cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
ggsave(
  "./cnetplot2.png",
  plot = cnet,
  width=10, height=8, dpi=300
)
