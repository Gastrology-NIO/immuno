if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "GSVA"))

library(DESeq2)
library(GSVA)
library(survival)


run_analyse <- function(x, dds, metadata, kegg) {

    kegg_sig <- kegg[kegg$p.adjust < 0.05, ]


   gene_sets_entrez <- lapply(
  kegg_sig$geneID,
  function(x) unique(unlist(strsplit(x, "/")))
)

names(gene_sets_entrez) <- kegg_sig$Description


all_entrez <- unique(unlist(gene_sets_entrez))

conversion <- bitr(
  all_entrez,
  fromType = "ENTREZID",
  toType = "ENSEMBL",
  OrgDb = org.Hs.eg.db
)

gene_sets <- lapply(
  gene_sets_entrez,
  function(genes) {
    conversion$ENSEMBL[
      match(genes, conversion$ENTREZID)
    ] |> 
      na.omit() |> 
      unique()
  }
)

  
  param <- gsvaParam(
      exprData = expr,
      geneSets = gene_sets
  )
  
  gsva_res <- gsva(param)
  
  dim(gsva_res)
  
  # connect from this point df with metadata
  score <- gsva_res["platelet", ]
  
  score_df <- data.frame(
      patient_id = names(platelet_score),
      score = as.numeric(score)
  )
  
  clinical <- merge(
      metadata,
      score_df,
      by.x = "patient_id",
      by.y = "patient_id"
  )
  
  
  model <- glm(
      wzrost.NtproBNP ~ score + sex + dni,
      data = clinical,
      family = binomial
  )
  
  summary(model)
  
  
  #death - zmarł, nie zmarł 1,0
  cox <- coxph(
      Surv(dni, death) ~
          score +
          sex,
      data = clinical
  )
  
  summary(cox)

  }
