if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "GSVA"))

library(DESeq2)
library(GSVA)
library(survival)


run_analyse <- function(x, dds, metadata) {

  rownames(metadata) <- metadata$patient_id
  metadata <- metadata[colnames(x), , drop = FALSE]
  metadata$s <- factor(metadata$s)
  metadata$research <- factor(metadata$research)
  vsd <- vst(dds, blind = TRUE)
  expr <- assay(vsd)
  
  # gene_sets <- list(
  #     platelet = platelet_genes,
  #     coagulation = coagulation_genes
  # )
  genes_set<-list()
  
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
