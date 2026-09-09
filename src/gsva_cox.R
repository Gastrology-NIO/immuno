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
vsd <- vst(dds, blind = TRUE)
  expr <- assay(vsd)


  param <- gsvaParam(
      exprData = expr,
      geneSets = gene_sets
  )
  
  gsva_res <- gsva(param)
  
  metadata$patient_id <- as.character(metadata$probe_name)

metadata <- metadata[
    match(colnames(gsva_res), metadata$patient_id),
]

rownames(metadata) <- NULL
metadata$wzrost.NtproBNP[
    metadata$wzrost.NtproBNP == ""
] <- NA
logistic_results <- data.frame()

for (pathway in rownames(gsva_res)) {
    
    clinical <- metadata
    
    clinical$score <- as.numeric(
        gsva_res[pathway, clinical$patient_id]
    )
    clinical$wzrost.NtproBNP <- as.numeric(clinical$wzrost.NtproBNP)

clinical$s <- factor(clinical$s)

clinical$score <- as.numeric(clinical$score)

clinical$dni <- as.numeric(clinical$dni)
    
    model <- glm(
        wzrost.NtproBNP ~ score + s + dni,
        data = clinical,
        family = binomial
    )
    
    coef_model <- summary(model)$coefficients
    
    logistic_results <- rbind(
        logistic_results,
        data.frame(
            pathway = pathway,
            beta = coef_model["score", "Estimate"],
            SE = coef_model["score", "Std. Error"],
            pvalue = coef_model["score", "Pr(>|z|)"],
            OR = exp(coef_model["score", "Estimate"]),
            CI_low = exp(
                coef_model["score", "Estimate"] -
                1.96 * coef_model["score", "Std. Error"]
            ),
            CI_high = exp(
                coef_model["score", "Estimate"] +
                1.96 * coef_model["score", "Std. Error"]
            )
        )
    )
}
    library(survival)

cox_results <- data.frame()

for (pathway in rownames(gsva_res)) {
    
    clinical <- metadata
    
    clinical$score <- as.numeric(
        gsva_res[pathway, clinical$patient_id]
    )
    
    model <- coxph(
        Surv(dni, death) ~ score + s,
        data = clinical
    )
    
    coef_model <- summary(model)$coefficients
    
    cox_results <- rbind(
        cox_results,
        data.frame(
            pathway = pathway,
            beta = coef_model["score", "coef"],
            HR = coef_model["score", "exp(coef)"],
            SE = coef_model["score", "se(coef)"],
            pvalue = coef_model["score", "Pr(>|z|)"]
        )
    )
}

cox_results <- data.frame()

for (pathway in rownames(gsva_res)) {
    
    clinical <- metadata
    
    clinical$score <- as.numeric(
        gsva_res[pathway, clinical$patient_id]
    )
    
    model <- coxph(
        Surv(dni, zgon_bool) ~ score + s,
        data = clinical
    )
    
    s <- summary(model)
    
    cox_results <- rbind(
        cox_results,
        data.frame(
            pathway = pathway,
            HR = s$coefficients["score", "exp(coef)"],
            CI_low = s$conf.int["score", "lower .95"],
            CI_high = s$conf.int["score", "upper .95"],
            pvalue = s$coefficients["score", "Pr(>|z|)"]
        )
    )
}

logistic_results$FDR <- p.adjust(
    logistic_results$pvalue,
    method = "BH"
)

cox_results$FDR <- p.adjust(
    cox_results$pvalue,
    method = "BH"
)
    

  }
