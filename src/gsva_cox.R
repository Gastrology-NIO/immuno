if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "GSVA"))

library(DESeq2)
library(GSVA)
library(survival)

library(httr)
library(dplyr)

run_analyse <- function(x, dds, metadata, kegg) {

    kegg_sig <- kegg[kegg$p.adjust < 0.05, ]


    for (pathway_id in rownames(kegg_sig)) {
        
        print(pathway_id)
        
        # pobranie genów przypisanych do pathway
        url <- paste0(
            "https://rest.kegg.jp/link/hsa/",
            pathway_id
        )
        
        kegg_path <- read.delim(
            url,
            header = FALSE,
            sep = "\t",
            stringsAsFactors = FALSE
        )
        
        colnames(kegg_path) <- c("pathway", "gene")
        
        # usunięcie prefiksu hsa:
        entrez_genes <- sub("^hsa:", "", kegg_path$gene)
        
        # Entrez -> Ensembl
        ensembl_genes <- mapIds(
            org.Hs.eg.db,
            keys = entrez_genes,
            keytype = "ENTREZID",
            column = "ENSEMBL",
            multiVals = "first"
        )
        
        ensembl_genes <- unique(na.omit(ensembl_genes))
        
        # zostawiamy tylko geny obecne w macierzy ekspresji
        ensembl_genes <- intersect(
            ensembl_genes,
            rownames(expr)
        )
        
        # zapisujemy gene set
        gene_sets[[pathway_id]] <- ensembl_genes
    }

gene_set_sizes <- sapply(gene_sets, length)

param <- gsvaParam(
    exprData = expr,
    geneSets = gene_sets
)

gsva_res <- gsva(param)

    
rownames(metadata) <- NULL
metadata$wzrost.NtproBNP[
    metadata$wzrost.NtproBNP == ""
] <- NA
logistic_results <- data.frame()

cox_results <- data.frame()

for (pathway in rownames(gsva_res)) {
    
    clinical <- metadata
    
    # GSVA score
    clinical$score <- as.numeric(
        gsva_res[pathway, clinical$patient_id]
    )
    
    # standaryzacja do 1 SD
    clinical$score <- as.numeric(
        scale(clinical$score)
    )
    
    # Cox model
    model <- coxph(
        Surv(dni, zgon_bool) ~ score +s,
        data = clinical
    )
    
    model_summary <- summary(model)
    
    cox_results <- rbind(
        cox_results,
        data.frame(
            pathway = pathway,
            HR = model_summary$coefficients[
                "score", "exp(coef)"
            ],
            CI_low = model_summary$conf.int[
                "score", "lower .95"
            ],
            CI_high = model_summary$conf.int[
                "score", "upper .95"
            ],
            pvalue = model_summary$coefficients[
                "score", "Pr(>|z|)"
            ]
        )
    )
}


cox_results$FDR <- p.adjust(
    cox_results$pvalue,
    method = "BH"
)
    

  }
