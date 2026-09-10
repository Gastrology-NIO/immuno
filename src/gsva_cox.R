if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "GSVA"))

library(DESeq2)
library(GSVA)
library(survival)

library(httr)
library(dplyr)
#As metadata, we should use all metadata, not only edge cases
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
        gsva_res[pathway, clinical$probe_name]
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
    
#     pathway       HR    CI_low  CI_high      pvalue        FDR
# 1  hsa04518 1.350671 1.0028659 1.819098 0.047838196 0.06909962
# 2  hsa04611 1.414180 1.0577650 1.890690 0.019335579 0.04189375
# 3  hsa04610 1.582251 1.1602771 2.157690 0.003740213 0.03507437
# 4  hsa04510 1.442771 1.0841247 1.920064 0.011940306 0.03507437
# 5  hsa05414 1.414371 1.0102204 1.980206 0.043467655 0.06909962
# 6  hsa04820 1.241436 0.9115343 1.690735 0.169987576 0.18415321
# 7  hsa05410 1.432295 1.0219356 2.007433 0.036983675 0.06868397
# 8  hsa04640 1.160365 0.8519573 1.580415 0.345397742 0.34539774
# 9  hsa04512 1.284622 0.9335886 1.767647 0.124051193 0.14660595
# 10 hsa04814 1.298775 0.9657504 1.746639 0.083734432 0.10885476
# 11 hsa04810 1.438994 1.0781184 1.920665 0.013490143 0.03507437
# 12 hsa05132 1.463849 1.0975174 1.952455 0.009509472 0.03507437
# 13 hsa05131 1.495538 1.1220164 1.993407 0.006047341 0.03507437

  }
