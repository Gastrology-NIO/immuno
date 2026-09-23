if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "GSVA"))

library(DESeq2)
library(GSVA)
library(survival)

library(httr)
library(dplyr)


# log-CPM


#As metadata, we should use all metadata, not only edge cases
run_cox <- function(kegg, name) {

    
    metadata<-load_data_meta("./data/metadata_jag.csv")
    conditions <- list(type='checkpoint 1')
    metadata_checkpoint1 <- cut_metadata(metadata, conditions)
    x<-load_DGE(metadata_checkpoint1,  "./htseq/") 
    dge <- calcNormFactors(x)

    expr <- cpm(
        dge,
        log = TRUE,
        prior.count = 1
    )
    rownames(expr) <- sub("\\..*$", "", rownames(expr))
    colnames(expr) <- basename(colnames(expr))
    
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
    clinical$zmiana.NTproBNP_10 <- as.numeric(clinical$zmiana.NTproBNP)
    
    # Cox model
    model <- coxph(
        Surv(dni, zgon_bool) ~ score +s +zmiana.NTproBNP_10,
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
    write.csv(cox_results, "cox.csv")
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



GSVA to GO to metadata_0 checkpoint1



metadata<-load_data_meta("./data/metadata_jag.csv")
conditions <- list(type='checkpoint 1')
metadata_checkpoint1 <- cut_metadata(metadata, conditions)
x<-load_DGE(metadata_checkpoint1,  "./htseq/") 
dge <- calcNormFactors(x)

expr <- cpm(
dge,
log = TRUE,
prior.count = 1
)
rownames(expr) <- sub("\\..*$", "", rownames(expr))
colnames(expr) <- basename(colnames(expr))
gsva_mat<-c()
go_sig <- go[go$p.adjust < 0.05, ]
gene_sets <- list()

# pobieramy poprawne mapowanie Entrez -> GO
go2gene <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
  keytype = "ENTREZID",
  columns = c("ENTREZID", "GO")
)

# usuwamy puste przypisania
go2gene <- go2gene[
  !is.na(go2gene$GO) &
  !is.na(go2gene$ENTREZID),
]

for (go_id in go_sig$ID) {

  print(go_id)

  # GO -> Entrez
  entrez_genes <- unique(
    go2gene$ENTREZID[
      go2gene$GO == go_id
    ]
  )

  # jeśli GO nie ma w bazie
  if (length(entrez_genes) == 0) {
    warning(
      paste("Brak genów dla", go_id)
    )
    gene_sets[[go_id]] <- character(0)
    next
  }

  # Entrez -> Ensembl
  ensembl_genes <- mapIds(
    org.Hs.eg.db,
    keys = entrez_genes,
    keytype = "ENTREZID",
    column = "ENSEMBL",
    multiVals = "first"
  )

  ensembl_genes <- unique(
    na.omit(ensembl_genes)
  )

  # tylko geny obecne w macierzy ekspresji
  ensembl_genes <- intersect(
    ensembl_genes,
    rownames(expr)
  )

  gene_sets[[go_id]] <- ensembl_genes
}

# rozmiary gene setów po dopasowaniu do RNA-seq
gene_set_sizes <- sapply(
  gene_sets,
  length
)


summary(gene_set_sizes)

# opcjonalnie usuwamy bardzo małe zestawy
gene_sets <- gene_sets[
  gene_set_sizes >= 10
]

# GSVA
param <- gsvaParam(
  exprData = expr,
  geneSets = gene_sets
)

gsva_res <- gsva(param)
gsva_mat <- gsva_res
rownames(metadata_checkpoint1) <- metadata_checkpoint1$probe_name

# dopasowanie próbek GSVA do metadata
common_samples <- intersect(
  colnames(gsva_mat),
  rownames(metadata_checkpoint1)
)

gsva_mat <- gsva_mat[, common_samples, drop = FALSE]
meta_cor <- metadata_checkpoint1[common_samples, , drop = FALSE]

# sprawdzenie kolejności
identical(colnames(gsva_mat), rownames(meta_cor))



cor_results <- data.frame()

for (go_id in rownames(gsva_mat)) {

  x <- as.numeric(meta_cor$zmiana.NTproBNP)
  y <- gsva_mat[go_id, ]

  keep <- complete.cases(x, y)

  if (sum(keep) >= 3) {

    test <- cor.test(
      x[keep],
      y[keep],
      method = "spearman"
    )

    cor_results <- rbind(
      cor_results,
      data.frame(
        GO = go_id,
        rho = unname(test$estimate),
        p_value = test$p.value,
        n = sum(keep)
      )
    )
  }
}

cor_results$padj <- p.adjust(
  cor_results$p_value,
  method = "BH"
)

cor_results <- cor_results[
  order(cor_results$padj, -abs(cor_results$rho)),
]




GSVA to GO to metadata_5 checkpoint2



metadata<-load_data_meta("./data/metadata_jag.csv")
conditions <- list(type='checkpoint 2')
metadata_checkpoint1 <- cut_metadata(metadata, conditions)
metadata_checkpoint1<- metadata_checkpoint1[metadata_checkpoint1$probe_name != "",]
metadata_checkpoint1<-metadata_checkpoint1[!c(metadata_checkpoint1$probe_name %in% c("63IM")),]
x<-load_DGE(metadata_checkpoint1,  "./htseq2/") 
dge <- calcNormFactors(x)

expr <- cpm(
dge,
log = TRUE,
prior.count = 1
)
rownames(expr) <- sub("\\..*$", "", rownames(expr))
colnames(expr) <- basename(colnames(expr))
gsva_mat<-c()
go<-enrich_MF
go_sig <- go[go$p.adjust < 0.05, ]
gene_sets <- list()

# pobieramy poprawne mapowanie Entrez -> GO
go2gene <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
  keytype = "ENTREZID",
  columns = c("ENTREZID", "GO")
)

# usuwamy puste przypisania
go2gene <- go2gene[
  !is.na(go2gene$GO) &
  !is.na(go2gene$ENTREZID),
]

for (go_id in go_sig$ID) {

  print(go_id)

  # GO -> Entrez
  entrez_genes <- unique(
    go2gene$ENTREZID[
      go2gene$GO == go_id
    ]
  )

  # jeśli GO nie ma w bazie
  if (length(entrez_genes) == 0) {
    warning(
      paste("Brak genów dla", go_id)
    )
    gene_sets[[go_id]] <- character(0)
    next
  }

  # Entrez -> Ensembl
  ensembl_genes <- mapIds(
    org.Hs.eg.db,
    keys = entrez_genes,
    keytype = "ENTREZID",
    column = "ENSEMBL",
    multiVals = "first"
  )

  ensembl_genes <- unique(
    na.omit(ensembl_genes)
  )

  # tylko geny obecne w macierzy ekspresji
  ensembl_genes <- intersect(
    ensembl_genes,
    rownames(expr)
  )

  gene_sets[[go_id]] <- ensembl_genes
}

# rozmiary gene setów po dopasowaniu do RNA-seq
gene_set_sizes <- sapply(
  gene_sets,
  length
)


summary(gene_set_sizes)

# opcjonalnie usuwamy bardzo małe zestawy
gene_sets <- gene_sets[
  gene_set_sizes >= 10
]

# GSVA
param <- gsvaParam(
  exprData = expr,
  geneSets = gene_sets
)

gsva_res <- gsva(param)
gsva_mat <- gsva_res
rownames(metadata_checkpoint1) <- metadata_checkpoint1$probe_name

# dopasowanie próbek GSVA do metadata
common_samples <- intersect(
  colnames(gsva_mat),
  rownames(metadata_checkpoint1)
)

gsva_mat <- gsva_mat[, common_samples, drop = FALSE]
meta_cor <- metadata_checkpoint1[common_samples, , drop = FALSE]

# sprawdzenie kolejności
identical(colnames(gsva_mat), rownames(meta_cor))



cor_results <- data.frame()

for (go_id in rownames(gsva_mat)) {

  x <- as.numeric(meta_cor$zmiana.NTproBNP)
  y <- gsva_mat[go_id, ]

  keep <- complete.cases(x, y)

  if (sum(keep) >= 3) {

    test <- cor.test(
      x[keep],
      y[keep],
      method = "spearman"
    )

    cor_results <- rbind(
      cor_results,
      data.frame(
        GO = go_id,
        rho = unname(test$estimate),
        p_value = test$p.value,
        n = sum(keep)
      )
    )
  }
}

cor_results$padj <- p.adjust(
  cor_results$p_value,
  method = "BH"
)

cor_results <- cor_results[
  order(cor_results$padj, -abs(cor_results$rho)),
]
write.csv2(cor_results, "corr_MF_Cardiac_biomarker_defined.csv")
