if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("immunedeconv")

library(immunedeconv)
library(EnsDb.Hsapiens.v86)
library(ensembldb)
library(tidyr)
library(dplyr)

gene_info <- genes(
    EnsDb.Hsapiens.v86,
    columns = c("gene_id", "gene_name")
)

gene_info$gene_length <- width(gene_info)

counts <- dge$counts

idx <- match(rownames(counts), gene_info$gene_id)

gene_length <- gene_info$gene_length[idx]
gene_symbol <- gene_info$gene_name[idx]

sum(is.na(gene_length))
sum(is.na(gene_symbol))

keep <- !is.na(gene_length) &
        !is.na(gene_symbol) &
        gene_symbol != ""

counts2 <- counts[keep, , drop = FALSE]
gene_length2 <- gene_length[keep]
gene_symbol2 <- gene_symbol[keep]

rpk <- counts2 / (gene_length2 / 1000)

tpm <- sweep(
    rpk,
    2,
    colSums(rpk) / 1e6,
    "/"
)

colSums(tpm)[1:6]

tpm_symbol <- rowsum(
    tpm,
    group = gene_symbol2,
    reorder = FALSE
)

res_quantiseq <- immunedeconv::deconvolute(
    tpm_symbol,
    method = "quantiseq"
)


# cell_fraction <- res_quantiseq %>%
#     pivot_wider(
#         names_from = cell_type,
#         values_from = fraction
#     )


# clinical <- clinical %>%
#     left_join(
#         cell_fraction,
#         by = "sample_id"
#     )


# all(colnames(expr) %in% clinical$sample_id)


# wilcox.test(
#     `Monocytes`,
#     data = subset(clinical, research %in% c("Control", "Searched"))
# )



# library(dplyr)
# library(tidyr)
# library(purrr)

# cell_cols <- setdiff(
#     colnames(cell_fraction),
#     "sample_id"
# )

# cell_stats <- map_dfr(cell_cols, function(cell) {

#     x <- clinical %>%
#         filter(research %in% c("Control", "Searched"))

#     wt <- wilcox.test(
#         x[[cell]] ~ x$research
#     )

#     data.frame(
#         cell_type = cell,
#         pvalue = wt$p.value
#     )
# })

# cell_stats$FDR <- p.adjust(
#     cell_stats$pvalue,
#     method = "BH"
# )

# cell_stats


# cor.test(
#     clinical$Platelet_activation,
#     clinical$Neutrophils,
#     method = "spearman"
# )



# cor_results <- expand.grid(
#     pathway = c(
#         "Platelet_activation",
#         "Coagulation"
#     ),
#     cell = cell_cols,
#     stringsAsFactors = FALSE
# )

# cor_results <- cor_results %>%
#     rowwise() %>%
#     mutate(
#         rho = cor.test(
#             clinical[[pathway]],
#             clinical[[cell]],
#             method = "spearman"
#         )$estimate,
#         pvalue = cor.test(
#             clinical[[pathway]],
#             clinical[[cell]],
#             method = "spearman"
#         )$p.value
#     ) %>%
#     ungroup()

# cor_results$FDR <- p.adjust(
#     cor_results$pvalue,
#     method = "BH"
# )


