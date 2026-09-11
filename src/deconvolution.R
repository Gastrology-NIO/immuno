if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("immunedeconv")

library(immunedeconv)

res_quantiseq <- deconvolute(
    expr,
    method = "quantiseq"
)

library(tidyr)
library(dplyr)

cell_fraction <- res_quantiseq %>%
    pivot_wider(
        names_from = cell_type,
        values_from = fraction
    )


clinical <- clinical %>%
    left_join(
        cell_fraction,
        by = "sample_id"
    )


all(colnames(expr) %in% clinical$sample_id)


wilcox.test(
    `Monocytes`,
    data = subset(clinical, research %in% c("Control", "Searched"))
)



library(dplyr)
library(tidyr)
library(purrr)

cell_cols <- setdiff(
    colnames(cell_fraction),
    "sample_id"
)

cell_stats <- map_dfr(cell_cols, function(cell) {

    x <- clinical %>%
        filter(research %in% c("Control", "Searched"))

    wt <- wilcox.test(
        x[[cell]] ~ x$research
    )

    data.frame(
        cell_type = cell,
        pvalue = wt$p.value
    )
})

cell_stats$FDR <- p.adjust(
    cell_stats$pvalue,
    method = "BH"
)

cell_stats


cor.test(
    clinical$Platelet_activation,
    clinical$Neutrophils,
    method = "spearman"
)



cor_results <- expand.grid(
    pathway = c(
        "Platelet_activation",
        "Coagulation"
    ),
    cell = cell_cols,
    stringsAsFactors = FALSE
)

cor_results <- cor_results %>%
    rowwise() %>%
    mutate(
        rho = cor.test(
            clinical[[pathway]],
            clinical[[cell]],
            method = "spearman"
        )$estimate,
        pvalue = cor.test(
            clinical[[pathway]],
            clinical[[cell]],
            method = "spearman"
        )$p.value
    ) %>%
    ungroup()

cor_results$FDR <- p.adjust(
    cor_results$pvalue,
    method = "BH"
)


# coxph(
#     Surv(dni, zgon_bool) ~
#         score +
#         s +
#         wzrost.NtproBNP +
#         Neutrophils +
#         Monocytes +
#         CD8_T_cells,
#     data = clinical
# )
# coxph(
#     Surv(dni, zgon_bool) ~
#         score +
#         s +
#         wzrost.NtproBNP +
#         Neutrophils,
#     data = clinical
# )
