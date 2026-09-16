if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("immunedeconv")

library(immunedeconv)
library(EnsDb.Hsapiens.v86)
library(ensembldb)
library(tidyr)
library(dplyr)
library(ggplot2)
library(rstatix)



############################

deconvolution_difference<-function(name, metadata){
    gene_info <- genes(
            EnsDb.Hsapiens.v86,
            columns = c("gene_id", "gene_name")
        )
        gene_info$gene_length <- width(gene_info)
        metadata<- metadata[metadata$probe_name != "",]
        
        x<-load_DGE(metadata,  "./htseq/") 
        dge <- calcNormFactors(x)
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
        

        
        # metadata tylko dla próbek, które są w Quantiseq
        meta <- metadata %>%
          select(probe_name, research) %>%
          filter(!is.na(research))
        
        # Quantiseq: wide -> long
        quant_long <- res_quantiseq %>%
          pivot_longer(
            cols = -cell_type,
            names_to = "probe_name",
            values_to = "fraction"
          ) %>%
          left_join(meta, by = "probe_name") %>%
          filter(!is.na(research))
        
        head(quant_long)
        
        results <- quant_long %>%
          group_by(cell_type) %>%
          wilcox_test(
            fraction ~ research
          ) %>%
          adjust_pvalue(method = "BH") %>%
          add_significance("p.adj")

        rbc_results <- quant_long %>%
          group_by(cell_type) %>%
          wilcox_effsize(
            fraction ~ research,
            paired = FALSE
          ) %>%
          select(cell_type, effsize)

        
        summary_wide <- summary_cells %>%
          select(cell_type, research, median) %>%
          pivot_wider(
            names_from = research,
            values_from = median
          ) %>%
          mutate(
            difference = Searched - Control
          )
        
        summary_wide
    
            

        # dodaj p.adj do danych wykresu
        results <- results %>%
              left_join(rbc_results, by = "cell_type") %>%
              rename(rank_biserial = effsize)
        quant_long_plot <- quant_long %>%
          left_join(
            results %>% 
              select(cell_type, p.adj, rank_biserial),
            by = "cell_type"
          )
        gg <- ggplot(
          quant_long_plot,
          aes(x = research, y = fraction, fill = research)
        ) +
          geom_boxplot(
            width = 0.6,
            outlier.shape = NA
          ) +
          geom_jitter(
            width = 0.12,
            size = 2,
            alpha = 0.8
          ) +
          facet_wrap(~ cell_type, scales = "free_y") +
          
          geom_text(
          data = quant_long_plot %>%
            group_by(cell_type) %>%
            summarise(
              p.adj = first(p.adj),
              rank_biserial = first(rank_biserial),
              .groups = "drop"
            ) %>%
            filter(!is.na(p.adj)),
          aes(
            x = 1.5,
            y = Inf,
            label = paste0(
              "adj. p = ",
              format.pval(p.adj, digits = 2, eps = 1e-4),
              "\nRBC = ",
              round(rank_biserial, 2)
            )
          ),
          inherit.aes = FALSE,
          vjust = 1.5,
          size = 3.5
        ) +
          
          theme_bw() +
          labs(
            x = NULL,
            y = "Cell fraction",
            fill = "Research"
          ) +
          theme(
            strip.text = element_text(face = "bold"),
            legend.position = "top"
          )
        
        ggsave(
          paste0("./",name,"_doconv.png"),
          plot = gg,
          width = 12,
          height = 8,
          dpi = 300
        )
    }
#######################


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


