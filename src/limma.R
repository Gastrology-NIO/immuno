
get_voom<-function(x, metadata){
    s<-metadata$s
    research <- metadata$research
    mm <- model.matrix(~ s + research, data=x)
    keep <- filterByExpr(x, design=mm)
    keep[c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")] <-FALSE
    d0<-x[keep,]
    y <- voom(d0 , mm, plot = F)
    return(y)
}
run_limma<-function(x, metadata, output_file){
    s<-metadata$s
    research <- metadata$research
    mm <- model.matrix(~ s + research, data=x)
    keep <- filterByExpr(x, design=mm)
    # keep[c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")] <-FALSE
    d0<-x[keep,]
    y <- voom(d0 , mm, plot = F)
    fit <- lmFit(y, mm)
    contr <- makeContrasts(researchSearched, levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp)
    top.table <- topTable(tmp, sort.by = "P", n = Inf)
    write.csv2(top.table, output_file)
    return(top.table)
}
run_deseq2<-function(x, metadata){
    s<-metadata$s
    research <- metadata$research
     mm <- model.matrix(~ s + research, data=x)
    keep <- filterByExpr(x, design=mm)
    keep[c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")] <-FALSE
    d0<-x[keep,]
    dds <- DESeqDataSetFromMatrix(countData = d0,
                                  colData = metadata,
                                  design= ~ s + research)
    diagdds = DESeq(dds , test="Wald", fitType="parametric")
    res = results(diagdds, cooksCutoff = FALSE)
        res <- results(diagdds, contrast = c("research", "Searched", "Control"))
    res<-as(res, "data.frame")
    return(res)
}



run_deseq2_paired<-function(x, metadata){
    s<-metadata$s
    research <- metadata$research
    patient_id <- metadata$patient_id
     mm <- model.matrix(~ patient_id + research, data=x)
    keep <- filterByExpr(x, design=mm)
    keep[c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")] <-FALSE
    d0<-x[keep,]
    dds <- DESeqDataSetFromMatrix(countData = d0,
                                  colData = metadata,
                                  design= ~ patient_id + research)
    diagdds = DESeq(dds , test="Wald", fitType="parametric")
    res = results(diagdds, cooksCutoff = FALSE)
        res <- results(diagdds, contrast = c("research", "Searched", "Control"))
    res<-as(res, "data.frame")
    return(res)
}



sensitivity_limma <- function(v, metadata, outlier, design_formula = ~ s + research) {

  library(limma)

  # ----------------------------
  # 1. INDEX OUTLIERA
  # ----------------------------
  keep <- metadata$probe_name != outlier
  
  v_full <- v
  v_sub  <- v[, keep]
  
  meta_full <- metadata
  meta_sub  <- metadata[keep, ]
  
  # ----------------------------
  # 2. DESIGN MATRICES
  # ----------------------------
  design_full <- model.matrix(design_formula, data = meta_full)
  design_sub  <- model.matrix(design_formula, data = meta_sub)
  
  # ----------------------------
  # 3. FIT MODELS (limma)
  # ----------------------------
  fit_full <- lmFit(v_full, design_full)
  contr <- makeContrasts(researchSearched, levels = colnames(coef(fit_full)))
  fit_full <- eBayes(contrasts.fit(fit_full, contr))
  res_full <- topTable(fit_full, number = Inf, sort.by = "none")
  
  fit_sub <- lmFit(v_sub, design_sub)
  contr <- makeContrasts(researchSearched, levels = colnames(coef(fit_sub)))
  fit_sub <- eBayes(contrasts.fit(fit_sub, contr))
  res_sub <- topTable(fit_sub, number = Inf, sort.by = "none")
  
  # ----------------------------
  # 4. COMMON GENES
  # ----------------------------
  common <- intersect(rownames(res_full), rownames(res_sub))
  
  # ----------------------------
  # 5. STABILITY METRICS
  # ----------------------------
  x <- res_full[common, "logFC"]
  y <- res_sub[common, "logFC"]
  
  ok <- is.finite(x) & is.finite(y)
  
  logfc_cor <- cor(x[ok], y[ok])
  
  sig_full <- rownames(res_full)[res_full$adj.P.Val < 0.05]
  sig_sub  <- rownames(res_sub)[res_sub$adj.P.Val < 0.05]
  
  overlap_deg <- length(intersect(sig_full, sig_sub))
  deg_loss <- length(sig_full) - overlap_deg
  
  # ----------------------------
  # 6. MDS SHIFT (INFLUENCE COMPONENT)
  # ----------------------------
  mds_full <- plotMDS(v_full, plot = FALSE)
  mds_sub  <- plotMDS(v_sub, plot = FALSE)

  mds_dist <- sqrt(
    (mean(mds_full$x) - mean(mds_sub$x))^2 +
    (mean(mds_full$y) - mean(mds_sub$y))^2
  )


  # ----------------------------
  # 7. INFLUENCE SCORE (PAPER-STYLE)
  # ----------------------------
  max_deg_loss <- max(deg_loss, 1)
logfc_diff <- abs(x - y)
    influence_score <-
  (1 - logfc_cor) +
  mean(logfc_diff[is.finite(logfc_diff)])

  # ----------------------------
  # 8. PLOT
  # ----------------------------
  pdf(paste0("sensitivity_analysis_", outlier, ".pdf"), width = 7, height = 7)

  plot(x[ok], y[ok],
       pch = 16, cex = 0.5,
       xlab = "Full model logFC",
       ylab = paste0("Without ", outlier, " logFC"),
       main = paste("Sensitivity analysis:", outlier))

  abline(0,1,col="red")

  dev.off()

  # ----------------------------
  # 9. OUTPUT
  # ----------------------------
  cat("\n============================\n")
  cat("Sensitivity analysis:", outlier, "\n")
  cat("============================\n")

  cat("Genes full:", nrow(res_full), "\n")
  cat("Genes reduced:", nrow(res_sub), "\n")
  cat("Correlation logFC:", round(logfc_cor, 3), "\n")
  cat("DEG overlap:", overlap_deg, "\n")
  cat("DEG loss:", deg_loss, "\n")
  cat("MDS shift:", round(mds_dist, 3), "\n")
  cat("INFLUENCE SCORE:", round(influence_score, 3), "\n")

  # ----------------------------
  # 10. RETURN
  # ----------------------------
  return(list(
    full = res_full,
    reduced = res_sub,
    correlation_logFC = logfc_cor,
    overlap_DEG = overlap_deg,
    deg_loss = deg_loss,
    mds_shift = mds_dist,
    influence_score = influence_score
  ))
}



test<-function(metadata, y){
results_list <- list()

    for (probe in metadata$probe_name) {

  cat("Processing:", probe, "\n")

  res <- sensitivity_limma(y, metadata, probe)

  results_list[[probe]] <- data.frame(
    sample = probe,

    logFC_cor = if (!is.null(res$correlation_logFC)) res$correlation_logFC else NA,
    DEG_overlap = if (!is.null(res$overlap_DEG)) res$overlap_DEG else NA,
    deg_loss = if (!is.null(res$deg_loss)) res$deg_loss else NA,

    influence_score = if (!is.null(res$influence_score)) res$influence_score else NA
  )
}

results <- do.call(rbind, results_list)

return(results)
}

