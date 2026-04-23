
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
  # 3. FIT MODELS
  # ----------------------------
  fit_full <- lmFit(v_full, design_full)
  fit_full <- eBayes(fit_full)
  res_full <- topTable(fit_full, number = Inf)

  fit_sub <- lmFit(v_sub, design_sub)
  fit_sub <- eBayes(fit_sub)
  res_sub <- topTable(fit_sub, number = Inf)

  # ----------------------------
  # 4. COMMON GENES
  # ----------------------------
  common <- intersect(rownames(res_full), rownames(res_sub))

  # ----------------------------
  # 5. STABILITY METRICS
  # ----------------------------
  logfc_cor <- cor(
    res_full[common, "logFC"],
    res_sub[common, "logFC"]
  )

  sig_full <- rownames(res_full)[res_full$adj.P.Val < 0.05]
  sig_sub  <- rownames(res_sub)[res_sub$adj.P.Val < 0.05]

  overlap_deg <- length(intersect(sig_full, sig_sub))

  # ----------------------------
  # 6. PLOT
  # ----------------------------
  plot(res_full[common, "logFC"],
       res_sub[common, "logFC"],
       pch = 16, cex = 0.5,
       xlab = "Full model logFC",
       ylab = paste0("Without ", outlier, " logFC"),
       main = paste("Sensitivity analysis:", outlier))

  abline(0,1,col="red")

  # ----------------------------
  # 7. OUTPUT SUMMARY
  # ----------------------------
  cat("\n============================\n")
  cat("Sensitivity analysis for:", outlier, "\n")
  cat("============================\n")

  cat("Genes in full model:", nrow(res_full), "\n")
  cat("Genes in reduced model:", nrow(res_sub), "\n")
  cat("Correlation logFC:", round(logfc_cor, 3), "\n")
  cat("Overlap DEGs:", overlap_deg, "\n")
  cat("Unique to full:", length(setdiff(sig_full, sig_sub)), "\n")
  cat("Unique to reduced:", length(setdiff(sig_sub, sig_full)), "\n")

  # ----------------------------
  # 8. RETURN RESULTS
  # ----------------------------
  return(list(
    full = res_full,
    reduced = res_sub,
    correlation_logFC = logfc_cor,
    overlap_DEG = overlap_deg
  ))
}
