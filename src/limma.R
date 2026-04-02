

run_limma<-function(x, metadata, output_file){
    s<-metadata$s
    research <- metadata$research
    mm <- model.matrix(~ s + research, data=x)
    keep <- filterByExpr(x, design=mm)
    # keep[c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")] <-FALSE
    d0<-x[keep,]
    y <- voom(d0 , mm, plot = F)
    fit <- lmFit(y, mm)
    contr <- makeContrasts(researchContol, levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp)
    top.table <- topTable(tmp, sort.by = "P", n = Inf)
    write.csv2(top.table, output_file)
    return(top.table)
}
