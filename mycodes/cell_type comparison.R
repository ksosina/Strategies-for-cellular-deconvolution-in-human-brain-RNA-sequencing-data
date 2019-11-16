# Replicate plots

# Preamble ----------------------------------------------------------------

# Load libraries
packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer", "matrixStats")
loaded_libs <- sapply(packs, library, character.only = T)




# Plot ---------------------------------------------------------------

# dt <- fread(file.path(".", "MatchedCellComp", "singleCell_iPSC_quake_coefEsts_calibration_Zscale_adultOnly.csv"))

## load data
load("./NAc_Nicotine_hg38_rseGene_rawCounts_postQCSamples_n223.rda")
load("./processed_objects.rda")
colnames(bVals) <- rownames(pd) <- pd$BrNum

## filter for gene expression
bVals <- bVals[,rse_gene$BrNum]
pd <- pd[rse_gene$BrNum,]

## get DNAm composition
load("./cellComp_estimates_cellLines_NeuNs.rda")
coefEsts_DNAm <- coefs[rownames(coefs) %in% rownames(bVals),]
cellCounts_DNAm <- minfi:::projectCellType(bVals[rownames(coefEsts_DNAm),],coefEsts_DNAm)
colnames(cellCounts_DNAm) <- paste0(colnames(cellCounts_DNAm) , "_DNAm")
pd_DNAm <- cbind(pd, cellCounts_DNAm)

## do PCA of DNAm data
ooMeth <- order(rowSds(bVals), decreasing=TRUE)[1:50000]
pca_DNAm <- prcomp(t(bVals[ooMeth,]))
# pcaVars_DNAm <- getPcaVars(pca_DNAm)
plot(pca_DNAm$x[,1], pd_DNAm$NeuN_neg_DNAm) # very correlated
pd_DNAm$PC1_cellType <- pca_DNAm$x[,1]

## do the RNA-based estiamtion
# load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/deconvolution/cell_type/singleCell_iPSC_quake_coefEsts_calibration_Zscale_adultOnly.rda")
load("./MatchedCellComp/singleCell_iPSC_quake_coefEsts_calibration_Zscale_adultOnly.rda")
yExprs <- log2(getRPKM(rse_gene, "Length")+1)

# project
yExprs_Z <- scale(yExprs[rownames(coefEsts),])
cellCounts_Exprs <- minfi:::projectCellType(yExprs_Z,coefEsts)
colnames(cellCounts_Exprs) <- paste0(colnames(cellCounts_Exprs) , "_RNA")
colData(rse_gene) <- cbind(colData(rse_gene), as.data.frame(cellCounts_Exprs))
rse_gene$NeuN_neg_RNA <- 1 - rse_gene$Neurons_RNA 

# colData(rse_gene) <- cbind(colData(rse_gene), pd_DNAm[,43:45])
colData(rse_gene) <- cbind(colData(rse_gene), pd_DNAm[,43:44])
save(rse_gene, file = "NAc_rse_gene_withCompEsts.rda")

###################
## check plots ####
load("NAc_rse_gene_withCompEsts.rda")
yExprs <- log2(getRPKM(rse_gene, "Length")+1)

pdf("check_comp_estimates_DNAmVsExprs_1.pdf")
palette(brewer.pal(4,"Dark2"))
# par(mar=c(5,6,3,2), cex.axis=1.8,cex.lab=1.8,cex.main=1.6, mfrow = c(2,2), pty = "m")
par(mfrow = c(2,2), pty = "m")
plot(rse_gene$NeuN_neg_RNA , rse_gene$NeuN_neg_DNAm,
     xlab = "RNA-based", ylab = "DNAm-based", 
     main = "Cell Composition Estimation",
     pch = 21, bg = 1, ylim = c(0,1),xlim=c(0,1))
abline(0,1,lty=2,lwd=2)
points(rse_gene$Neurons_RNA , rse_gene$NeuN_pos_DNAm,
       pch = 21, bg = 2)
legend("topleft", c("NeuN-", "NeuN+"), pch = 15, col=1:2)

plot(rse_gene$NeuN_neg_DNAm , rse_gene$NeuN_neg_RNA,
     ylab = "RNA-based", xlab = "DNAm-based", 
     main = "Cell Composition Estimation",
     pch = 21, bg = 1, ylim = c(0,1),xlim=c(0,1))
abline(0,1,lty=2,lwd=2)
points(rse_gene$Neurons_RNA , rse_gene$NeuN_pos_DNAm,
       pch = 21, bg = 2)
legend("topleft", c("NeuN-", "NeuN+"), pch = 15, col=1:2)


plot(density((rse_gene$NeuN_neg_RNA)), 
     main = "NeuN- cell Composition Estimation",
     pch = 21, bg = 1, col = "#1B9E77", ylim = c(0,5))
lines(density((rse_gene$NeuN_neg_DNAm)), col = "#D95F02")
legend("topleft", c("NeuN- RNA", "NeuN- DNAm"), pch = 15, col=1:2)

plot(density((rse_gene$Neurons_RNA)), 
     main = "NeuN+ cell Composition Estimation",
     pch = 21, bg = 1, col = "#1B9E77", ylim = c(0,6))
lines(density((rse_gene$NeuN_pos_DNAm)), col = "#D95F02")
legend("topright", c("NeuN+ RNA", "NeuN+ DNAm"), pch = 15, col=1:2)

dev.off()

cor(rse_gene$NeuN_neg_RNA , rse_gene$NeuN_neg_DNAm)

## not so much quality
cor(rse_gene$NeuN_neg_RNA , rse_gene$totalAssignedGene)
plot(rse_gene$NeuN_neg_RNA , rse_gene$totalAssignedGene,
     xlab = "RNA-Non-Neuronal Proportion Estimation", 
     ylab = "Gene Mapping Rate", 
     pch = 21, bg = "grey")
abline(0,1,col="red")


cell_props_rnavsdnam <- data.table(NeuN_neg_RNA = rse_gene$NeuN_neg_RNA , NeuN_neg_DNAm = rse_gene$NeuN_neg_DNAm,
           Neurons_RNA = rse_gene$Neurons_RNA , NeuN_pos_DNAm = rse_gene$NeuN_pos_DNAm)

cell_props_rnavsdnam[, `:=`(total_rna = NeuN_neg_RNA + Neurons_RNA,
                            total_dnam = NeuN_neg_DNAm + NeuN_pos_DNAm)]

summary(cell_props_rnavsdnam[, .(total_rna, total_dnam)])
