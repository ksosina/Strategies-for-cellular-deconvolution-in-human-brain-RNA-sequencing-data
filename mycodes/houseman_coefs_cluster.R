
# Connect -----------------------------------------------------------------

# Open ssh and enter passwd
jhpce

# Start interactive session
qrsh

# Open R
R



# Preamble ----------------------------------------------------------------

# Load libraries
packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer", "rafalib")
sapply(packs, library, character.only = T)


pal <- brewer.pal(10, "Set3")
pal[2] <- "#8B795E" # darken NPCs

MAINDIR <- "/dcl01/lieber/ajaffe/lab/libd_stem_timecourse"





# Load Files --------------------------------------------------------------

# load gene-level
load("/dcl01/ajaffe/data/lab/singleCell/Darmanis/rna-seq-pipeline_run2/rse_gene_Darmanis_scRNASeq_Darmanis_n466.Rdata")
rse_geneQuake <- rse_gene
rm(rse_gene)

# Darmanis et al data
# phenotype
load("/dcl01/lieber/ajaffe/PublicData/Brain/Darmanis_Quake/rdas/phenotype_quake_processed.rda")
colData(rse_geneQuake) <- cbind(colData(rse_geneQuake), 
                               pd[match(colnames(rse_geneQuake), pd$RunName),2:14])

rse_geneQuake$Cell_type[rse_geneQuake$Cell_type == "hybrids"] <- "hybrid"
## exclude hybrid
rse_geneQuake <- rse_geneQuake[,rse_geneQuake$Cell_type != "hybrid"]

# get expression
yExprsQuake <- log2(getRPKM(rse_geneQuake, "Length")+1)


# Song et al -> Gore et al data
## bring in YEO iPSC data
load(file.path(MAINDIR,"yeo_singleCell/rse_gene_yeo_n214.Rdata"))

rse_geneYeo <- rse_gene[,rse_gene$cell_type %in% c("iPSC", "NPC") &
                         rse_gene$sample_type=="Single_Cell"]
yExprsYeo <- log2(getRPKM(rse_geneYeo, "Length")+1)

##########
## merge #

## merge w/ scorecard
yExprs_Merge <- cbind(yExprsYeo, yExprsQuake)	
group <- c(as.character(rse_geneYeo$cell_type), 
          as.character(rse_geneQuake$Cell_type))

group <- factor(group,levels =c("iPSC", "NPC","Fetal_replicating",
                               "Fetal_quiescent", "OPC", "Neurons", "Astrocytes",
                               "Oligodendrocytes", "Microglia", "Endothelial"))

## #split by age cat and region					
tIndexes <- splitit(group)

tstatList <- lapply(tIndexes, function(i) {
  x <- rep(0, ncol(yExprs_Merge))
  x[i] <- 1
  return(genefilter::rowttests(yExprs_Merge, factor(x)))
})

numProbes=25
probeList <- lapply(tstatList, function(x) {
  y <- x[which(x[, "p.value"] < 1e-15), ]
  yUp <- y[order(y[, "dm"], decreasing = FALSE),] # signs are swapped
  rownames(yUp)[1:numProbes]
})

## filter
trainingProbes <- unique(unlist(probeList))
trainingProbes = trainingProbes[!is.na(trainingProbes)]

mergeMarkerExprs <- yExprs_Merge[trainingProbes, ]

mergeMarkerMeanExprs <- colMeans(mergeMarkerExprs)

form <- as.formula(sprintf("y ~ %s - 1", paste(levels(group),collapse = "+")))
phenoDF <- as.data.frame(model.matrix(~group - 1))
colnames(phenoDF) <- sub("^group", "", colnames(phenoDF))

## try z-score
mergeMarkerExprsZ = scale(mergeMarkerExprs)
mergeMarkerMeanExprsZ = colMeans(mergeMarkerExprsZ)

## do calibration
coefEsts <- minfi:::validationCellType(Y = mergeMarkerExprsZ, 
                                       pheno = phenoDF, modelFix = form)$coefEsts

save(coefEsts, mergeMarkerMeanExprsZ, 
     file = "singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda")
write.csv(coefEsts, file = "singleCell_iPSC_quake_coefEsts_calibration_Zscale.csv")
