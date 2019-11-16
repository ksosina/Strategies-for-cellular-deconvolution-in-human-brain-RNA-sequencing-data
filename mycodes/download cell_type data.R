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




# Merge/Save --------------------------------------------------------------

##########
## merge #

## merge w/ scorecard
yExprs_Merge <- cbind(yExprsYeo, yExprsQuake)	
group <- c(as.character(rse_geneYeo$cell_type), 
           as.character(rse_geneQuake$Cell_type))

group <- factor(group,levels =c("iPSC", "NPC","Fetal_replicating",
                                "Fetal_quiescent", "OPC", "Neurons", "Astrocytes",
                                "Oligodendrocytes", "Microglia", "Endothelial"))
group2 <- group == "Neurons"

group2 <- data.table(cells = colnames(yExprs_Merge), Neurons = group2, all = group)

fwrite(group2, file.path(getwd(), "cell_type_comp", "cell_type_info.txt"))
fwrite(data.table(genes = row.names(yExprs_Merge),yExprs_Merge), file.path(getwd(), "cell_type_comp", "rna_gene_x_cells.txt"))
