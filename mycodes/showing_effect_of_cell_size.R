# Preamble ----------------------------------------------------------------

packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer", 
           "mixtools","matrixStats", "MuSiC", "convert", "xbioc", "ggplot2", "sva", "plotly")
libs_loaded <- sapply(packs, library, character.only = T)

# MuSic -------------------------------------------------------------------

rna_gene_x_cells <- fread(file.path(".", "cell_data", "rna_gene_x_cells.txt"))
cell_type_data <- fread(file.path(".", "cell_data", "cell_type_info.txt"))
cell_type_data <- cell_type_data[ all %in% c("Neurons", "OPC",  "Astrocytes",
                                             "Oligodendrocytes", "Microglia", "Endothelial"), ]


# # Subset to informative genes
# my_vals <- est_ref_params(rna_gene_x_cells, cell_type_data)
# my_params <- my_vals$out
# # my_params <- my_params[genes %in% pv[prop.v_e > 0.1, genes], ]
# params_out_m <- my_vals$mean_rna
# # params_out_m <- params_out_m[genes %in% pv[prop.v_e > 0.2, genes], ]
# params_out_l <- my_vals$lambda
# # params_out_l <- params_out_l[genes %in% pv[prop.v_e > 0.2, genes], ]

load("./NAc_Nicotine_hg38_rseGene_rawCounts_postQCSamples_n223.rda")
load("./MatchedCellComp/singleCell_iPSC_quake_coefEsts_calibration_Zscale_adultOnly.rda")
yExprs <- log2(getRPKM(rse_gene, "Length")+1)


# project
yExprs <- yExprs[rownames(coefEsts),]
yExprs_Z <- scale(yExprs[rownames(coefEsts),])


my_params <- my_params[genes %in% rownames(coefEsts),]

gg <- rownames(coefEsts)
# gg2 <- gg[gg %in% params_out_m$genes]
# gg <- gg[gg %in% my_params$genes]
# 
# my_params <- my_params[match(gg, genes),]
# params_out_m <- params_out_m[match(gg2, genes),]
# params_out_l <- params_out_l[match(gg2, genes),]


analysis_dt <- which(row.names(yExprs) %in% gg)
analysis_dt <- yExprs[analysis_dt, ]


sc_data <- rna_gene_x_cells[genes %in% row.names(analysis_dt),]
sc_data <- sc_data[match(row.names(analysis_dt), genes),]
rn <- sc_data$genes
cn <- colnames(sc_data)
cn <- cn[cn %in% cell_type_data$cells]
sc_data <- sc_data[ , cn, with = F]
sc_data <- data.matrix(sc_data)
row.names(sc_data) <- rn
sc_data <- ExpressionSet(assayData = sc_data)
bulk_data <- ExpressionSet(assayData = yExprs)
music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells)
music_est_nnls <- data.table(samples  = row.names(music_est$Est.prop.allgene), music_est$Est.prop.allgene)
music_est <- data.table(samples  = row.names(music_est$Est.prop.allgene),music_est$Est.prop.weighted) 

music_est_nnls[, RNum := sapply(strsplit(samples, "_"), "[", 1)]
music_est[, RNum := sapply(strsplit(samples, "_"), "[", 1)]
load("./NAc_rse_gene_withCompEsts.rda")
music_est_nnls <- music_est_nnls[match(rse_gene$RNum, RNum),]
music_est <- music_est[match(rse_gene$RNum, RNum),]

# Estimate S_K ------------------------------------------------------------

rna_gene_x_cells <- data.table::fread(file.path(".", "cell_data", "rna_gene_x_cells.txt"))
cell_type_data <- data.table::fread(file.path(".", "cell_data", "cell_type_info.txt"))
cell_type_data <- cell_type_data[ all %in% c("Neurons", "OPC",  "Astrocytes",
                                             "Oligodendrocytes", "Microglia", "Endothelial"), ]

types_cell <- unique(cell_type_data$all)
my_s_k <- sapply(types_cell, function(g){
  c_g <- cell_type_data$cells[cell_type_data$all == g]
  print(length(c_g))
  anal_data <- rna_gene_x_cells[genes %in% row.names(analysis_dt), c_g, with = F]
  anal_out <- colSums(anal_data)
  s <- mean(anal_out)
  return(s)
})


sc.markers <- row.names(bulk_data)
music_basis_ests <- music_basis(sc_data, non.zero = TRUE, markers = sc.markers, 
                                clusters = cell_type_data$all, samples = cell_type_data$cells, select.ct = NULL, 
                                ct.cov = FALSE, verbose = TRUE)



# CIBERSORT ---------------------------------------------------------------

# source(file.path(".", "cell_data", "cibersort_data", "CIBERSORT.R"))
# results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN, absolute, abs_method)


# osmfish loom ------------------------------------------------------------

# Install package
my_packs <- installed.packages(fields = "Package")

if(sum(my_packs[,1] %in% "loomR") != 1){
  devtools::install_github(repo = "hhoeflin/hdf5r")
  devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
}{
  library(loomR)
}
# Connect to the loom file in read/write mode
lfile <- connect(filename = file.path("cell_data", "osmfish data", "osmFISH_SScortex_mouse_all_cells.loom"), mode = "r+")
lfile
lfile$col.attrs$size_um2[1:5]

# Pull three bits of metadata from the column attributes
my_cell_atts <- data.table("CellID" = lfile$col.attrs$CellID[],
          "ClusterName" = lfile$col.attrs$ClusterName[],
          "CellArea" = lfile$col.attrs$size_um2[],
          "valid" = lfile$col.attrs$Valid[],
          "totalmolecules" = lfile$col.attrs$Total_molecules[])

cell_types_aj <- fread("celltypes_ajEdits.csv")
cell_types_aj <- cell_types_aj[!is.na(Class)]

osmfihs_ests <- my_cell_atts[, .(mean(CellArea), .N), by = .(ClusterName, valid)]

osmfihs_ests <- inner_join(osmfihs_ests, cell_types_aj[, .(ClusterName, Celltype = Class, Subclass)]) %>% data.table

# osmfihs_ests[, Celltype := sapply(strsplit(ClusterName, " "), "[", 1)]
# osmfihs_ests[, Celltype := ifelse(Celltype %in% c("Inhibitory", "Pyramidal", "pyramidal"), "Neuron", Celltype)] 

my_ests <- osmfihs_ests[, .("w_mean" = weighted.mean(V1, N)), by = Celltype]
my_ests[, rel_size := w_mean/sum(w_mean)]

# fwrite(unique(osmfihs_ests[, .(ClusterName, Celltype)]), "celltypes.txt")

inner_join(osmfihs_ests[, .(ClusterName, "osmFISH_est" = V1, Celltype)], 
           data.table("music_est" = my_s_k, Celltype = c("Astrocyte", "Endothelial", "Microglia", "Neuron", "Oligodendrocyte", "OPC")))

pd_data <- inner_join(my_ests[, .(Celltype, "osmFISH_est" = w_mean)], 
           data.table("music_est" = my_s_k, Celltype = c("Astrocyte", "Endothelial", "Microglia", "Neuron", "Oligodendrocyte", "OPC"))) %>% data.table

pd_data_unique <- pd_data[, .(osmFISH_est, Celltype, music_est)] %>% unique

plot(pd_data$osmFISH_est, pd_data$music_est, xlab = "osmFISH", ylab = "MUSiC")


plot_dat_oshm <- my_cell_atts
plot_dat_oshm[, Celltype := NULL]
# plot_dat_oshm[, Celltype := sapply(strsplit(ClusterName, " "), "[", 1)]
plot_dat_oshm <- inner_join(plot_dat_oshm, cell_types_aj[, .(ClusterName, Celltype = Class, Subclass)]) %>% data.table


rna_ests <- plot_dat_oshm[, .(mean(totalmolecules), .N, Celltype = unique(Celltype)), by = .(ClusterName, valid)]
rna_ests <- rna_ests[, .("w_mean" = weighted.mean(V1, N)), by = Celltype]
rna_ests[, rel_abun := w_mean/sum(w_mean)]

w_mole_count <- rna_ests$w_mean/my_ests$w_mean
w_mole_count/sum(w_mole_count)

out <- inner_join(rna_ests[, .(Celltype, "Rel_abun RNA in mouse" = rel_abun, total_mole = w_mean)], my_ests[, .(Celltype, "Rel_size in mouse" = rel_size, size = w_mean)]) %>% 
  inner_join(data.table("music_est" = my_s_k/sum(my_s_k), Celltype = c("Astrocyte", "Endothelial", "Microglia", "Neuron", "Oligodendrocyte", "OPC"))) %>% 
  inner_join(data.table(Celltype = rna_ests$Celltype, scaled_rna = w_mole_count/sum(w_mole_count))) %>% data.table

out[, Celltype := c("Neurons", "Astrocytes", "OPC",  "Oligodendrocytes", "Microglia", "Endothelial") ]
colMeans(music_basis_ests$S, na.rm = T)
out[match(names(colMeans(music_basis_ests$S, na.rm = T)), out$Celltype),]



{
  
  
  my_ests <- osmfihs_ests
  my_ests[, Neurons := Celltype == "Neuron"]
  my_ests <- my_ests[, .("w_mean" = weighted.mean(V1, N)), by = Neurons]
  my_ests[, rel_size := w_mean/sum(w_mean)]
  
  
  plot_dat_oshm <- my_cell_atts
  plot_dat_oshm[, Celltype := NULL]
  plot_dat_oshm <- inner_join(plot_dat_oshm, cell_types_aj[, .(ClusterName, Celltype = Class, Subclass)]) %>% data.table
  plot_dat_oshm[, Neurons := Celltype == "Neuron"]
  
  rna_ests <- plot_dat_oshm[, .(mean(totalmolecules), .N, Celltype = unique(Celltype), Neurons = unique(Neurons)), by = .(ClusterName, valid)]
  rna_ests <- rna_ests[, .("w_mean" = weighted.mean(V1, N)), by = Neurons]
  rna_ests[, rel_abun := w_mean/sum(w_mean)]
  
  w_mole_count <- rna_ests$w_mean/my_ests$w_mean
  w_mole_count/sum(w_mole_count)
  
  my_out <- inner_join(rna_ests[, .(Neurons, "Rel_abun RNA in mouse" = rel_abun, total_mole = w_mean)], my_ests[, .(Neurons, "Rel_size in mouse" = rel_size, size = w_mean)]) %>% 
    inner_join(data.table(Neurons = rna_ests$Neurons, scaled_rna = w_mole_count/sum(w_mole_count))) %>% data.table
  
 
  my_music_est <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
                                my_abundance = my_out[, .(Celltype = Neurons, w_mean = total_mole )])
  my_music_est_nnls_mole <- data.table(samples  = row.names(my_music_est$Est.prop.allgene), my_music_est$Est.prop.allgene)
  my_music_est_mole <- data.table(samples  = row.names(my_music_est$Est.prop.allgene),my_music_est$Est.prop.weighted) 
  
  
  my_music_est <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
                                my_abundance = my_out[, .(Celltype = Neurons, w_mean = size )])
  my_music_est_nnls_size <- data.table(samples  = row.names(my_music_est$Est.prop.allgene), my_music_est$Est.prop.allgene)
  my_music_est_size <- data.table(samples  = row.names(my_music_est$Est.prop.allgene),my_music_est$Est.prop.weighted) 
}



ggplot(plot_dat_oshm, aes(y = log10(totalmolecules), x = ClusterName, fill = `ClusterName` )) + 
  geom_boxplot() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, size = 8, hjust = .9, face = "bold"))

ggplot(plot_dat_oshm, aes(y = log10(totalmolecules), x = log10(CellArea), color = `ClusterName` )) + 
  geom_point() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.direction = "horizontal")
  
ggsave(file.path(".", "figs", "molecule_vs_cell_type.png"))

lm(scale(totalmolecules) ~ scale(CellArea) + as.factor(Celltype), data = plot_dat_oshm[Celltype != "Excluded",]) %>% summary
glm(round(totalmolecules, 0) ~  scale(CellArea) + as.factor(Celltype), 
    data = plot_dat_oshm[Celltype != "Excluded",], family = "quasipoisson") %>% summary

# Cell type data September 9th --------------------------------------------

# Load cell type data (9th September 2019)
# Gene x Cell matrix
count_mat <- fread(file.path("cell_type_data", "countMatrix_n4169-NAc-nuclei.csv"))
cell_mat <- fread(file.path("cell_type_data", "pd-cellTypeAssignment_n4169.csv"))

names(count_mat)[1] <- "Gene"
names(cell_mat)[1] <- "Cells"

# Check
identical(names(count_mat)[-1], cell_mat$Cells) # Ordering is also the same

# nCount_RNA per Cell is the same as the column sum of count_mat for each Cell
# There are 24,189 Genes

my_out_cell_mat <- cell_mat[, .(size = mean(nCount_RNA)), by = nucleusCellType]


my_out_cell_mat$CellType <- sapply(my_out_cell_mat$nucleusCellType, switch, 
                                   "Oligo" = "Oligodendrocytes", "Micro" = "Microglia", 
                                   "Astro" = "Astrocytes", "Neuron" = "Neurons", "OPC" = "OPC")



rna_gene_x_cells <- fread(file.path(".", "cell_data", "rna_gene_x_cells.txt"))
cell_type_data <- fread(file.path(".", "cell_data", "cell_type_info.txt"))
cell_type_data <- cell_type_data[ all %in% c("Neurons", "OPC",  "Astrocytes",
                                             "Oligodendrocytes", "Microglia"), ]


# # Subset to informative genes

load("./NAc_Nicotine_hg38_rseGene_rawCounts_postQCSamples_n223.rda")
load("./MatchedCellComp/singleCell_iPSC_quake_coefEsts_calibration_Zscale_adultOnly.rda")
yExprs <- log2(getRPKM(rse_gene, "Length")+1)


# project
yExprs <- yExprs[rownames(coefEsts),]
yExprs_Z <- scale(yExprs[rownames(coefEsts),])


my_params <- my_params[genes %in% rownames(coefEsts),]

gg <- rownames(coefEsts)



analysis_dt <- which(row.names(yExprs) %in% gg)
analysis_dt <- yExprs[analysis_dt, ]


sc_data <- rna_gene_x_cells[genes %in% row.names(analysis_dt),]
sc_data <- sc_data[match(row.names(analysis_dt), genes),]
rn <- sc_data$genes
cn <- colnames(sc_data)
cn <- cn[cn %in% cell_type_data$cells]
sc_data <- sc_data[ , cn, with = F]
sc_data <- data.matrix(sc_data)
row.names(sc_data) <- rn
sc_data <- ExpressionSet(assayData = sc_data)
bulk_data <- ExpressionSet(assayData = yExprs)
music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells)
music_est_nnls <- data.table(samples  = row.names(music_est$Est.prop.allgene), music_est$Est.prop.allgene)
music_est <- data.table(samples  = row.names(music_est$Est.prop.allgene),music_est$Est.prop.weighted) 

music_est_nnls[, RNum := sapply(strsplit(samples, "_"), "[", 1)]
music_est[, RNum := sapply(strsplit(samples, "_"), "[", 1)]
load("./NAc_rse_gene_withCompEsts.rda")
music_est_nnls <- music_est_nnls[match(rse_gene$RNum, RNum),]
music_est <- music_est[match(rse_gene$RNum, RNum),]






my_music_est_nw <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells,
                              my_abundance = my_out_cell_mat[, .(Celltype = CellType, w_mean = size )])

my_music_est_nw <- data.table(samples  = row.names(my_music_est$Est.prop.allgene), my_music_est_nw$Est.prop.weighted) 

round(cor(music_est$`Neurons`, rse_gene$NeuN_pos_DNAm)^2, 4)
round(cor(my_music_est_nw$`Neurons`, rse_gene$NeuN_pos_DNAm)^2, 4)

# Campbell ----------------------------------------------------------------


campbell <- readRDS(file.path("cell_data", "Campbell", "campbell.rds"))

# Abbas data --------------------------------------------------------------

mixes_abbas <- fread(file.path(".", "cell_data", "abbas_et_al_data", "GSE11103_matrix_mixtures.txt"))
cell_type_data <- fread(file.path(".", "cell_data", "abbas_et_al_data", "GSE11103_matrix_pure.txt"))
eval_mix <- function(dt){
  if(!is.data.table(dt)){
    dt <- data.table("ID" = sapply(strsplit(row.names(dt), " "), "[", 1),
                     "Mix" = sapply(strsplit(row.names(dt), " "), "[", 2),
                     dt)
  }
  # truth
  MixA <- c("Jurkat" = 2.5, "IM-9" = 1.25, "Raji" =  2.5,"THP-1" = 3.75)/10
  MixB <-  c("Jurkat" = 0.5, "IM-9" = 3.17, "Raji" = 4.75, "THP-1" = 1.58)/10
  MixC <-  c("Jurkat" = 0.1, "IM-9" = 4.95, "Raji" = 1.65, "THP-1" = 3.3)/10
  MixD  <- c("Jurkat" = 0.02, "IM-9" = 3.33, "Raji" = 3.33, "THP-1" = 3.33)/10
  Mix <- rbind(MixA, MixB, MixC, MixD)
  Mix <- data.table("Mix" = paste0("Mix", LETTERS[1:4]), Mix)
  my_set <-which( names(Mix) %in% names(dt))
  Mix <- Mix[, my_set, with = F]
  
  pear <- sapply(1:nrow(dt), function(i){
    dt <- dt[i,]
    m <- dt$Mix
    dt <- dt[, -c(1:2), with = F]
    dt <- unlist(dt)
    t <- Mix[Mix == m, -1, with = F]
    t <- unlist(t)
    r <- cor(dt, t)
    return(r)
  })
  
  pear_rmse <- sapply(1:nrow(dt), function(i){
    dt <- dt[i,]
    m <- dt$Mix
    dt <- dt[, -c(1:2), with = F]
    dt <- unlist(dt)
    t <- Mix[Mix == m, -1, with = F]
    t <- unlist(t)
    r <- crossprod(dt, t)/length(t)
    r <- c(r)
    return(r)
  })
  
  pear_ratio <- sapply(1:nrow(dt), function(i){
    dt <- dt[i,]
    m <- dt$Mix
    dt <- dt[, -c(1:2), with = F]
    dt <- unlist(dt)
    t <- Mix[Mix == m, -1, with = F]
    t <- unlist(t)
    r <- dt/t
    r <- mean(r)
    r <- c(r)
    return(r)
  })
  
  pear_ratio <- tapply(pear_ratio, rep(paste0("Mix", LETTERS[1:4]), each = 3), mean)
  pear_rmse <- tapply(pear_rmse, rep(paste0("Mix", LETTERS[1:4]), each = 3), mean)
  pear <- tapply(pear, rep(paste0("Mix", LETTERS[1:4]), each = 3), mean)
  pear <- data.table(type = names(pear), pearson = (pear), rmse = pear_rmse, mean_ratio = pear_ratio)
  return(pear)
}

# Normalize
{
  yExprs_normfactors <- edgeR::calcNormFactors(cbind(data.matrix(mixes_abbas[, -1]), cell_type_data[, -c(1), with = F]), method = "TMM" )
  a <- ncol(mixes_abbas)  - 1 + 1
  b <- a + ncol(cell_type_data[, -c(1), with = F]) - 1
  q_out <- yExprs_normfactors[a:b]
  q_out <- replicate(nrow(mixes_abbas)  , q_out)
  q_out <- t(q_out)
  raw_dt <- cbind(cell_type_data[, 1], cell_type_data[, -c(1), with = F]*q_out)
  
  yExprs_normfactors <- yExprs_normfactors[-c(a:b)]
  yExprs_normfactors <- replicate(nrow(mixes_abbas)  , yExprs_normfactors)
  yExprs_normfactors <- t(yExprs_normfactors)
  all_dt <- cbind(mixes_abbas[, 1], mixes_abbas[,-1]*yExprs_normfactors)
}

# subest to informative genes
# p <- .75
# cutoff <- quantile(cell_type_data_var$prop, probs = c(1-p,p))
# info_idx <- which(cell_type_data_var$prop >= cutoff[2] & cell_type_data_var$`E(Var(Y|X))` > 0)

# Analyze
sc_data <- raw_dt
# sc_data <- raw_dt[info_idx, ]
rn <- sc_data$`!Sample_title`
cn <- colnames(sc_data)[-1]
# cn <- cn[cn %in% cell_type_data$cells]
# sc_data <- sc_data[ , cn, with = F]
sc_data <- data.matrix(sc_data[, -1])
row.names(sc_data) <- rn
music_dt <- data.matrix(all_dt[, -1])
# music_dt <- data.matrix(all_dt[info_idx, -1])
row.names(music_dt) <- rn
sc_data <- ExpressionSet(assayData = sc_data)
bulk_data <- ExpressionSet(assayData = music_dt)
music_est_abbas <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = rep(c( "Jurkat", "IM-9", "Raji", "THP-1"), each = 3), samples = cn, normalize = T)
music_est_nnls_abbas <- music_est_abbas$Est.prop.allgene
music_est_abbas <- music_est_abbas$Est.prop.weighted


eval_mix(music_est_nnls_abbas)
eval_mix(music_est_abbas)


music_est_nnls_abbas <- data.table("ID" = sapply(strsplit(row.names(music_est_nnls_abbas), " "), "[", 1),
                             "Mix" = sapply(strsplit(row.names(music_est_nnls_abbas), " "), "[", 2),
                             music_est_nnls_abbas)
music_est_abbas <- data.table("ID" = sapply(strsplit(row.names(music_est_abbas), " "), "[", 1),
                        "Mix" = sapply(strsplit(row.names(music_est_abbas), " "), "[", 2),
                        music_est_abbas)
# truth
MixA <- c("Jurkat" = 2.5, "IM-9" = 1.25, "Raji" =  2.5,"THP-1" = 3.75)/10
MixB <-  c("Jurkat" = 0.5, "IM-9" = 3.17, "Raji" = 4.75, "THP-1" = 1.58)/10
MixC <-  c("Jurkat" = 0.1, "IM-9" = 4.95, "Raji" = 1.65, "THP-1" = 3.3)/10
MixD  <- c("Jurkat" = 0.02, "IM-9" = 3.33, "Raji" = 3.33, "THP-1" = 3.33)/10

Mix <- rbind(MixA, MixB, MixC, MixD)
Mix <- data.table("Mix" = paste0("Mix", LETTERS[1:4]), Mix)

plot_dt_music <- rbind(inner_join(music_est_abbas[, .(Mix, Estimate = Jurkat)],  Mix[, .(Mix, Truth = Jurkat, type = "Jurkat")]),
                     inner_join(music_est_abbas[, .(Mix, Estimate = `IM-9`)],  Mix[, .(Mix, Truth = `IM-9`, type = "IM-9")]),
                     inner_join(music_est_abbas[, .(Mix, Estimate = Raji)],  Mix[, .(Mix, Truth = Raji, type = "Raji")]),
                     inner_join(music_est_abbas[, .(Mix, Estimate = `THP-1`)],  Mix[, .(Mix, Truth = `THP-1`, type = "THP-1")])) %>% data.table


plot_dt_nnls <- rbind(inner_join(music_est_nnls_abbas[, .(Mix, Estimate = Jurkat)],  Mix[, .(Mix, Truth = Jurkat, type = "Jurkat")]),
                       inner_join(music_est_nnls_abbas[, .(Mix, Estimate = `IM-9`)],  Mix[, .(Mix, Truth = `IM-9`, type = "IM-9")]),
                       inner_join(music_est_nnls_abbas[, .(Mix, Estimate = Raji)],  Mix[, .(Mix, Truth = Raji, type = "Raji")]),
                       inner_join(music_est_nnls_abbas[, .(Mix, Estimate = `THP-1`)],  Mix[, .(Mix, Truth = `THP-1`, type = "THP-1")])) %>% data.table
with(plot_dt_music, cor(x = Truth, y = Estimate))
with(plot_dt_nnls, cor(x = Truth, y = Estimate))




cibersort <- fread(file.path(".", "cell_data", "cibersort_data","CIBERSORT.Output_Job2.txt"))

cibersort <- data.table("ID" = sapply(strsplit(cibersort$`Input Sample`, " "), "[", 1),
                        "Mix" = sapply(strsplit(cibersort$`Input Sample`, " "), "[", 2),
                        cibersort[,-1])

cibersort_e <- data.matrix(cibersort[, .(Jurkat, `IM-9`, `Raji`, `THP-1`)])
rownames(cibersort_e) <- rownames(music_est_nnls_abbas)

eval_mix(cibersort_e)

plot_dt_cibersort <- rbind(inner_join(cibersort[, .(Mix, Estimate = Jurkat)],  Mix[, .(Mix, Truth = Jurkat, type = "Jurkat")]),
                           inner_join(cibersort[, .(Mix, Estimate = `IM-9`)],  Mix[, .(Mix, Truth = `IM-9`, type = "IM-9")]),
                           inner_join(cibersort[, .(Mix, Estimate = Raji)],  Mix[, .(Mix, Truth = Raji, type = "Raji")]),
                           inner_join(cibersort[, .(Mix, Estimate = `THP-1`)],  Mix[, .(Mix, Truth = `THP-1`, type = "THP-1")]))

# Plot --------------------------------------------------------------------

plot_dat <- rbind(data.table("MUSiC" = music_est$`FALSE`, "type" = "NeuN-"),
                  data.table("MUSiC" = music_est$`TRUE`, "type" = "NeuN+"))

plot_dat2 <- rbind(data.table("Houseman DNAm-based" = rse_gene$NeuN_neg_DNAm),
                   data.table("Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm))

plot_dat3 <- rbind(data.table("Houseman RNA-based" = rse_gene$NeuN_neg_RNA),
                   data.table("Houseman RNA-based" = rse_gene$Neurons_RNA))

plot_dat <- cbind(plot_dat, plot_dat2, plot_dat3)




transparent_legend =  theme(
  legend.background = element_rect(fill ="transparent"),
  legend.key = element_rect(fill = "transparent",
                            color = "transparent")
)

remove_grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))

no_x_axis_label <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())


cols <- RColorBrewer::brewer.pal(length(unique(plot_dat$type)), "Dark2")
cols <- cols[1:2]
names(cols) <- unique(plot_dat$type)

c1 <- round(cor(music_est$`FALSE`, rse_gene$NeuN_neg_DNAm)^2, 2)
p1 <- ggplot(data = plot_dat, aes(x = MUSiC, y = `Houseman DNAm-based`, color = type)) + 
  geom_point(size = 3) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = cols, name = "Cell type" ) +
  annotate("text", x = 0.2, y=.9, label = paste0("italic(R) ^ 2 == ", c1), parse = T) +
  transparent_legend + remove_grid +
  ggtitle("Relationship between Houseman Methylation and MUSiC estimates") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.9,0.35),
        legend.title = element_text(face="bold"))
ggplotly(p1)

c2 <- round(cor(rse_gene$NeuN_neg_RNA, rse_gene$NeuN_neg_DNAm)^2, 2)
p2 <- ggplot(data = plot_dat, aes(x = `Houseman RNA-based`, y = `Houseman DNAm-based`, color = type)) + 
  geom_point(size = 3) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = cols, name = "Cell type" ) +
  annotate("text", x = 0, y=.9, label = paste0("italic(R) ^ 2 == ", c2), parse = T) +
  transparent_legend + remove_grid +
  ggtitle("Relationship between Houseman Methylation and  RNA estimates") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.9,0.35),
        legend.title = element_text(face="bold"))


library(SCnorm)
rna_gene_x_cells <- data.table::fread(file.path(".", "cell_data", "rna_gene_x_cells.txt"))
cell_type_data <- data.table::fread(file.path(".", "cell_data", "cell_type_info.txt"))
cell_type_data <- cell_type_data[ all %in% c("Neurons", "OPC",  "Astrocytes",
                                             "Oligodendrocytes", "Microglia", "Endothelial"), ]

scnorm_data <- data.matrix(rna_gene_x_cells[, -1])
row.names(scnorm_data) <- rna_gene_x_cells$genes
plotCountDepth(scnorm_data)


library(scran)
scale_factors_by_cells <- computeSumFactors(scnorm_data)
