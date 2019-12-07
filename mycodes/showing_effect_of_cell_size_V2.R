# Preamble ----------------------------------------------------------------

packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer", 
           "mixtools","matrixStats", "MuSiC", "convert", "xbioc", "ggplot2", "sva", "plotly",
           "doParallel", "parallel")
libs_loaded <- sapply(packs, library, character.only = T)

type_anal <- "Neurons" # or Celltype

scrna <- "single" # both

# Get common set of Genes -------------------------------------------------------------------

rna_gene_x_cells <- fread(file.path(".", "cell_data", "rna_gene_x_cells.txt"))
cell_type_data <- fread(file.path(".", "cell_data", "cell_type_info.txt"))
cell_type_data <- cell_type_data[ all %in% c("Neurons", "OPC",  "Astrocytes",
                                             "Oligodendrocytes", "Microglia"), ]


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
# yExprs <- yExprs[rownames(coefEsts),]
# yExprs_Z <- scale(yExprs[rownames(coefEsts),])


# gg <- rownames(coefEsts)

gg <- rownames(yExprs)
gg <- gg[gg%in%rna_gene_x_cells$genes]

# Load cell type data (9th September 2019)
# Gene x Cell matrix
count_mat <- fread(file.path("cell_type_data", "countMatrix_n4169-NAc-nuclei.csv"))
cell_mat <- fread(file.path("cell_type_data", "pd-cellTypeAssignment_n4169.csv"))
cell_mat[, Neurons:= nucleusCellType == "Neuron"]
names(count_mat)[1] <- "Gene"
names(cell_mat)[1] <- "Cells"

cell_mat_dis_25 <- fread(file.path("~", "cell_type_data", "clusterMarkers_25-per-cluster_n4169.csv"))
cell_mat_dis_50 <- fread(file.path("~", "cell_type_data", "clusterMarkers_50-per-cluster_n4169.csv"))

count_mat_50 <- count_mat[Gene %in% cell_mat_dis_50$gene]
count_mat_25 <- count_mat[Gene %in% cell_mat_dis_25$gene]





names(count_mat_50)[1] <- "Gene"
names(count_mat_25)[1] <- "Gene"

# Check
identical(names(count_mat)[-1], cell_mat$Cells) # Ordering is also the same
identical(names(count_mat_50)[-1], cell_mat$Cells) # Ordering is also the same
identical(names(count_mat_25)[-1], cell_mat$Cells) # Ordering is also the same

# nCount_RNA per Cell is the same as the column sum of count_mat for each Cell
# There are 24,189 Genes

if(type_anal == "Neurons"){
  my_out_cell_mat <- cell_mat[, .(size = mean(nCount_RNA)), by = Neurons]
  
  my_out_cell_mat_50 <- data.table(Cells = names(count_mat_50)[-1],
                                   nCount_RNA = apply(count_mat_50[, -1], 2, sum))
  my_out_cell_mat_50 <- inner_join(my_out_cell_mat_50, cell_mat[, .(Cells, Neurons)]) %>% data.table
  my_out_cell_mat_50 <- my_out_cell_mat_50[, .(size = mean(nCount_RNA)), by = Neurons]
  
  my_out_cell_mat_25 <- data.table(Cells = names(count_mat_25)[-1],
                                   nCount_RNA = apply(count_mat_25[, -1], 2, sum))
  my_out_cell_mat_25 <- inner_join(my_out_cell_mat_25, cell_mat[, .(Cells, Neurons)]) %>% data.table
  my_out_cell_mat_25 <- my_out_cell_mat_25[, .(size = mean(nCount_RNA)), by = Neurons]
}else{
  my_out_cell_mat <- cell_mat[, .(size = mean(nCount_RNA)), by = nucleusCellType]
  my_out_cell_mat$CellType <- sapply(my_out_cell_mat$nucleusCellType, switch, 
                                     "Oligo" = "Oligodendrocytes", "Micro" = "Microglia", 
                                     "Astro" = "Astrocytes", "Neuron" = "Neurons", "OPC" = "OPC")
  
  cell_mat[, CellType := sapply(nucleusCellType, switch, 
                                "Oligo" = "Oligodendrocytes", "Micro" = "Microglia", 
                                "Astro" = "Astrocytes", "Neuron" = "Neurons", "OPC" = "OPC") ]
  
  
  my_out_cell_mat_50 <- data.table(Cells = names(count_mat_50)[-1],
                                   nCount_RNA = apply(count_mat_50[, -1], 2, sum))
  my_out_cell_mat_50 <- inner_join(my_out_cell_mat_50, cell_mat[, .(Cells, CellType)]) %>% data.table
  my_out_cell_mat_50 <- my_out_cell_mat_50[, .(size = mean(nCount_RNA)), by = CellType]
  
  my_out_cell_mat_25 <- data.table(Cells = names(count_mat_25)[-1],
                                   nCount_RNA = apply(count_mat_25[, -1], 2, sum))
  my_out_cell_mat_25 <- inner_join(my_out_cell_mat_25, cell_mat[, .(Cells, CellType)]) %>% data.table
  my_out_cell_mat_25 <- my_out_cell_mat_25[, .(size = mean(nCount_RNA)), by = CellType]
}








# Change HGNC name to ensembl IDs

id_dt <- as.data.table(rowData(rse_gene))
id_dt <- id_dt[, .(gencodeID, Gene = Symbol)]

# gencodeGTF <- rtracklayer::import(con= file.path("~", "cell_type_data", "gencode.v25.annotationGRCh38.gtf"), format="gtf")
# gencodeGENES <- mcols(gencodeGTF)[which(gencodeGTF$type=="gene"),c("gene_id","type","gene_type","gene_name")]
gencodeGTF <- rtracklayer::import(con= file.path("~", "cell_type_data", "genes.gtf"), format="gtf")
gencodeGENES <- mcols(gencodeGTF)[which(gencodeGTF$type=="gene"),c("gene_id","gene_version","gene_name")]
# rownames(gencodeGENES) <- gencodeGENES$gene_id

id_dt <- as.data.table(gencodeGENES)
id_dt <- id_dt[, .(gene_id, Gene = gene_name)]
id_dt <- inner_join(id_dt, 
                    as.data.table(rowData(rse_gene))[, .(gene_id = ensemblID, gencodeID)]) %>%  data.table




setdiff(count_mat$Gene, id_dt$Gene) %>% data.table


count_mat_id <- inner_join(id_dt[, .(gencodeID, Gene)], count_mat) %>%  data.table
count_mat_id <- inner_join(data.table( gencodeID = gg), count_mat_id) %>% data.table
count_mat_id[, Gene:=gencodeID]
count_mat_id[, `:=`(gencodeID = NULL)]


count_mat_id_50 <- inner_join(id_dt[, .(gencodeID, Gene)], count_mat_50) %>%  data.table
count_mat_id_50 <- inner_join(data.table( gencodeID = gg), count_mat_id_50) %>% data.table
count_mat_id_50[, Gene:=gencodeID]
count_mat_id_50[, `:=`(gencodeID = NULL)]


count_mat_id_25 <- inner_join(id_dt[, .(gencodeID, Gene)], count_mat_25) %>%  data.table
count_mat_id_25 <- inner_join(data.table( gencodeID = gg), count_mat_id_25) %>% data.table
count_mat_id_25[, Gene:=gencodeID]
count_mat_id_25[, `:=`(gencodeID = NULL)]



gg <- count_mat_id$Gene
gg_50 <- count_mat_id_50$Gene
gg_25 <- count_mat_id_25$Gene

analysis_dt <- which(row.names(yExprs) %in% gg)
analysis_dt <- yExprs[analysis_dt, ]



# MuSiC + Darmanis -------------------------------------------------------------------




sc_data <- rna_gene_x_cells[genes %in% row.names(analysis_dt),]
sc_data <- inner_join(data.table(genes = row.names(analysis_dt)), sc_data) %>% data.table  #sc_data[match(row.names(analysis_dt), genes),]
rn <- sc_data$genes
cn <- colnames(sc_data)
cn <- inner_join(data.table(cells = cell_type_data$cells), data.table(cells = cn)) %>% data.table #cn <- cn[cn %in% cell_type_data$cells]
sc_data <- sc_data[ , cn$cells, with = F] -> darmanis_gene_cells
darmanis_gene_cells$genes <- rn
# sc_data <- sc_data[ , cn, with = F]
sc_data <- data.matrix(sc_data)
row.names(sc_data) <- rn
sc_data <- ExpressionSet(assayData = sc_data)
bulk_data <- ExpressionSet(assayData = analysis_dt)

if(type_anal == "Neurons"){
  music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells)
}else{
  music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells)
}

music_est_nnls <- data.table(samples  = row.names(music_est$Est.prop.allgene), music_est$Est.prop.allgene)
music_est <- data.table(samples  = row.names(music_est$Est.prop.allgene),music_est$Est.prop.weighted) 

music_est_nnls[, RNum := sapply(strsplit(samples, "_"), "[", 1)]
music_est[, RNum := sapply(strsplit(samples, "_"), "[", 1)]
# load("./NAc_rse_gene_withCompEsts.rda")
load(file.path("~", "cell_type_data", "NAc_rse_gene_withCompEsts_update.rda"))
music_est_nnls <- music_est_nnls[match(rse_gene$RNum, RNum),]
music_est <- music_est[match(rse_gene$RNum, RNum),]



# Estimate S_K ------------------------------------------------------------

types_cell <- unique(cell_type_data$all)
my_s_k <- sapply(types_cell, function(g){
  c_g <- cell_type_data$cells[cell_type_data$all == g]
  print(length(c_g))
  anal_data <- rna_gene_x_cells[genes %in% row.names(analysis_dt), c_g, with = F]
  anal_out <- colSums(anal_data)
  s <- mean(anal_out)
  return(s)
})
names(my_s_k) <- types_cell

sc.markers <- row.names(bulk_data)
music_basis_ests <- music_basis(sc_data, non.zero = TRUE, markers = sc.markers, 
                                clusters = cell_type_data$all, samples = cell_type_data$cells, select.ct = NULL, 
                                ct.cov = FALSE, verbose = TRUE)

music_basis_ests$M.S
my_s_k

types_cell <- unique(cell_type_data$Neurons)
my_s_k <- sapply(types_cell, function(g){
  c_g <- cell_type_data$cells[cell_type_data$Neurons == g]
  print(length(c_g))
  anal_data <- rna_gene_x_cells[genes %in% row.names(analysis_dt), c_g, with = F]
  anal_out <- colSums(anal_data)
  s <- mean(anal_out)
  return(s)
})
names(my_s_k) <- types_cell

sc.markers <- row.names(bulk_data)
music_basis_ests <- music_basis(sc_data, non.zero = TRUE, markers = sc.markers, 
                                clusters = cell_type_data$Neurons, samples = cell_type_data$cells, select.ct = NULL, 
                                ct.cov = FALSE, verbose = TRUE)

music_basis_ests$M.S
my_s_k

darmanis_sizes <- data.table(ct = names(my_s_k),
                             cs = my_s_k) %>% data.frame

# osmfish loom ------------------------------------------------------------

# Install package
my_packs <- installed.packages(fields = "Package")

if(sum(my_packs[,1] %in% "loomR") != 1){
  devtools::install_github(repo = "hhoeflin/hdf5r")
  devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
}else{
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
cell_types_aj <- cell_types_aj[Class != "Endothelial"]

osmfish_test_dt <- inner_join(my_cell_atts, cell_types_aj[, .(ClusterName, Celltype = Class, Subclass)]) %>% data.table 

osmfihs_ests <- my_cell_atts[, .(mean(CellArea), .N, var(CellArea)), by = .(ClusterName, valid)]

osmfihs_ests <- inner_join(osmfihs_ests, cell_types_aj[, .(ClusterName, Celltype = Class, Subclass)]) %>% data.table

# osmfihs_ests[, Celltype := sapply(strsplit(ClusterName, " "), "[", 1)]
# osmfihs_ests[, Celltype := ifelse(Celltype %in% c("Inhibitory", "Pyramidal", "pyramidal"), "Neuron", Celltype)] 

my_ests <- osmfihs_ests[Celltype != "Endothelial", .("w_mean" = weighted.mean(V1, N)), by = Celltype]
my_ests[, rel_size := w_mean/sum(w_mean)]

# fwrite(unique(osmfihs_ests[, .(ClusterName, Celltype)]), "celltypes.txt")

inner_join(osmfihs_ests[, .(ClusterName, "osmFISH_est" = V1, Celltype)], 
           data.table("music_est" = my_s_k, Celltype = c("Astrocyte", "Microglia", "Neuron", "Oligodendrocyte", "OPC")))

pd_data <- inner_join(my_ests[, .(Celltype, "osmFISH_est" = w_mean)], 
                      data.table("music_est" = my_s_k, Celltype = c("Astrocyte", "Microglia", "Neuron", "Oligodendrocyte", "OPC"))) %>% data.table

pd_data_unique <- pd_data[, .(osmFISH_est, Celltype, music_est)] %>% unique

# plot(pd_data$osmFISH_est, pd_data$music_est, xlab = "osmFISH", ylab = "MUSiC")


plot_dat_oshm <- my_cell_atts
plot_dat_oshm[, Celltype := NULL]
# plot_dat_oshm[, Celltype := sapply(strsplit(ClusterName, " "), "[", 1)]
plot_dat_oshm <- inner_join(plot_dat_oshm, cell_types_aj[Class != "Endothelial", .(ClusterName, Celltype = Class, Subclass)]) %>% data.table


rna_ests <- plot_dat_oshm[, .(mean(totalmolecules), .N, Celltype = unique(Celltype)), by = .(ClusterName, valid)]
rna_ests <- rna_ests[, .("w_mean" = weighted.mean(V1, N)), by = Celltype]
rna_ests[, rel_abun := w_mean/sum(w_mean)]

w_mole_count <- rna_ests$w_mean/my_ests$w_mean
w_mole_count/sum(w_mole_count)

out <- inner_join(rna_ests[, .(Celltype, "Rel_abun RNA in mouse" = rel_abun, total_mole = w_mean)], my_ests[, .(Celltype, "Rel_size in mouse" = rel_size, size = w_mean)]) %>% 
  inner_join(data.table("music_est" = my_s_k/sum(my_s_k), Celltype = c("Astrocyte", "Microglia", "Neuron", "Oligodendrocyte", "OPC"))) %>% 
  inner_join(data.table(Celltype = rna_ests$Celltype, scaled_rna = w_mole_count/sum(w_mole_count))) %>% data.table

out[, Celltype := c("Neurons", "Astrocytes", "OPC",  "Oligodendrocytes", "Microglia") ]
colMeans(music_basis_ests$S, na.rm = T)
inner_join(data.table(Celltype = names(colMeans(music_basis_ests$S, na.rm = T))),
                      out) %>% data.table

# out[match(names(colMeans(music_basis_ests$S, na.rm = T)), out$Celltype),]


{
  
  my_ests <- osmfihs_ests[Celltype != "Endothelial"]
  my_ests[, Neurons := Celltype == "Neuron"]
  
  osmfish_test_dt[, Neurons := Celltype == "Neuron"]
  if(type_anal == "Neurons"){
    my_ests[, `:=`(tot_w = sum(N)), by = Neurons]
    my_ests[, w_prime := N/tot_w]
    my_ests[, w_prime_sqrd := w_prime*w_prime]
    my_ests[, var_w_mean := sum(w_prime_sqrd*V3),  by = Neurons]
    
    my_ests <- my_ests[, .("w_mean" = weighted.mean(V1, N), var_w_mean = mean(var_w_mean), tot_n = mean(tot_w)), by = Neurons]
    
    t_stat <- diff(my_ests$w_mean)/sqrt(sum(my_ests$var_w_mean))
    tot_df <- sum(my_ests$tot_n) - 2
    
    2*pt(abs(t_stat), df = tot_df, lower.tail = F)
    
    t.test(formula = CellArea ~ Neurons, data = osmfish_test_dt, var.equal = FALSE)
    wilcox.test(formula = CellArea ~ Neurons, data = osmfish_test_dt)
   
    
    t.test(formula = totalmolecules ~ Neurons, data = osmfish_test_dt, var.equal = FALSE)
    wilcox.test(formula = totalmolecules ~ Neurons, data = osmfish_test_dt)
    
  }else{
    my_ests[, `:=`(tot_w = sum(N)), by = Celltype]
    my_ests[, w_prime := N/tot_w]
    my_ests[, w_prime_sqrd := w_prime*w_prime]
    my_ests[, var_w_mean := sum(w_prime_sqrd*V3),  by = Celltype]
    my_ests <- my_ests[, .("w_mean" = weighted.mean(V1, N), var_w_mean = mean(var_w_mean)), by = Celltype]
  }
  
  my_ests[, rel_size := w_mean/sum(w_mean)]
  
  
  plot_dat_oshm <- my_cell_atts
  plot_dat_oshm[, Celltype := NULL]
  plot_dat_oshm <- inner_join(plot_dat_oshm, cell_types_aj[Class != "Endothelial", .(ClusterName, Celltype = Class, Subclass)]) %>% data.table
  plot_dat_oshm[, Neurons := Celltype == "Neuron"]
  
  rna_ests <- plot_dat_oshm[, .(mean(totalmolecules), .N, Celltype = unique(Celltype), Neurons = unique(Neurons)), by = .(ClusterName, valid)]
  if(type_anal == "Neurons"){
    rna_ests <- rna_ests[, .("w_mean" = weighted.mean(V1, N)), by = Neurons]
  }else{
    rna_ests <- rna_ests[, .("w_mean" = weighted.mean(V1, N)), by = Celltype]
  }
  
  rna_ests[, rel_abun := w_mean/sum(w_mean)]
  
  w_mole_count <- rna_ests$w_mean/my_ests$w_mean
  
  
  
  if(type_anal == "Celltype"){
    my_out <- inner_join(rna_ests[, .(Celltype, "Rel_abun RNA in mouse" = rel_abun, total_mole = w_mean)], my_ests[, .(Celltype, "Rel_size in mouse" = rel_size, size = w_mean)]) %>% 
      inner_join(data.table(Celltype = rna_ests$Celltype, scaled_rna = w_mole_count/sum(w_mole_count))) %>% data.table
    my_out[, Celltype := c("Neurons", "Astrocytes", "OPC",  "Oligodendrocytes", "Microglia") ]
  }else{
    my_out <- inner_join(rna_ests[, .(Neurons, "Rel_abun RNA in mouse" = rel_abun, total_mole = w_mean)], my_ests[, .(Neurons, "Rel_size in mouse" = rel_size, size = w_mean)]) %>% 
      inner_join(data.table(Neurons = rna_ests$Neurons, scaled_rna = w_mole_count/sum(w_mole_count))) %>% data.table
  }
 
  
  {
    # source(file.path(".", "mycodes", "edited_music.R"))
  }
  
  if(type_anal == "Neurons"){
    # my_music_est <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
    #                               my_abundance = my_out[, .(Celltype = Neurons, w_mean = total_mole )])
    my_music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
                               cell_size = data.frame(my_out[, .(Celltype = Neurons, w_mean = total_mole )]))
  }else{
    # my_music_est <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells,
    #                               my_abundance = my_out[, .(Celltype, w_mean = total_mole )])
    my_music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells,
                               cell_size = data.frame(my_out[, .(Celltype, w_mean = total_mole )]))
  }
  my_music_est_nnls_mole <- data.table(samples  = row.names(my_music_est$Est.prop.allgene), my_music_est$Est.prop.allgene)
  my_music_est_mole <- data.table(samples  = row.names(my_music_est$Est.prop.allgene),my_music_est$Est.prop.weighted) 
  
  if(type_anal == "Neurons"){
    # my_music_est <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
    #                               my_abundance = my_out[, .(Celltype = Neurons, w_mean = size )])
    my_music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
                                  cell_size = data.frame(my_out[, .(Celltype = Neurons, w_mean = size )]))
  }else{
    # my_music_est <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells,
    #                               my_abundance = my_out[, .(Celltype, w_mean = size )])
    my_music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells,
                                  cell_size = data.frame(my_out[, .(Celltype, w_mean = size )]))
  }
  
  my_music_est_nnls_size <- data.table(samples  = row.names(my_music_est$Est.prop.allgene), my_music_est$Est.prop.allgene)
  my_music_est_size <- data.table(samples  = row.names(my_music_est$Est.prop.allgene),my_music_est$Est.prop.weighted) 
  
  
  
}

osmFISH_sizes <- my_out

# ggplot(plot_dat_oshm, aes(y = log10(totalmolecules), x = ClusterName, fill = `ClusterName` )) + 
#   geom_boxplot() +
#   theme(legend.position = "none",
#         legend.title = element_blank(),
#         axis.text.x = element_text(angle = 45, size = 8, hjust = .9, face = "bold"))
# 
# ggplot(plot_dat_oshm, aes(y = log10(totalmolecules), x = log10(CellArea), color = `ClusterName` )) + 
#   geom_point() +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         legend.direction = "horizontal")
# 
# ggsave(file.path(".", "figs", "molecule_vs_cell_type.png"))

lm(scale(totalmolecules) ~ scale(CellArea) + as.factor(Celltype), data = plot_dat_oshm[Celltype != "Excluded",]) %>% summary
glm(round(totalmolecules, 0) ~  scale(CellArea) + as.factor(Celltype), 
    data = plot_dat_oshm[Celltype != "Excluded",], family = "quasipoisson") %>% summary

# Cell type data September 9th --------------------------------------------


if(type_anal == "Neurons"){
  # my_music_est_nw <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
  #                                  my_abundance = my_out_cell_mat[, .(Celltype = Neurons, w_mean = size )])
  # my_music_est_nw_50 <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
  #                                  my_abundance = my_out_cell_mat_50[, .(Celltype = Neurons, w_mean = size )])
  # my_music_est_nw_25 <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
  #                                  my_abundance = my_out_cell_mat_25[, .(Celltype = Neurons, w_mean = size )])
  # my_music_est_v3 <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
  #                                  my_abundance = my_out_cell_mat[, .(Celltype = Neurons, w_mean = 1 )])
  
  my_music_est_nw <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
                                   cell_size = data.frame(my_out_cell_mat[, .(Celltype = Neurons, w_mean = size )]))
  my_music_est_nw_50 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
                                      cell_size = data.frame(my_out_cell_mat_50[, .(Celltype = Neurons, w_mean = size )]))
  my_music_est_nw_25 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
                                      cell_size =data.frame( my_out_cell_mat_25[, .(Celltype = Neurons, w_mean = size )]))
  my_music_est_v3 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells,
                                   cell_size = data.frame(my_out_cell_mat[, .(Celltype = Neurons, w_mean = 1 )]))
}else{
  # my_music_est_nw <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells,
  #                                  my_abundance = my_out_cell_mat[, .(Celltype = CellType, w_mean = size )])
  # my_music_est_nw_50 <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells,
  #                                  my_abundance = my_out_cell_mat_50[, .(Celltype = CellType, w_mean = size )])
  # my_music_est_nw_25 <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells,
  #                                  my_abundance = my_out_cell_mat_25[, .(Celltype = CellType, w_mean = size )])
  # my_music_est_v3 <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells,
  #                                  my_abundance = my_out_cell_mat[, .(Celltype = CellType, w_mean = 1 )])
  
  my_music_est_nw <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells,
                                cell_size = data.frame(my_out_cell_mat[, .(Celltype = CellType, w_mean = size )]))
  my_music_est_nw_50 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells,
                                   cell_size = data.frame(my_out_cell_mat_50[, .(Celltype = CellType, w_mean = size )]))
  my_music_est_nw_25 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells,
                                   cell_size = data.frame(my_out_cell_mat_25[, .(Celltype = CellType, w_mean = size )]))
  my_music_est_v3 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$all, samples = cell_type_data$cells,
                                cell_size = data.frame(my_out_cell_mat[, .(Celltype = CellType, w_mean = 1 )]))
}

my_music_est_nw <- data.table(samples  = row.names(my_music_est_nw$Est.prop.allgene), my_music_est_nw$Est.prop.weighted) 
my_music_est_nw_50 <- data.table(samples  = row.names(my_music_est_nw_50$Est.prop.allgene), my_music_est_nw_50$Est.prop.weighted) 
my_music_est_nw_25 <- data.table(samples  = row.names(my_music_est_nw_25$Est.prop.allgene), my_music_est_nw_25$Est.prop.weighted) 
my_music_est_v3 <- data.table(samples  = row.names(my_music_est_v3$Est.prop.allgene), my_music_est_v3$Est.prop.weighted) 


NAc_sizes <- rbind(my_out_cell_mat_50[, type := "50"], my_out_cell_mat_25[, type := "25"])





# Use NAc ScRNA ref


sc_data <- count_mat_id[Gene %in% row.names(analysis_dt),]
sc_data <- inner_join(data.table(Gene = row.names(analysis_dt)), sc_data) %>% data.table
rn <- sc_data$Gene
cn <- colnames(sc_data)
cn <- inner_join(data.table(cells = cell_mat$Cells), data.table(cells = cn)) %>% data.table #cn[cn %in% cell_mat$Cells]
sc_data <- sc_data[ , cn$cells, with = F] -> nac_gene_cells
nac_gene_cells$Gene <- rn
sc_data <- data.matrix(sc_data)
row.names(sc_data) <- rn
sc_data <- ExpressionSet(assayData = sc_data)
bulk_data <- ExpressionSet(assayData = analysis_dt)





if(type_anal == "Neurons"){
  music_est_v2 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells)
  # my_music_est_nw_V2 <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
  #                                     my_abundance = my_out_cell_mat[, .(Celltype = Neurons, w_mean = size )])
  # my_music_est_nac_size <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
  #                                   my_abundance = my_out[, .(Celltype = Neurons, w_mean = size )])
  # my_music_est_nac_nrna <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
  #                                   my_abundance = my_out[, .(Celltype = Neurons, w_mean = total_mole )])
  # 
  # my_music_est_nac_none <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
  #                                        my_abundance = my_out[, .(Celltype = Neurons, w_mean = 1 )])
  
  my_music_est_nw_V2 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
                                   cell_size = data.frame(my_out_cell_mat[, .(Celltype = Neurons, w_mean = size )]))
  my_music_est_nac_size <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
                                      cell_size = data.frame(my_out[, .(Celltype = Neurons, w_mean = size )]))
  my_music_est_nac_nrna <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
                                      cell_size = data.frame(my_out[, .(Celltype = Neurons, w_mean = total_mole )]))
  
  my_music_est_nac_none <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
                                      cell_size = data.frame(my_out[, .(Celltype = Neurons, w_mean = 1 )]))
  
  # my_music_est_nw_V2_50 <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
  #                                     my_abundance = my_out_cell_mat_50[, .(Celltype = Neurons, w_mean = size )])
  # my_music_est_nw_V2_25 <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
  #                                     my_abundance = my_out_cell_mat_25[, .(Celltype = Neurons, w_mean = size )])
  
  my_music_est_nw_V2_50 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
                                      cell_size = data.frame( my_out_cell_mat_50[, .(Celltype = Neurons, w_mean = size )]))
  my_music_est_nw_V2_25 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
                                      cell_size = data.frame(my_out_cell_mat_25[, .(Celltype = Neurons, w_mean = size )]))
  
  
  my_music_est_nw_darmanis <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
                                   cell_size = darmanis_sizes)
}else{
  music_est_v2 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$CellType, samples = cell_mat$Cells)
  # my_music_est_nw_V2 <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$CellType, samples = cell_mat$Cells,
  #                                     my_abundance = my_out_cell_mat[, .(Celltype = CellType, w_mean = size )])
  # 
  # my_music_est_nac_size <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$CellType, samples = cell_mat$Cells,
  #                                        my_abundance = my_out[, .(Celltype, w_mean = size )])
  # my_music_est_nac_nrna <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$CellType, samples = cell_mat$Cells,
  #                                       my_abundance = my_out[, .(Celltype, w_mean = size )])
  # 
  # my_music_est_nac_none <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
  #                                        my_abundance = my_out[, .(Celltype = Neurons, w_mean = 1 )])
  # 
  # my_music_est_nw_V2_50 <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
  #                                        my_abundance = my_out_cell_mat_50[, .(Celltype, w_mean = size )])
  # my_music_est_nw_V2_25 <- my_music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
  #                                        my_abundance = my_out_cell_mat_25[, .(Celltype, w_mean = size )])
  
  my_music_est_nw_V2 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$CellType, samples = cell_mat$Cells,
                                   cell_size = data.frame(my_out_cell_mat[, .(Celltype = CellType, w_mean = size )]))
  
  my_music_est_nac_size <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$CellType, samples = cell_mat$Cells,
                                      cell_size = data.frame(my_out[, .(Celltype, w_mean = size )]))
  my_music_est_nac_nrna <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$CellType, samples = cell_mat$Cells,
                                      cell_size = data.frame(my_out[, .(Celltype, w_mean = size )]))
  
  my_music_est_nac_none <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$CellType, samples = cell_mat$Cells,
                                      cell_size = data.frame(my_out[, .(Celltype = Neurons, w_mean = 1 )]))
  
  my_music_est_nw_V2_50 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$CellType, samples = cell_mat$Cells,
                                      cell_size = data.frame(my_out_cell_mat_50[, .(Celltype, w_mean = size )]))
  my_music_est_nw_V2_25 <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$CellType, samples = cell_mat$Cells,
                                      cell_size = data.frame(my_out_cell_mat_25[, .(Celltype, w_mean = size )]))
  
  my_music_est_nw_darmanis <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$CellType, samples = cell_mat$Cells,
                                         cell_size = darmanis_sizes)
}


music_est_nnls_v2 <- data.table(samples  = row.names(music_est_v2$Est.prop.allgene), music_est_v2$Est.prop.allgene)
music_est_v2 <- data.table(samples  = row.names(music_est_v2$Est.prop.allgene), music_est_v2$Est.prop.weighted) 
my_music_est_nw_V2 <- data.table(samples  = row.names(my_music_est_nw_V2$Est.prop.allgene), my_music_est_nw_V2$Est.prop.weighted) 

my_music_est_size_nac <- data.table(samples  = row.names(my_music_est_nac_size$Est.prop.allgene), my_music_est_nac_size$Est.prop.weighted) 
my_music_est_nrna_nac <- data.table(samples  = row.names(my_music_est_nac_nrna$Est.prop.allgene), my_music_est_nac_nrna$Est.prop.weighted) 
my_music_est_nac_none <- data.table(samples  = row.names(my_music_est_nac_none$Est.prop.allgene), my_music_est_nac_none$Est.prop.weighted) 

my_music_est_nw_V2_50 <- data.table(samples  = row.names(my_music_est_nw_V2_50$Est.prop.allgene), my_music_est_nw_V2_50$Est.prop.weighted) 
my_music_est_nw_V2_25 <- data.table(samples  = row.names(my_music_est_nw_V2_25$Est.prop.allgene), my_music_est_nw_V2_25$Est.prop.weighted) 

my_music_est_nw_darmanis <- data.table(samples  = row.names(my_music_est_nw_darmanis$Est.prop.allgene), my_music_est_nw_darmanis$Est.prop.weighted) 

if(type_anal == "Celltype"){
  fwrite(music_est_v2, file = file.path("cell_type_data", "music_est_nac_default.txt"))
}


if(type_anal == "Neurons"){
  round(cor(music_est$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_nw$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_nw_25$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_nw_50$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_nw_V2$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(music_est_v2$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_v3$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_mole$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_size$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  
  round(cor(my_music_est_nrna_nac$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_size_nac$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  
  round(cor(my_music_est_nw_V2_25$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_nw_V2_50$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  
  
  round(cor(my_music_est_nw_darmanis$`TRUE`, rse_gene$NeuN_pos_DNAm)^2, 4)
  
  # RMSE
  round(sqrt(mean(music_est$`TRUE`-rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_nw$`TRUE`- rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_nw_25$`TRUE` - rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_nw_50$`TRUE` - rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_nw_V2$`TRUE` - rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(music_est_v2$`TRUE` - rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_v3$`TRUE` - rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_mole$`TRUE` - rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_size$`TRUE` - rse_gene$NeuN_pos_DNAm)^2), 4)
  
  round(sqrt(mean(my_music_est_nrna_nac$`TRUE` - rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_size_nac$`TRUE` - rse_gene$NeuN_pos_DNAm)^2), 4)
  
  round(sqrt(mean((my_music_est_nw_V2_25$`TRUE` - rse_gene$NeuN_pos_DNAm)^2)), 4)
  round(sqrt(mean((my_music_est_nw_V2_50$`TRUE` - rse_gene$NeuN_pos_DNAm)^2)), 4)
  
  round(sqrt(mean((my_music_est_nw_darmanis$`TRUE` - rse_gene$NeuN_pos_DNAm)^2)), 4)
  
  
}else{
  round(cor(music_est$`Neurons`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_nw$`Neurons`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_nw_V2$`Neurons`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(music_est_v2$`Neurons`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_v3$`Neurons`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_mole$`Neurons`, rse_gene$NeuN_pos_DNAm)^2, 4)
  round(cor(my_music_est_size$`Neurons`, rse_gene$NeuN_pos_DNAm)^2, 4)
  
  # RMSE
  round(sqrt(mean(music_est$`Neurons`-rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_nw$`Neurons`- rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_nw_25$`Neurons` - rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_nw_50$`Neurons` - rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_nw_V2$`Neurons` - rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(music_est_v2$`Neurons` - rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_v3$`Neurons` - rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_mole$`Neurons` - rse_gene$NeuN_pos_DNAm)^2), 4)
  round(sqrt(mean(my_music_est_size$`Neurons` - rse_gene$NeuN_pos_DNAm)^2), 4)
}



# Plot
{
  if(type_anal == "Neurons"){
    
    out_pd <- rbind(inner_join(music_est[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_mole[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[osmFISH totalRNA]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_size[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[osmFISH cellsize]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(music_est_v2[, .(samples, est = `TRUE`, type = "Ref:[NAc];Size:[Default]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_nw[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default:Nac all genes]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_nw_50[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default:Nac top 50 genes]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_nw_25[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default:Nac top 25 genes]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_nw_V2[, .(samples, est = `TRUE`, type = "Ref:[NAc];Size:[Default:Nac all genes]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_v3[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[None]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm))) %>% data.table
    
    out_pd_v2 <- out_pd[type != "Ref:[Darmanis];Size:[None]"]
    out_pd <- out_pd[!type %in% c("Ref:[Darmanis];Size:[Nac top 50 genes]",
                                  "Ref:[Darmanis];Size:[Nac top 25 genes]")]
    
    out_pd_v3 <- rbind(inner_join(music_est[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_mole[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[osmFISH totalRNA]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_size[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[osmFISH cellsize]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_nw[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default:Nac all genes]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_nw_50[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default:Nac top 50 genes]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_nw_25[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default:Nac top 25 genes]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_v3[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[None]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm))) %>% data.table
    
    
    out_pd_v4 <- rbind(inner_join(music_est_v2[, .(samples, est = `TRUE`, type = "Ref:[NAc];Size:[Default]")],
                                 data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                      inner_join(my_music_est_nw_V2[, .(samples, est = `TRUE`, type = "Ref:[NAc];Size:[Default:Nac all genes]")],
                                 data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                      inner_join(my_music_est_nrna_nac[, .(samples, est = `TRUE`, type = "Ref:[NAc];Size:[osmFISH totalRNA]")],
                                 data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                      inner_join(my_music_est_size_nac[, .(samples, est = `TRUE`, type = "Ref:[NAc];Size:[osmFISH cellsize]")],
                                 data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                      inner_join(my_music_est_nw_V2_50[, .(samples, est = `TRUE`, type = "Ref:[NAc];Size:[Default:Nac top 50 genes]")],
                                 data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                      inner_join(my_music_est_nw_V2_25[, .(samples, est = `TRUE`, type = "Ref:[NAc];Size:[Default:Nac top 25 genes]")],
                                 data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                      inner_join(my_music_est_nac_none[, .(samples, est = `TRUE`, type = "Ref:[NAc];Size:[None]")],
                                 data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm))) %>% data.table
    
    
    
    
  }else{
    out_pd <- rbind(inner_join(music_est[, .(samples, est = Neurons, type = "Ref:[Darmanis];Size:[Default]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_mole[, .(samples, est = Neurons, type = "Ref:[Darmanis];Size:[osmFISH totalRNA]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_size[, .(samples, est = Neurons, type = "Ref:[Darmanis];Size:[osmFISH cellsize]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(music_est_v2[, .(samples, est = Neurons, type = "Ref:[NAc];Size:[Default]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_nw[, .(samples, est = Neurons, type = "Ref:[Darmanis];Size:[Default:Nac all genes]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_nw_50[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default:Nac top 50 genes]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_nw_25[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default:Nac top 25 genes]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_nw_V2[, .(samples, est = Neurons, type = "Ref:[NAc];Size:[Default:Nac all genes]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                    inner_join(my_music_est_v3[, .(samples, est = Neurons, type = "Ref:[Darmanis];Size:[None]")],
                               data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm))) %>% data.table 
    out_pd_v2 <- out_pd[type != "Ref:[Darmanis];Size:[None]"]
    out_pd <- out_pd[!type %in% c("Ref:[Darmanis];Size:[Nac top 50 genes]",
                                  "Ref:[Darmanis];Size:[Nac top 25 genes]")]
    
    out_pd_v3 <- rbind(inner_join(music_est[, .(samples, est = Neurons, type = "Ref:[Darmanis];Size:[Default]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_mole[, .(samples, est = Neurons, type = "Ref:[Darmanis];Size:[osmFISH totalRNA]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_size[, .(samples, est = Neurons, type = "Ref:[Darmanis];Size:[osmFISH cellsize]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_nw[, .(samples, est = Neurons, type = "Ref:[Darmanis];Size:[Default:Nac all genes]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_nw_50[, .(samples, est = Neurons, type = "Ref:[Darmanis];Size:[Default:Nac top 50 genes]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_nw_25[, .(samples, est = Neurons, type = "Ref:[Darmanis];Size:[Default:Nac top 25 genes]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_v3[, .(samples, est = Neurons, type = "Ref:[Darmanis];Size:[None]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm))) %>% data.table
    
    
    out_pd_v4 <- rbind(inner_join(music_est_v2[, .(samples, est = Neurons, type = "Ref:[NAc];Size:[Default]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_nw_V2[, .(samples, est = Neurons, type = "Ref:[NAc];Size:[Default:Nac all genes]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_nrna_nac[, .(samples, est = Neurons, type = "Ref:[NAc];Size:[osmFISH totalRNA]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_size_nac[, .(samples, est = Neurons, type = "Ref:[NAc];Size:[osmFISH cellsize]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_nw_V2_50[, .(samples, est = Neurons, type = "Ref:[NAc];Size:[Default:Nac top 50 genes]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_nw_V2_25[, .(samples, est = Neurons, type = "Ref:[NAc];Size:[Default:Nac top 25 genes]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                       inner_join(my_music_est_nac_none[, .(samples, est = Neurons, type = "Ref:[NAc];Size:[None]")],
                                  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm))) %>% data.table
    
  }
  
  
  mytable <- out_pd[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                        RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = type]
  
  mytable[, `:=`(R_sqrd = round(R_sqrd, 4),
                 RMSE = round(RMSE, 4))]
  names(mytable)[1] <- "Approach"
  mytable <- mytable[order(RMSE)]
  mytable <- mytable[Approach != "Ref:[NAc];Size:[Default:Nac all genes]"]
  
  mytable_v2 <- out_pd_v2[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                        RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = type]
  
  mytable_v2[, `:=`(R_sqrd = round(R_sqrd, 4),
                 RMSE = round(RMSE, 4))]
  names(mytable_v2)[1] <- "Approach"
  mytable_v2 <- mytable_v2[order(RMSE)]
  
  
  out_pd <- out_pd[type != "Ref:[NAc];Size:[Default:Nac all genes]"]
  
  
  p_save_cols <- RColorBrewer::brewer.pal(length(unique(out_pd$type)), "Dark2")
  names(p_save_cols) <- unique(out_pd$type)
  
  p_save_shape <- c(15:19, 22:24)
  names(p_save_shape) <- unique(out_pd$type)
  
  transparent_legend =  theme(
    legend.background = element_rect(fill ="transparent"),
    legend.key = element_rect(fill = "transparent",
                              color = "transparent")
  )
  
  remove_grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  no_x_axis_label <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
  out_pd$type <- factor(out_pd$type, levels = mytable$Approach)
  
  set <- c("Ref:[Darmanis];Size:[Default]",
           "Ref:[Darmanis];Size:[osmFISH totalRNA]",
           "Ref:[Darmanis];Size:[osmFISH cellsize]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd[type %in%set ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 3) + 
    geom_abline(slope = 1, intercept = 0) +
    xlab("MuSiC estimates") +
    scale_color_manual(values = p_save_cols, name = "" ) +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx,], rows = NULL), xmin=0.4, xmax=.4, ymin=.35, ymax=.3) +
    transparent_legend + remove_grid +
    ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          text = element_text(size = 12),
          axis.title = element_text(face="bold", size = 9),
          axis.text.y=element_text(size = 8, face="bold"),
          axis.text.x=element_text(size = 8, face="bold"),
          legend.position = c(0.85,0.10),
          legend.title = element_text(face="bold"))
  ggsave(file.path(".", "model", "bias_correction.png"), plot = p_save, dpi = "retina", width = 40, height = 30, units = "cm")
  
  
  
  
  p_save <- ggplot(data = out_pd, aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 6) + 
    geom_abline(slope = 1, intercept = 0) +
    xlab("MuSiC estimates") +
    scale_color_manual(values = p_save_cols, name = "") +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable, rows = NULL), xmin=0.24, xmax=.1, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          text = element_text(size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text.y=element_text(size = 12, face="bold"),
          axis.text.x=element_text(size = 12, face="bold"),
          legend.position = c(0.85,0.10),
          legend.title = element_text(face="bold"), 
          legend.text = element_text(size = 12, face="bold"))
  ggsave(file.path(".", "model", "bias_correction_v2.png"), plot = p_save, dpi = "retina", width = 60, height = 50, units = "cm")
  
  
  set <- c("Ref:[NAc];Size:[Default]",
           "Ref:[NAc];Size:[Default:Nac all genes]")
  set_idx <- which(!mytable$Approach %in% set)
  
  p_save <- ggplot(data = out_pd[!type %in% set ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 3) + 
    geom_abline(slope = 1, intercept = 0) +
    xlab("MuSiC estimates") +
    scale_color_manual(values = p_save_cols, name = "" ) +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx], rows = NULL), xmin=0.22, xmax=.22, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          text = element_text(size = 12),
          axis.title = element_text(face="bold", size = 9),
          axis.text.y=element_text(size = 8, face="bold"),
          axis.text.x=element_text(size = 8, face="bold"),
          legend.position = c(0.85, 0.10),
          legend.title = element_text(face="bold"))
  ggsave(file.path(".", "model", "bias_correction_v3_Darmanis_ref.png"), plot = p_save, dpi = "retina", width = 40, height = 30, units = "cm")
  
  set <- c("Ref:[NAc];Size:[Default]",
           "Ref:[NAc];Size:[Default:Nac all genes]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd[type %in% set], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 3) + 
    geom_abline(slope = 1, intercept = 0) +
    xlab("MuSiC estimates") +
    scale_color_manual(values = p_save_cols, name = "" ) +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx], rows = NULL), xmin=0.04, xmax=.04, ymin=.4, ymax=.4) +
    transparent_legend + remove_grid +
    ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          text = element_text(size = 12),
          axis.title = element_text(face="bold", size = 9),
          axis.text.y=element_text(size = 8, face="bold"),
          axis.text.x=element_text(size = 8, face="bold"),
          legend.position = c(0.85, 0.10),
          legend.title = element_text(face="bold"))
  ggsave(file.path(".", "model", "bias_correction_v3_nac_ref.png"), plot = p_save, dpi = "retina", width = 40, height = 30, units = "cm")
  
  
  t <- out_pd$type %>% unique
  sapply(1:n_distinct(out_pd$type), function(i){
    p_save <- ggplot(data = out_pd[type == t[i], ], aes(x = est, y = `Houseman DNAm-based`, color = type)) + 
      geom_point(size = 3) + 
      geom_abline(slope = 1, intercept = 0) +
      xlab("MuSiC estimates") +
      scale_color_manual(values = p_save_cols, name = "" ) + 
      transparent_legend + remove_grid +
      ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
      theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            text = element_text(size = 12),
            axis.title = element_text(face="bold", size = 9),
            axis.text.y=element_text(size = 8, face="bold"),
            axis.text.x=element_text(size = 8, face="bold"),
            legend.position = c(0.2,0.85),
            legend.title = element_text(face="bold"))
    ggsave(file.path(".", "model", paste0("bias_correction_type", i, ".png")), plot = p_save, dpi = "retina", width = 30, height = 20, units = "cm")
  })
  
  # Main figs


  if(!dir.exists(file.path(".", "model", "manuscript"))){
    dir.create(file.path(".", "model", "manuscript"))
  }
  # •	A: Houseman RNA with Darmanis
  plt_dt <-  data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)
  # Get Houseman RNA estimates
  e <- new.env()
  load("./NAc_rse_gene_withCompEsts.rda", envir = e)
  plt_dt <- inner_join(plt_dt, data.table(samples = e$rse_gene$SampleID, "Houseman RNA-based" = e$rse_gene$Neurons_RNA)) %>% data.table
  
  tble <- plt_dt[, .(R_sqrd = cor(`Houseman RNA-based`, `Houseman DNAm-based`)^2,
                        RMSE = sqrt(mean((`Houseman RNA-based`-`Houseman DNAm-based`)^2)))]
  
  tble[, `:=`(R_sqrd = round(R_sqrd, 4),
                 RMSE = round(RMSE, 4))]
  
  mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 2.0)),
    colhead = list(fg_params=list(cex = 2.0)),
    rowhead = list(fg_params=list(cex =2.0)))
  
  
  
  p_save <- ggplot(data = plt_dt, aes(x = `Houseman RNA-based`,y = `Houseman DNAm-based`)) + 
    geom_point(size = 10, col = "skyblue") + 
    geom_abline(slope = 1, intercept = 0) +
    labs(x = "Houseman RNA-based", size = rel(1.5)) + 
    transparent_legend + remove_grid +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0,.8)) +
    annotation_custom(gridExtra::tableGrob(tble, rows = NULL, theme = mytheme), xmin=0.65, xmax=.65, ymin=.05, ymax=.05) +
    ggtitle("Houseman Methylation vs Houseman RNA based\n cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = c(0.2,0.85),
          legend.title = element_text(face="bold"),
          plot.margin = margin(5,10,10,2))
  ggsave(file.path(".", "model", "manuscript", "houseman_RNA_DNA.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")
  
  # •	B: Music RNA Default with Darmanis
  set <- c("Ref:[Darmanis];Size:[Default]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd[type %in%set, ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0,.8)) +
    labs(x = "MuSiC estimates", size = rel(1.5)) +
    scale_color_manual(values = p_save_cols, name = "" ) +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx, .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.65, xmax=.65, ymin=.05, ymax=.05) +
    transparent_legend + remove_grid +
    ggtitle("Houseman Methylation vs MuSiC \ncell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = "none",
          legend.title = element_text(face="bold"),
          plot.margin = margin(5,10,10,2))
  ggsave(file.path(".", "model", "manuscript", "dna_v_music.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")
  
  
  # •	C: CIBERSORT with Darmanis
  
  DNam_darmanis <- data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)
  ciber_darmanis <- fread(file.path("cell_type_data", "CIBERSORT.Output_darmanis_no_quantile_V2.txt"))
  
  ciber_darmanis[, type:="Ref:[Darmanis]"]
  out_ciber <- inner_join(ciber_darmanis[, .(samples = `Input Sample`, Neurons_pos, type)], DNam_darmanis) %>% data.table
  
  mytable_ciber <- out_ciber[, .(R_sqrd = cor(Neurons_pos, `Houseman DNAm-based`)^2,
                                 RMSE = sqrt(mean((Neurons_pos -`Houseman DNAm-based`)^2)))]
  
  mytable_ciber[, `:=`(R_sqrd = round(R_sqrd, 4),
                       RMSE = round(RMSE, 4))]
  mytable_ciber <- mytable_ciber[order(RMSE)]
  
  
  # Plot
  p_save_cols <- RColorBrewer::brewer.pal(3, "Dark2")[2]
  names(p_save_cols) <- unique(out_ciber$type)
  
  p_save_shape <- c(16)
  names(p_save_shape) <- unique(out_ciber$type)
  
  p_save_ciber <- ggplot(data = out_ciber, aes(x = Neurons_pos, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    labs(x = "CIBERSORT", size = rel(1.5)) +
    scale_color_manual(values = p_save_cols, name = "" ) +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0,.5)) +
    annotation_custom(gridExtra::tableGrob(mytable_ciber, rows = NULL, theme = mytheme), xmin=0.4, xmax=.4, ymin=.05, ymax=.05) +
    transparent_legend + remove_grid +
    ggtitle("Houseman Methylation vs CIBERSORT\n cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = "none",
          legend.title = element_text(face="bold"),
          plot.margin = margin(5,10,10,2))
  ggsave(file.path(".", "model",  "manuscript", "dna_v_cibersort.pdf"), plot = p_save_ciber, dpi = "retina", width = 20, height = 20, units = "cm")
  
  
  # Darmanis (osmfish)
  p_save_cols <- RColorBrewer::brewer.pal(length(unique(out_pd$type)), "Dark2")
  names(p_save_cols) <- unique(out_pd$type)
  p_save_shape <- c(15:19, 22:24)
  names(p_save_shape) <- unique(out_pd$type)
  
  
  set <- c("Ref:[Darmanis];Size:[osmFISH cellsize]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd[type %in%set, ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .8)) +
    labs(x = "MuSiC + osmFISH cell area", size = rel(1.5)) +
    scale_color_manual(values = p_save_cols, name = "" ) +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx, .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.15, xmax=.15, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    # ggtitle("Houseman Methylation vs MuSiC \ncell-type proportion estimates\n using osmFISH cellsize\n (Neurons only)") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = "none",
          legend.title = element_text(face="bold"))
  ggsave(file.path(".", "model", "manuscript", "osmFISH_cellsize_v_dna.png"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")
  
  
  
  set <- c("Ref:[Darmanis];Size:[osmFISH totalRNA]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd[type %in%set, ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .8)) +
    labs(x = "MuSiC + osmFISH totalRNA", size = rel(1.5)) +
    scale_color_manual(values = p_save_cols, name = "" ) +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx, .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.15, xmax=.15, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    # ggtitle("Houseman Methylation\ vs MuSiC\n cell-type proportion estimates\n using total RNA count as cellsize\n (Neurons only)") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = "none",
          legend.title = element_text(face="bold"))
  ggsave(file.path(".", "model", "manuscript", "osmFISH_mole_v_dna.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")
  
  # set <- c("Ref:[Darmanis];Size:[osmFISH cellsize]",
  #          "Ref:[Darmanis];Size:[osmFISH totalRNA]")
  # set_idx <- which(mytable$Approach %in% set)
  # p_save <- ggplot(data = out_pd[type %in%set, ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
  #   geom_point(size = 10) + 
  #   geom_abline(slope = 1, intercept = 0) +
  #   scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .8)) +
  #   labs(x = "MuSiC estimates", size = rel(1.5)) +
  #   scale_color_manual(values = p_save_cols, name = "" ) +
  #   scale_shape_manual(values = p_save_shape, name = "" ) +
  #   annotation_custom(gridExtra::tableGrob(mytable[set_idx, .(`Cellsize estimate` = Approach, R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.15, xmax=.15, ymin=.38, ymax=.38) +
  #   transparent_legend + remove_grid +
  #   ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
  #   theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
  #         text = element_text(size = 20),
  #         axis.title = element_text(face="bold", size = 30),
  #         axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
  #         axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
  #         legend.position = "none",
  #         legend.title = element_text(face="bold"))
  # ggsave(file.path(".", "model", "manuscript", "osmFISH_v_dna.pdf"), plot = p_save, dpi = "retina", width = 40, height = 40, units = "cm")
  
  
  set <- c("Ref:[Darmanis];Size:[Default:Nac top 50 genes]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd[type %in%set, ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
    labs(x = "MuSiC + top 50 genes", size = rel(1.5)) +
    scale_color_manual(values = p_save_cols, name = "" ) +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx, .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.1, xmax=.1, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = "none",
          legend.title = element_text(face="bold"))
  ggsave(file.path(".", "model", "manuscript", "top50_v_dna.png"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")
  
  set <- c("Ref:[Darmanis];Size:[Default:Nac top 25 genes]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd[type %in%set, ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
    labs(x = "MuSiC + top 25 genes", size = rel(1.5)) +
    scale_color_manual(values = p_save_cols, name = "" ) +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx, .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.1, xmax=.1, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = "none",
          legend.title = element_text(face="bold"))
  ggsave(file.path(".", "model", "manuscript", "top25_v_dna.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")
  
  set <- c("Ref:[Darmanis];Size:[Default:Nac all genes]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd[type %in%set, ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
    xlab("MuSiC + all genes") +
    scale_color_manual(values = p_save_cols, name = "" ) +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx, .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.1, xmax=.1, ymin=.4, ymax=.4) +
    transparent_legend + remove_grid +
    # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = "none",
          legend.title = element_text(face="bold"))
  ggsave(file.path(".", "model", "manuscript", "dar_nacgenes_v_dna.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")
  
  set <- c("Ref:[Darmanis];Size:[None]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd[type %in%set, ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .8)) +
    xlab("MuSiC + no cell size genes") +
    scale_color_manual(values = p_save_cols, name = "" ) +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx, .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.1, xmax=.1, ymin=.4, ymax=.4) +
    transparent_legend + remove_grid +
    # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = "none",
          legend.title = element_text(face="bold"))
  ggsave(file.path(".", "model", "manuscript", "dar_no_cell_size_v_dna.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")
  
  
  set <- c("Ref:[Darmanis];Size:[Default:Nac all genes]",
           "Ref:[Darmanis];Size:[Default:Nac top 25 genes]",
           "Ref:[Darmanis];Size:[Default:Nac top 50 genes]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd[type %in%set, ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
    xlab("MuSiC estimates") +
    scale_color_manual(values = p_save_cols, name = "" ) +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx, .(`Cellsize estimate` = Approach, R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.1, xmax=.1, ymin=.4, ymax=.4) +
    transparent_legend + remove_grid +
    ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
          text = element_text(size = 32),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = 25, face="bold", hjust = 0.5),
          axis.text.x=element_text(size = 25, face="bold", hjust = 0.5),
          legend.position = c(.8, .1),
          legend.title = element_text(face="bold"))
  ggsave(file.path(".", "model", "manuscript", "all_nacgenes_sizes_v_music.pdf"), plot = p_save, dpi = "retina", width = 40, height = 40, units = "cm")
  
  
  # NAc reference
  
  
  mytable <- out_pd_v4[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                           RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = type]
  
  mytable[, `:=`(R_sqrd = round(R_sqrd, 4),
                 RMSE = round(RMSE, 4))]
  names(mytable)[1] <- "Approach"
  mytable <- mytable[order(RMSE)]
  
  out_pd_v4$type <- factor(out_pd_v4$type, levels = mytable$Approach)
  
  p_save_cols <- RColorBrewer::brewer.pal(length(unique(out_pd_v4$type)), "Dark2")
  names(p_save_cols) <- unique(out_pd_v4$type)
  
  p_save_shape <- c(15:19, 22:23)
  names(p_save_shape) <- unique(out_pd_v4$type)
  
  set <- c("Ref:[NAc];Size:[osmFISH cellsize]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd_v4[type %in%set, ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
    labs(x = "MuSiC + osmFISH cell area", size = rel(1.5)) +
    scale_color_manual(values = p_save_cols, name = "") +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx, .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.1, xmax=.1, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = "none",
          legend.title = element_text(face="bold"), 
          legend.text = element_text(size = 12, face="bold"))
  ggsave(file.path(".", "model", "manuscript", "osmFISH_cellsize_v_dna_nac_ref.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")
  
  
  set <- c("Ref:[NAc];Size:[osmFISH totalRNA]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd_v4[type %in%set, ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
    labs(x = "MuSiC + osmFISH totalRNA", size = rel(1.5)) +
    scale_color_manual(values = p_save_cols, name = "") +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx, .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.1, xmax=.1, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = "none",
          legend.title = element_text(face="bold"), 
          legend.text = element_text(size = 12, face="bold"))
  ggsave(file.path(".", "model", "manuscript", "osmFISH_mole_v_dna_nac_ref.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")
  
  
  set <- c("Ref:[NAc];Size:[Default]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd_v4[type %in%set, ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
    labs(x = "MuSiC default", size = rel(1.5)) +
    scale_color_manual(values = p_save_cols, name = "") +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx, .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.24, xmax=.1, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = "none",
          legend.title = element_text(face="bold"), 
          legend.text = element_text(size = 12, face="bold"))
  ggsave(file.path(".", "model", "manuscript", "nac_ref_v_dna.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")
  
  set <- c("Ref:[NAc];Size:[None]")
  set_idx <- which(mytable$Approach %in% set)
  p_save <- ggplot(data = out_pd_v4[type %in%set, ], aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
    labs(x = "MuSiC no cell size", size = rel(1.5)) +
    scale_color_manual(values = p_save_cols, name = "") +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable[set_idx, .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.24, xmax=.1, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = "none",
          legend.title = element_text(face="bold"), 
          legend.text = element_text(size = 12, face="bold"))
  ggsave(file.path(".", "model", "manuscript", "nac_ref_no_cell_size_v_dna.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")
  
  
  pt_dt <- inner_join(my_music_est_nw_darmanis[, .(samples, est = `TRUE`, type = "Ref:[NAc];Size:[ Darmanis cell size]")],
                      data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)) %>% data.table
  
  mytable <- pt_dt[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                           RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2)))]
  
  mytable[, `:=`(R_sqrd = round(R_sqrd, 4),
                 RMSE = round(RMSE, 4))]

  
  
  p_save_cols <- RColorBrewer::brewer.pal(3, "Dark2")[1]
  names(p_save_cols) <- unique(pt_dt$type)
  
  
  p_save <- ggplot(data = pt_dt, aes(x = est, y = `Houseman DNAm-based`, color = type)) + 
    geom_point(size = 10) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
    labs(x = "MuSiC + Darmanis Cell size", size = rel(1.5)) +
    scale_color_manual(values = p_save_cols, name = "") +
    annotation_custom(gridExtra::tableGrob(mytable[, .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.08, xmax=0.08, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          text = element_text(size = 20),
          axis.title = element_text(face="bold", size = 30),
          axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
          legend.position = "none",
          legend.title = element_text(face="bold"), 
          legend.text = element_text(size = 12, face="bold"))
  ggsave(file.path(".", "model", "manuscript", "nac_ref_darmanis_v_dna.png"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")
  
  
  
  
  
  # More comparisons/remove cell size == 1
  p_save_cols <- RColorBrewer::brewer.pal(length(unique(out_pd_v2$type)), "Dark2")
  names(p_save_cols) <- unique(out_pd_v2$type)
  
  p_save_shape <- c(14:19, 22:23)
  names(p_save_shape) <- unique(out_pd_v2$type)
  
  out_pd_v2$type <- factor(out_pd_v2$type, levels = mytable_v2$Approach)
  p_save <- ggplot(data = out_pd_v2, aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 6) + 
    geom_abline(slope = 1, intercept = 0) +
    xlab("MuSiC estimates") +
    scale_color_manual(values = p_save_cols, name = "") +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable_v2, rows = NULL), xmin=0.24, xmax=.1, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          text = element_text(size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text.y=element_text(size = 12, face="bold"),
          axis.text.x=element_text(size = 12, face="bold"),
          legend.position = c(0.85,0.10),
          legend.title = element_text(face="bold"), 
          legend.text = element_text(size = 12, face="bold"))
  ggsave(file.path(".", "model", "bias_correction_new.png"), plot = p_save, dpi = "retina", width = 60, height = 50, units = "cm")
  
  
  # Darmanis ref only
  mytable <- out_pd_v3[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                        RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = type]
  
  mytable[, `:=`(R_sqrd = round(R_sqrd, 4),
                 RMSE = round(RMSE, 4))]
  names(mytable)[1] <- "Approach"
  mytable <- mytable[order(RMSE)]
  
  out_pd_v3$type <- factor(out_pd_v3$type, levels = mytable$Approach)
  
  p_save_cols <- RColorBrewer::brewer.pal(length(unique(out_pd_v3$type)), "Dark2")
  names(p_save_cols) <- unique(out_pd_v3$type)
  
  p_save_shape <- c(15:19, 22:23)
  names(p_save_shape) <- unique(out_pd_v3$type)
  
  
  p_save <- ggplot(data = out_pd_v3, aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 6) + 
    geom_abline(slope = 1, intercept = 0) +
    xlab("MuSiC estimates") +
    scale_color_manual(values = p_save_cols, name = "") +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable, rows = NULL), xmin=0.24, xmax=.1, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          text = element_text(size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text.y=element_text(size = 12, face="bold"),
          axis.text.x=element_text(size = 12, face="bold"),
          legend.position = c(0.85,0.10),
          legend.title = element_text(face="bold"), 
          legend.text = element_text(size = 12, face="bold"))
  ggsave(file.path(".", "model", "bias_correction_new_darmanis_ref.png"), plot = p_save, dpi = "retina", width = 60, height = 50, units = "cm")
  
  
  # NAc reference only
  mytable <- out_pd_v4[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                           RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = type]
  
  mytable[, `:=`(R_sqrd = round(R_sqrd, 4),
                 RMSE = round(RMSE, 4))]
  names(mytable)[1] <- "Approach"
  mytable <- mytable[order(RMSE)]
  
  out_pd_v4$type <- factor(out_pd_v4$type, levels = mytable$Approach)
  
  p_save_cols <- RColorBrewer::brewer.pal(length(unique(out_pd_v4$type)), "Dark2")
  names(p_save_cols) <- unique(out_pd_v4$type)
  
  p_save_shape <- c(15:19, 22:23)
  names(p_save_shape) <- unique(out_pd_v4$type)
  
  
  p_save <- ggplot(data = out_pd_v4, aes(x = est, y = `Houseman DNAm-based`, color = type, shape = type)) + 
    geom_point(size = 6) + 
    geom_abline(slope = 1, intercept = 0) +
    xlab("MuSiC estimates") +
    scale_color_manual(values = p_save_cols, name = "") +
    scale_shape_manual(values = p_save_shape, name = "" ) +
    annotation_custom(gridExtra::tableGrob(mytable, rows = NULL), xmin=0.24, xmax=.1, ymin=.38, ymax=.38) +
    transparent_legend + remove_grid +
    ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          text = element_text(size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text.y=element_text(size = 12, face="bold"),
          axis.text.x=element_text(size = 12, face="bold"),
          legend.position = c(0.85,0.10),
          legend.title = element_text(face="bold"), 
          legend.text = element_text(size = 12, face="bold"))
  ggsave(file.path(".", "model", "bias_correction_new_NAc_ref.png"), plot = p_save, dpi = "retina", width = 60, height = 50, units = "cm")
  
  
}


# NNLS --------------------------------------------------------------------

# Darmanis

pb <- txtProgressBar(min = 0, max = nrow(darmanis_gene_cells) , width = NA, style = 3)
pb_i <- 0

n_t_idx <- which(cell_type_data$Neurons == T)
n_f_idx <- which(cell_type_data$Neurons == F)
darmanis_gene_ave <- apply(data.matrix(darmanis_gene_cells[, -266]), 1 ,function(x){
  
  out <- data.table(gene = NA, exp_neu_n = mean(x[n_f_idx]),
                    exp_neu_p  = mean(x[n_t_idx]))
  
  pb_i <<- pb_i+1
  setTxtProgressBar(pb, pb_i, title = paste(round(pb_i/nrow(darmanis_gene_cells))*100,"% done"))
  
  return(out)
})

darmanis_gene_ave <- do.call(rbind, darmanis_gene_ave)

darmanis_gene_ave$gene <- darmanis_gene_cells$genes

# NAc
# Use NAc ScRNA ref
pb <- txtProgressBar(min = 0, max = nrow(nac_gene_cells) , width = NA, style = 3)
pb_i <- 0
n_t_idx <- which(cell_mat$Neurons == T)
n_f_idx <- which(cell_mat$Neurons == F)
nac_gene_ave <- apply(data.matrix(nac_gene_cells[, -4170]), 1, function(x){
  
  
  out <- data.table(gene = NA, exp_neu_n = mean(x[n_f_idx]),
                    exp_neu_p  = mean(x[n_t_idx]))
  
  pb_i <<- pb_i+1
  setTxtProgressBar(pb, pb_i, title = paste(round(pb_i/nrow(nac_gene_cells))*100,"% done"))
  
  return(out)
  
})
beepr::beep()
nac_gene_ave <- do.call(rbind, nac_gene_ave)
nac_gene_ave$gene <- nac_gene_cells$Gene


identical(row.names(analysis_dt), darmanis_gene_ave$gene)
identical(row.names(analysis_dt), nac_gene_ave$gene)

dar_m <- data.matrix(darmanis_gene_ave[, 2:3])
dar_nnls <- apply(analysis_dt, 2, function(b){
  ft <- nnls(dar_m, b)
  ft$x/sum(ft$x)
})

dar_nnls <- data.table(t(dar_nnls))
names(dar_nnls) <- c("neu_n", "neu_p")
dar_nnls$samples <- colnames(analysis_dt)


nac_m <- data.matrix(nac_gene_ave[, 2:3])
nac_nnls <- apply(analysis_dt, 2, function(b){
  ft <- nnls(nac_m, b)
  ft$x/sum(ft$x)
})

nac_nnls <- data.table(t(nac_nnls))
names(nac_nnls) <- c("neu_n", "neu_p")
nac_nnls$samples <- colnames(analysis_dt)


pt_dt <- rbind(inner_join(dar_nnls[, .(samples, est = neu_p, ref = "Darmanis")],
                          data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
               inner_join(nac_nnls[, .(samples, est = neu_p, ref = "NAc")],
                          data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm))) %>% data.table

mytable <- pt_dt[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                     RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = ref]

mytable[, `:=`(R_sqrd = round(R_sqrd, 4),
               RMSE = round(RMSE, 4))]



p_save_cols <- RColorBrewer::brewer.pal(3, "Dark2")[1:2]
names(p_save_cols) <- unique(pt_dt$ref)

p_save <- ggplot(data = pt_dt[ref == "Darmanis",], aes(x = est, y = `Houseman DNAm-based`, color = ref)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .8)) +
  labs(x = "MuSiC + Darmanis Cell size", size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  annotation_custom(gridExtra::tableGrob(mytable[ref == "Darmanis", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.08, xmax=0.08, ymin=.38, ymax=.38) +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "darmanis_ref_nnls_v_dna.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")


p_save <- ggplot(data = pt_dt[ref == "NAc",], aes(x = est, y = `Houseman DNAm-based`, color = ref)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
  labs(x = "MuSiC + Darmanis Cell size", size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  annotation_custom(gridExtra::tableGrob(mytable[ref == "NAc", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.08, xmax=0.08, ymin=.38, ymax=.38) +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "nac_ref_nnls_v_dna.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")


# Matrix panel plot data --------------------------------------------------

dar_pd <- rbind(inner_join(music_est[, .(samples, est = `TRUE`, size = "Darmanis")],
                           data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                inner_join(my_music_est_mole[, .(samples, est = `TRUE`, size = "osmFISH totalRNA")],
                           data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                inner_join(my_music_est_size[, .(samples, est = `TRUE`, size = "osmFISH cellsize")],
                           data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                inner_join(my_music_est_nw[, .(samples, est = `TRUE`, size = "NAc")],
                           data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                inner_join(my_music_est_v3[, .(samples, est = `TRUE`, size = "None")],
                           data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm))) %>% data.table
dar_pd$ref <- "Darmanis"



nac_pd <- rbind(inner_join(music_est_v2[, .(samples, est = `TRUE`, size = "NAc")],
                           data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                inner_join(my_music_est_nrna_nac[, .(samples, est = `TRUE`, size = "osmFISH totalRNA")],
                           data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                inner_join(my_music_est_size_nac[, .(samples, est = `TRUE`, size = "osmFISH cellsize")],
                           data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                inner_join(my_music_est_nac_none[, .(samples, est = `TRUE`, size = "None")],
                           data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)),
                inner_join(my_music_est_nw_darmanis[, .(samples, est = `TRUE`, size = "Darmanis")],
                           data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)) %>% data.table) %>% data.table

nac_pd$ref <- "NAc"

all_pd <- rbind(dar_pd, nac_pd)

mytable <- all_pd[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                      RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = .(ref, size)]

mytable[, `:=`(R_sqrd = round(R_sqrd, 4),
               RMSE = round(RMSE, 4))]


p_save_cols <- RColorBrewer::brewer.pal(3, "Dark2")[1:2]
names(p_save_cols) <- unique(all_pd$ref)
p_save_shape <- c(15:19)
names(p_save_shape) <- unique(all_pd$size)

transparent_legend =  theme(
  legend.background = element_rect(fill ="transparent"),
  legend.key = element_rect(fill = "transparent",
                            color = "transparent")
)

remove_grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 2.0)),
  colhead = list(fg_params=list(cex = 2.0)),
  rowhead = list(fg_params=list(cex =2.0)))

# Matrix panel plots ------------------------------------------------------

# Darmanis

p_save <- ggplot(data = all_pd[size == "None" & ref == "Darmanis", ], aes(x = est, y = `Houseman DNAm-based`, color = ref, shape = size)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .8)) +
  xlab("MuSiC + no cell size genes") +
  scale_color_manual(values = p_save_cols, name = "" ) +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable[size == "None" & ref == "Darmanis", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.1, xmax=.1, ymin=.4, ymax=.4) +
  transparent_legend + remove_grid +
  # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"))
ggsave(file.path(".", "model", "manuscript", "dar_no_cell_size_v_dna_V2.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")


# •	B: Music RNA Default with Darmanis
p_save <- ggplot(data = all_pd[size == "Darmanis" & ref == "Darmanis", ], aes(x = est, y = `Houseman DNAm-based`, color = ref, shape = size)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0,.8)) +
  labs(x = "MuSiC estimates", size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "" ) +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable[size == "Darmanis" & ref == "Darmanis", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.65, xmax=.65, ymin=.05, ymax=.05) +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"),
        plot.margin = margin(5,10,10,2))
ggsave(file.path(".", "model", "manuscript", "dar_dna_v_music_V2.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")



p_save <- ggplot(data = all_pd[size == "NAc" & ref == "Darmanis", ], aes(x = est, y = `Houseman DNAm-based`, color = ref, shape = size)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
  xlab("MuSiC + all genes") +
  scale_color_manual(values = p_save_cols, name = "" ) +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable[size == "NAc" & ref == "Darmanis", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.1, xmax=.1, ymin=.4, ymax=.4) +
  transparent_legend + remove_grid +
  # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"))
ggsave(file.path(".", "model", "manuscript", "dar_nacgenes_v_dna_V2.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")




p_save <- ggplot(data = all_pd[size == "osmFISH cellsize" & ref == "Darmanis", ], aes(x = est, y = `Houseman DNAm-based`, color = ref, shape = size)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .8)) +
  labs(x = "MuSiC + osmFISH cell area", size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "" ) +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable[size == "osmFISH cellsize" & ref == "Darmanis", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.15, xmax=.15, ymin=.38, ymax=.38) +
  transparent_legend + remove_grid +
  # ggtitle("Houseman Methylation vs MuSiC \ncell-type proportion estimates\n using osmFISH cellsize\n (Neurons only)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"))
ggsave(file.path(".", "model", "manuscript", "dar_osmFISH_cellsize_v_dna_V2.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")




p_save <- ggplot(data = all_pd[size == "osmFISH totalRNA" & ref == "Darmanis", ], aes(x = est, y = `Houseman DNAm-based`, color = ref, shape = size)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .8)) +
  labs(x = "MuSiC + osmFISH totalRNA", size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "" ) +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable[size == "osmFISH totalRNA" & ref == "Darmanis", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.15, xmax=.15, ymin=.38, ymax=.38) +
  transparent_legend + remove_grid +
  # ggtitle("Houseman Methylation\ vs MuSiC\n cell-type proportion estimates\n using total RNA count as cellsize\n (Neurons only)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"))
ggsave(file.path(".", "model", "manuscript", "dar_osmFISH_mole_v_dna_V2.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")


# NAc

p_save <- ggplot(data = all_pd[size == "osmFISH cellsize" & ref == "NAc", ], aes(x = est, y = `Houseman DNAm-based`, color = ref, shape = size)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
  labs(x = "MuSiC + osmFISH cell area", size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable[size == "osmFISH cellsize" & ref == "NAc", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.1, xmax=.1, ymin=.38, ymax=.38) +
  transparent_legend + remove_grid +
  # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "nac_osmFISH_cellsize_v_dna_V2.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")



p_save <- ggplot(data = all_pd[size == "osmFISH totalRNA" & ref == "NAc", ], aes(x = est, y = `Houseman DNAm-based`, color = ref, shape = size)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
  labs(x = "MuSiC + osmFISH totalRNA", size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable[size == "osmFISH totalRNA" & ref == "NAc", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.1, xmax=.1, ymin=.38, ymax=.38) +
  transparent_legend + remove_grid +
  # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "nac_osmFISH_mole_v_dna_V2.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")



p_save <- ggplot(data = all_pd[size == "Darmanis" & ref == "NAc", ], aes(x = est, y = `Houseman DNAm-based`, color = ref, shape = size)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
  labs(x = "MuSiC + Darmanis Cell size", size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  annotation_custom(gridExtra::tableGrob(mytable[size == "Darmanis" & ref == "NAc", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.08, xmax=0.08, ymin=.38, ymax=.38) +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "nac_ref_darmanis_v_dna_V2.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")



p_save <- ggplot(data = all_pd[size == "NAc" & ref == "NAc", ], aes(x = est, y = `Houseman DNAm-based`, color = ref, shape = size)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
  labs(x = "MuSiC default", size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable[size == "NAc" & ref == "NAc", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.24, xmax=.1, ymin=.38, ymax=.38) +
  transparent_legend + remove_grid +
  # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "nac_ref_v_dna_V2.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")


p_save <- ggplot(data = all_pd[size == "None" & ref == "NAc", ], aes(x = est, y = `Houseman DNAm-based`, color = ref, shape = size)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
  labs(x = "MuSiC no cell size", size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable[size == "None" & ref == "NAc", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.24, xmax=.1, ymin=.38, ymax=.38) +
  transparent_legend + remove_grid +
  # ggtitle("Houseman Methylation vs MuSiC\n cell-type proportion estimates\n(Neurons only)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "nac_ref_no_cell_size_v_dna_V2.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")



# Tables ------------------------------------------------------------------

mytable_nac <- out_pd_v4[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                             RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = type]

mytable_nac[, `:=`(R_sqrd = round(R_sqrd, 4),
                   RMSE = round(RMSE, 4))]
names(mytable_nac)[1] <- "Approach"
mytable_nac <- mytable_nac[order(RMSE)]

mytable_dar <- out_pd_v3[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                             RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = type]

mytable_dar[, `:=`(R_sqrd = round(R_sqrd, 4),
                   RMSE = round(RMSE, 4))]
names(mytable_dar)[1] <- "Approach"
mytable_dar <- mytable_dar[order(RMSE)]


fwrite(rbind(mytable_nac, mytable_dar), 
       file = file.path("cell_type_data", "all_genes.txt"))

# Cibersort ----------------------------------------------------------------


darmanis_ref <- rna_gene_x_cells[genes %in% row.names(analysis_dt), c("genes", cell_type_data$cells), with = F]
darmanis_ref <- darmanis_ref[match(row.names(analysis_dt), genes)]


darmanis_phen <- cell_type_data[match(names(darmanis_ref)[-1], cells), .(cells, Neurons)]
darmanis_phen[, `:=`(N_pos = ifelse(Neurons, 1, 2),
                      N_neg = ifelse(!Neurons, 1, 2))]
identical(darmanis_phen$cells, names(darmanis_ref)[-1])


darmanis_phen_out <- t(darmanis_phen[, 3:4])
row.names(darmanis_phen_out) <- c("Neurons_pos", "Neurons_neg")

names(darmanis_ref)[1] <- "!Gene"

write.table(x = darmanis_ref, file.path("cell_type_data", "darmanis_ref.txt"), quote = F, row.names = F, sep = "\t")
write.table(x = darmanis_phen_out, file.path("cell_type_data", "darmanis_phen.txt"), quote = F, col.names = F, row.names = T, sep = "\t")
write.table(x = data.table("!Gene" = row.names(analysis_dt),analysis_dt), file.path("cell_type_data", "bulk.txt"), quote = F, col.names = T, row.names = F, sep = "\t")



nac_ref <- count_mat_id[Gene %in% row.names(analysis_dt), c("Gene", cell_mat$Cells), with = F]
nac_ref <- nac_ref[match(row.names(analysis_dt), Gene)]


nac_phen <- cell_mat[match(names(nac_ref)[-1], Cells), .(Cells, Neurons)]
nac_phen[, `:=`(N_pos = ifelse(Neurons, 1, 2),
                     N_neg = ifelse(!Neurons, 1, 2))]
identical(nac_phen$Cells, names(nac_ref)[-1])


nac_phen_out <- t(nac_phen[, 3:4])
row.names(nac_phen_out) <- c("Neurons_pos", "Neurons_neg")

write.table(x = nac_ref, file.path("cell_type_data", "nac_ref.txt"), quote = F, row.names = F, sep = "\t")
write.table(x = nac_phen_out, file.path("cell_type_data", "nac_phen.txt"), quote = F, col.names = F, row.names = T, sep = "\t")

# Load results
DNam_darmanis <- data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)
ciber_darmanis <- fread(file.path("cell_type_data", "CIBERSORT.Output_darmanis.txt"))
ciber_nac <- fread(file.path("cell_type_data", "CIBERSORT.Output_nac.txt"))

ciber_darmanis[, type:="Ref:[Darmanis]"]
ciber_nac[, type:="Ref:[Nac]"]


out_ciber <- rbind(inner_join(ciber_darmanis[, .(samples = `Input Sample`, Neurons_pos, type)], DNam_darmanis),
                   inner_join(ciber_nac[, .(samples = `Input Sample`, Neurons_pos, type)], DNam_darmanis)) %>% data.table

mytable_ciber <- out_ciber[, .(R_sqrd = cor(Neurons_pos, `Houseman DNAm-based`)^2,
                      RMSE = sqrt(mean((Neurons_pos -`Houseman DNAm-based`)^2))), by = type]

mytable_ciber[, `:=`(R_sqrd = round(R_sqrd, 4),
               RMSE = round(RMSE, 4))]
names(mytable_ciber)[1] <- "Approach"
mytable_ciber <- mytable_ciber[order(RMSE)]


# Plot
p_save_cols <- RColorBrewer::brewer.pal(3, "Dark2")[1:2]
names(p_save_cols) <- unique(out_ciber$type)

p_save_shape <- c(15:16)
names(p_save_shape) <- unique(out_ciber$type)

p_save_ciber <- ggplot(data = out_ciber, aes(x = Neurons_pos, y = `Houseman DNAm-based`, color = type, shape = type)) + 
  geom_point(size = 3) + 
  geom_abline(slope = 1, intercept = 0) +
  xlab("Cibersort estimates") +
  scale_color_manual(values = p_save_cols, name = "" ) +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable_ciber, rows = NULL), xmin=0.4, xmax=.4, ymin=.01, ymax=.05) +
  transparent_legend + remove_grid +
  ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.85, 0.10),
        legend.title = element_text(face="bold"))
ggsave(file.path(".", "model", "bias_correction_v4_cibersort.png"), plot = p_save_ciber, dpi = "retina", width = 40, height = 30, units = "cm")

set <- "Ref:[Nac]"
set_idx <- which(mytable_ciber$Approach %in% set)
p_save_ciber <- ggplot(data = out_ciber[type %in% set], aes(x = Neurons_pos, y = `Houseman DNAm-based`, color = type, shape = type)) + 
  geom_point(size = 3) + 
  geom_abline(slope = 1, intercept = 0) +
  xlab("Cibersort estimates") +
  scale_color_manual(values = p_save_cols, name = "" ) +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable_ciber[set_idx], rows = NULL), xmin=0.4, xmax=.4, ymin=.01, ymax=.05) +
  transparent_legend + remove_grid +
  ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.85, 0.10),
        legend.title = element_text(face="bold"))
ggsave(file.path(".", "model", "bias_correction_v4_cibersort_nac.png"), plot = p_save_ciber, dpi = "retina", width = 40, height = 30, units = "cm")

set_idx <- which(!mytable_ciber$Approach %in% set)
p_save_ciber <- ggplot(data = out_ciber[!type %in% set], aes(x = Neurons_pos, y = `Houseman DNAm-based`, color = type, shape = type)) + 
  geom_point(size = 3) + 
  geom_abline(slope = 1, intercept = 0) +
  xlab("Cibersort estimates") +
  scale_color_manual(values = p_save_cols, name = "" ) +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable_ciber[set_idx], rows = NULL), xmin=0.4, xmax=.4, ymin=.01, ymax=.05) +
  transparent_legend + remove_grid +
  ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.85, 0.10),
        legend.title = element_text(face="bold"))
ggsave(file.path(".", "model", "bias_correction_v4_cibersort_darmanis.png"), plot = p_save_ciber, dpi = "retina", width = 40, height = 30, units = "cm")


# No quantile normalisation
# Load results
DNam_darmanis <- data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)
ciber_darmanis <- fread(file.path("cell_type_data", "CIBERSORT.Output_darmanis_no_quantile.txt"))
ciber_nac <- fread(file.path("cell_type_data", "CIBERSORT.Output_nac_no_quantile.txt"))

ciber_darmanis[, type:="Ref:[Darmanis]; No quantile norm"]
ciber_nac[, type:="Ref:[Nac]; No quantile norm"]


out_ciber <- rbind(inner_join(ciber_darmanis[, .(samples = `Input Sample`, Neurons_pos, type)], DNam_darmanis),
                   inner_join(ciber_nac[, .(samples = `Input Sample`, Neurons_pos, type)], DNam_darmanis)) %>% data.table

mytable_ciber <- out_ciber[, .(R_sqrd = cor(Neurons_pos, `Houseman DNAm-based`)^2,
                               RMSE = sqrt(mean((Neurons_pos -`Houseman DNAm-based`)^2))), by = type]

mytable_ciber[, `:=`(R_sqrd = round(R_sqrd, 4),
                     RMSE = round(RMSE, 4))]
names(mytable_ciber)[1] <- "Approach"
mytable_ciber <- mytable_ciber[order(RMSE)]


# Plot
p_save_cols <- RColorBrewer::brewer.pal(3, "Dark2")[1:2]
names(p_save_cols) <- unique(out_ciber$type)

p_save_shape <- c(15:16)
names(p_save_shape) <- unique(out_ciber$type)

p_save_ciber <- ggplot(data = out_ciber, aes(x = Neurons_pos, y = `Houseman DNAm-based`, color = type, shape = type)) + 
  geom_point(size = 3) + 
  geom_abline(slope = 1, intercept = 0) +
  xlab("Cibersort estimates") +
  scale_color_manual(values = p_save_cols, name = "" ) +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable_ciber, rows = NULL), xmin=0.4, xmax=.4, ymin=.01, ymax=.05) +
  transparent_legend + remove_grid +
  ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.85, 0.10),
        legend.title = element_text(face="bold"))
ggsave(file.path(".", "model", "bias_correction_v4_cibersort_nq.png"), plot = p_save_ciber, dpi = "retina", width = 40, height = 30, units = "cm")

set <- "Ref:[Nac]; No quantile norm"
set_idx <- which(mytable_ciber$Approach %in% set)
p_save_ciber <- ggplot(data = out_ciber[type %in% set], aes(x = Neurons_pos, y = `Houseman DNAm-based`, color = type, shape = type)) + 
  geom_point(size = 3) + 
  geom_abline(slope = 1, intercept = 0) +
  xlab("Cibersort estimates") +
  scale_color_manual(values = p_save_cols, name = "" ) +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable_ciber[set_idx], rows = NULL), xmin=0.4, xmax=.4, ymin=.01, ymax=.05) +
  transparent_legend + remove_grid +
  ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.85, 0.10),
        legend.title = element_text(face="bold"))
ggsave(file.path(".", "model", "bias_correction_v4_cibersort_nac_nq.png"), plot = p_save_ciber, dpi = "retina", width = 40, height = 30, units = "cm")

set_idx <- which(!mytable_ciber$Approach %in% set)
p_save_ciber <- ggplot(data = out_ciber[!type %in% set], aes(x = Neurons_pos, y = `Houseman DNAm-based`, color = type, shape = type)) + 
  geom_point(size = 3) + 
  geom_abline(slope = 1, intercept = 0) +
  xlab("Cibersort estimates") +
  scale_color_manual(values = p_save_cols, name = "" ) +
  scale_shape_manual(values = p_save_shape, name = "" ) +
  annotation_custom(gridExtra::tableGrob(mytable_ciber[set_idx], rows = NULL), xmin=0.4, xmax=.4, ymin=.01, ymax=.05) +
  transparent_legend + remove_grid +
  ggtitle("Relationship between Houseman Methylation and cell-type proportion estimates\n(Neurons only)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.85, 0.10),
        legend.title = element_text(face="bold"))
ggsave(file.path(".", "model", "bias_correction_v4_cibersort_darmanis_nq.png"), plot = p_save_ciber, dpi = "retina", width = 40, height = 30, units = "cm")





# Misc --------------------------------------------------------------------

test_dt <- count_mat_id[Gene %in% row.names(analysis_dt),]
test_dt <- data.table(Cells = names(test_dt)[-1],
                      tot_rna = apply(test_dt[, -1], 2, sum ))

test_dt <- inner_join(test_dt, cell_mat) %>% data.table


darmanis_sizes <- data.table(darmanis_sizes)
NAc_all_genes <- test_dt[, .(N = mean(tot_rna)), by = Neurons]

all_sizes_r <- rbind(NAc_sizes,
                     darmanis_sizes[, .(Neurons = ct, size = cs, type = "Darmanis")],
                     osmFISH_sizes[, .(Neurons, size, type = "osmFISH cell Area")],
                     osmFISH_sizes[, .(Neurons, size = total_mole, type = "osmFISH nRNA")],
                     NAc_all_genes[, .(Neurons, size = N, type = "NAc all genes")])

all_sizes_c <- inner_join(NAc_sizes[type == 50, .(Neurons, "NAc 50 genes" = size)],
                          NAc_sizes[type == 25, .(Neurons, "NAc 25 genes" = size)]) %>% 
  inner_join(NAc_all_genes[, .(Neurons, "NAc all genes" = N)]) %>% 
  inner_join(darmanis_sizes[, .(Neurons = as.logical(ct)  , Darmanis = cs)]) %>% 
  inner_join(osmFISH_sizes[, .(Neurons, "osmFISH cell Area" = size)]) %>% 
  inner_join(osmFISH_sizes[, .(Neurons, "osmFISH nRNA" = total_mole)])

fwrite(all_sizes_c, file = "all_sizes_c.txt")

t.test( tot_rna ~ Neurons,data = test_dt)
wilcox.test(tot_rna ~ Neurons,data = test_dt)


# Comparing Gene expression profiles

dar_ttests <- lapply(1:nrow(darmanis_gene_cells), function(i){
  x <- darmanis_gene_cells[i, -266]
  nm <- names(x)
  
  test_dt <- inner_join(data.table(cells = nm, "rna" = as.numeric(t(x))),
                        cell_type_data[, .(cells, Neurons )], by = "cells") %>% data.table
  
  res <- t.test( rna ~ Neurons,data = test_dt)
  
  data.table(gene = darmanis_gene_cells[i, genes], t_stat = res$statistic, darmanis_pval = res$p.value)
})
dar_ttests <- do.call(rbind, dar_ttests)
names(dar_ttests)[2] <- "dar_t_stat"


# nac_ttests <- lapply(1:nrow(nac_gene_cells), function(i){
#   x <- nac_gene_cells[i, -4170]
#   nm <- names(x)
#   
#   test_dt <- inner_join(data.table(Cells = nm, "rna" = as.numeric(t(x))),
#                         cell_mat[, .(Cells, Neurons )], by = "Cells") %>% data.table
#   
#   res <- t.test( rna ~ Neurons, data = test_dt)
#   
#   data.table(gene = nac_gene_cells[i, Gene], t_stat = res$statistic, darmanis_pval = res$p.value)
# })
# nac_ttests <- do.call(rbind, nac_ttests)
# names(nac_ttests)[2] <- "nac_t_stat"

# For all Genes
ncpus <- 4
cl <- makeCluster(ncpus, outfile="")
registerDoParallel(cl)
clusterExport(cl, c("nac_gene_cells", "cell_mat"))
# clusterExport(cl, ls(), envir = .GlobalEnv)

nac_ttests <- parLapply(cl, 1:nrow(nac_gene_cells), function(i){
  require(data.table)
  require(dplyr)
  x <- nac_gene_cells[i, -4170]
  nm <- names(x)
  
  test_dt <- inner_join(data.table(Cells = nm, "rna" = as.numeric(t(x))),
                        cell_mat[, .(Cells, Neurons )], by = "Cells") %>% data.table
  
  res <- t.test( rna ~ Neurons, data = test_dt)
  
  data.table(gene = nac_gene_cells[i, Gene], t_stat = res$statistic, nac_pval = res$p.value)
})
beepr::beep()
nac_ttests <- do.call(rbind, nac_ttests)
names(nac_ttests)[2] <- "nac_t_stat"
stopCluster(cl)


all_ttests <- inner_join(dar_ttests, nac_ttests) %>%  data.table

# fwrite(all_ttests, file.path("all_ttests.txt"))
# fread("all_ttests.txt") -> all_ttests

transparent_legend =  theme(
  legend.background = element_rect(fill ="transparent"),
  legend.key = element_rect(fill = "transparent",
                            color = "transparent")
)

remove_grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))

with(all_ttests[is.finite(nac_t_stat) & is.finite(dar_t_stat)], round((cor(nac_t_stat, dar_t_stat)), 4))

p_save_ttests <- ggplot(data = all_ttests, aes(x = nac_t_stat, y = dar_t_stat)) + 
  geom_point(size = 4) + 
  # geom_abline(slope = 1, intercept = 0, size = 3) +
  geom_smooth(method = "auto", size = 3) +
  xlab("NAc T-statistics") +
  ylab("Darmanis T-statistics") +
  scale_x_continuous(breaks = seq(-60,60, by = 10), limits = c(-55,55)) +
  scale_y_continuous(breaks = seq(-60,60, by = 10), limits = c(-55,55)) +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 14),
        axis.text.y=element_text(size = 12, face="bold"),
        axis.text.x=element_text(size = 12, face="bold"),
        legend.position = c(0.85, 0.10),
        legend.title = element_text(face="bold"))
ggsave(file.path(".", "model", "manuscript", "nac_v_dar_ttests.pdf"), plot = p_save_ttests, dpi = "retina", width = 20, height = 20, units = "cm")


with(all_ttests[is.finite(nac_t_stat) & is.finite(dar_t_stat)], qqplot(nac_t_stat, dar_t_stat))
abline(a = 0, b = 1)

plt_dt <- rbind(all_ttests[, .(gene, t_stat = dar_t_stat, data = "Darmanis")],
                all_ttests[, .(gene, t_stat = nac_t_stat, data = "NAc")])

p_save_ttests <- ggplot(data = plt_dt, aes(x = data, y = abs(t_stat))) + 
  geom_violin(fill = "skyblue") +
  geom_boxplot(width = .05) +
  xlab("Reference dataset") +
  ylab(expression(abs("T-statistics"))) +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 14),
        axis.text.y=element_text(size = 12, face="bold"),
        axis.text.x=element_text(size = 12, face="bold"),
        legend.position = c(0.85, 0.10),
        legend.title = element_text(face="bold"))
ggsave(file.path(".", "model", "manuscript", "nac_v_dar_ttests_bxplt.pdf"), plot = p_save_ttests, dpi = "retina", width = 20, height = 20, units = "cm")



# Test the count of cells -------------------------------------------------

# Use NAc ScRNA ref

sc_data <- count_mat_id[Gene %in% row.names(analysis_dt),]
sc_data <- inner_join(data.table(Gene = row.names(analysis_dt)), sc_data) %>% data.table
rn <- sc_data$Gene
cn <- colnames(sc_data)
bulk_data <- ExpressionSet(assayData = analysis_dt)

M <- 1e3
k <- 265


exp_nac <- apply(sc_data[, -1], 2, sum)
exp_dar <- apply(darmanis_gene_cells[, -266], 2, sum)

# my_cn_sample <- matrix(0, ncol = M, nrow = k)

# p_dt_srt <- p_dt[order(Neurons)]
# cell_mat2 <- inner_join(p_dt_srt[, .(Cells)], cell_mat) %>% data.table


# s1 <- sampling::strata(p_dt_srt, stratanames = "Neurons", size = c(134, 131), pik = p_dt$p, method = "srswor") %>% data.table
# s1[, .N/265, by = Neurons]

# my_cn_sample[, 1] <- p_dt[s1$ID_unit, Cells]

cell_mat[, idx := 1:.N]
my_cn_sample <- replicate(M, c(sample(cell_mat[Neurons == F, idx], 134),
                                 sample(cell_mat[Neurons == T, idx], 131)) )
my_cn_sample <- apply(my_cn_sample, 2, sort)
my_cn_sample <- unique(my_cn_sample, MARGIN = 2)

my_cn_sample_cells <- matrix(0, ncol = M, nrow = k)


pb <- txtProgressBar(min = 0, max = M , width = NA, style = 3)
pb_i <- 0
music_est_nac_sample_size_2 <- lapply(1:M, function(r){
  
  if(ncol(my_cn_sample) != M){
    stop("Unique columns less than desired size M")
  }
  
  my_cn_sample_2 <- cell_mat[my_cn_sample[, r], .(Cells, Neurons)]
  my_cn_sample_cells[, r] <<- my_cn_sample_2$Cells
  my_sc_data <- sc_data[ , my_cn_sample_2$Cells, with = F]
  my_sc_data <- data.matrix(my_sc_data)
  row.names(my_sc_data) <- rn
  my_sc_data <- ExpressionSet(assayData = my_sc_data)
  
  
  music_est_v2 <- suppressMessages(music_prop(bulk.eset = bulk_data, sc.eset = my_sc_data, clusters = my_cn_sample_2$Neurons, samples = my_cn_sample_2$Cells))
  music_est_v2 <- data.table(samples  = row.names(music_est_v2$Est.prop.allgene), music_est_v2$Est.prop.weighted,
                             Iteration = r) 
  
  pb_i <<- pb_i+1
  setTxtProgressBar(pb, pb_i, title = paste(round(pb_i/M)*100,"% done"))
  
  return(music_est_v2)
  
})

music_est_nac_sample_size_2 <- do.call(rbind, music_est_nac_sample_size_2)

fwrite(music_est_nac_sample_size_2, "music_est_nac_sample_size_2.txt")
fwrite(my_cn_sample_cells, "my_cn_sample_cells.txt")

k <- 1e3
music_est_nac_sample_size <- rbind(music_est_nac_sample_size_2[Iteration %in% 1:k],
                                   music_est_v2[, .(samples, `FALSE`, `TRUE`, Iteration = "all" )])

cell_type_data[, .N/265, by = Neurons]
cell_mat[, .N/4169, by = Neurons]
osmfish_test_dt[, .N/3953, by = Neurons]

# Plot
load(file.path("~", "cell_type_data", "NAc_rse_gene_withCompEsts_update.rda"))
pt_dt <- inner_join(music_est_nac_sample_size[, .(samples, est = `TRUE`, type = Iteration)],
                    data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)) %>% data.table

mytable <- pt_dt[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                     RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = type]

mytable[, `:=`(R_sqrd = round(R_sqrd, 4),
               RMSE = round(RMSE, 4))]



p_save_cols <- c(rep("gray", k), "black")
names(p_save_cols) <- unique(pt_dt$type)
transparent_legend =  theme(
  legend.background = element_rect(fill ="transparent"),
  legend.key = element_rect(fill = "transparent",
                            color = "transparent")
)

remove_grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))

no_x_axis_label <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

p_save <- ggplot(data = pt_dt, aes(x = est, y = `Houseman DNAm-based`, color = type)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, .5)) +
  labs(x = "MuSiC:Darmanis default", size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  # annotation_custom(gridExtra::tableGrob(mytable[type == "all", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.08, xmax=0.08, ymin=.38, ymax=.38) +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "cell_sample_size.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")


p_save <- ggplot(data = mytable, aes(x = R_sqrd, y = RMSE, color = type)) + 
  geom_point(size = 10) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0,1, by = .25), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0,.25, by = 0.05), limits = c(0, .25)) +
  labs(x = expression("R"^2), size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  # annotation_custom(gridExtra::tableGrob(mytable[type == "all", .(R_sqrd, RMSE)], rows = NULL, theme = mytheme), xmin=0.08, xmax=0.08, ymin=.38, ymax=.38) +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "r_sqrd_v_rmse.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")



# Proportion of expressed cells per cell type across genes (per dataset) --------

pb <- txtProgressBar(min = 0, max = ncol(my_cn_sample_cells) , width = NA, style = 3)
pb_i <- 0

nac_gene_ave_dt <- apply(my_cn_sample_cells, 2, function(my_cells){
  
  my_cell_mat <- cell_mat[Cells %in% my_cells]
  n_t_idx <- which(my_cell_mat$Neurons == T)
  n_f_idx <- which(my_cell_mat$Neurons == F)
  
  
  nac_gene_ave <- apply(data.matrix(nac_gene_cells[, my_cells, with = F]), 1, function(x){
    
    
    out <- data.table(gene = NA, prop_exp_neu_n = mean(x[n_f_idx] > 0),
                      prop_exp_neu_p  = mean(x[n_t_idx] > 0))
    
    return(out)
    
  })
  
  nac_gene_ave <- do.call(rbind, nac_gene_ave)
  
  nac_gene_ave <- nac_gene_ave[, .(sc_ave_prop_n = mean(prop_exp_neu_n)/sd(prop_exp_neu_n),
                                   sc_ave_prop_p = mean(prop_exp_neu_p)/sd(prop_exp_neu_p),
                                   ave_prop_n = mean(prop_exp_neu_n),
                                   ave_prop_ = mean(prop_exp_neu_p))]
  
  pb_i <<- pb_i+1
  setTxtProgressBar(pb, pb_i, title = paste(round(pb_i/nrow(nac_gene_cells))*100,"% done"))
  
  return(nac_gene_ave)
})


nac_gene_ave_dt <- do.call(rbind, nac_gene_ave_dt)

fwrite(nac_gene_ave_dt, "nac_gene_ave_dt.txt")




# Cell size per dataset ---------------------------------------------------

pb <- txtProgressBar(min = 0, max = ncol(my_cn_sample_cells) , width = NA, style = 3)
pb_i <- 0

nac_cs_dt <- apply(my_cn_sample_cells, 2, function(my_cells){
  
  my_cell_mat <- cell_mat[Cells %in% my_cells]
  n_t_idx <- which(my_cell_mat$Neurons == T)
  n_f_idx <- which(my_cell_mat$Neurons == F)
  
  
  nac_cs <- apply(data.matrix(nac_gene_cells[, my_cells, with = F]), 2, sum)
  
  out <- data.table(cell_size_f = mean(nac_cs[n_f_idx]),
                             cell_size_t = mean(nac_cs[n_t_idx]))
  
  pb_i <<- pb_i+1
  setTxtProgressBar(pb, pb_i, title = paste(round(pb_i/nrow(nac_gene_cells))*100,"% done"))
  
  return(out)
})


nac_cs_dt <- do.call(rbind, nac_cs_dt)

fwrite(nac_cs_dt, "nac_cs_dt.txt")


# Proportion of expressed cells per cell type across genes (all NAc cells) ----------------





pb <- txtProgressBar(min = 0, max = nrow(darmanis_gene_cells) , width = NA, style = 3)
pb_i <- 0

n_t_idx <- which(cell_type_data$Neurons == T)
n_f_idx <- which(cell_type_data$Neurons == F)
darmanis_gene_ave <- apply(data.matrix(darmanis_gene_cells[, -266]), 1 ,function(x){
  
  out <- data.table(gene = NA, prop_exp_neu_n = mean(x[n_f_idx] > 0),
                    prop_exp_neu_p  = mean(x[n_t_idx] > 0))
  
  pb_i <<- pb_i+1
  setTxtProgressBar(pb, pb_i, title = paste(round(pb_i/nrow(darmanis_gene_cells))*100,"% done"))
  
  return(out)
})

darmanis_gene_ave <- do.call(rbind, darmanis_gene_ave)

darmanis_gene_ave$gene <- darmanis_gene_cells$genes

pb <- txtProgressBar(min = 0, max = nrow(nac_gene_cells) , width = NA, style = 3)
pb_i <- 0
n_t_idx <- which(cell_mat$Neurons == T)
n_f_idx <- which(cell_mat$Neurons == F)
nac_gene_ave <- apply(data.matrix(nac_gene_cells[, -4170]), 1, function(x){
  
  
  out <- data.table(gene = NA, prop_exp_neu_n = mean(x[n_f_idx] > 0),
                    prop_exp_neu_p  = mean(x[n_t_idx] > 0))
  
  pb_i <<- pb_i+1
  setTxtProgressBar(pb, pb_i, title = paste(round(pb_i/nrow(nac_gene_cells))*100,"% done"))
  
  return(out)
  
})
beepr::beep()

nac_gene_ave <- do.call(rbind, nac_gene_ave)

nac_gene_ave$gene <- nac_gene_cells$Gene


plt_dt <- rbind(darmanis_gene_ave[, .(gene, prop = prop_exp_neu_n, type = "non-neuronal", data = "Darmanis")],
                darmanis_gene_ave[, .(gene, prop = prop_exp_neu_p, type = "neuronal", data = "Darmanis")],
                nac_gene_ave[, .(gene, prop = prop_exp_neu_n, type = "non-neuronal", data = "NAc")],
                nac_gene_ave[, .(gene, prop = prop_exp_neu_p, type = "neuronal", data = "NAc")])

t.test(prop ~ data, data = plt_dt[type == "non-neuronal"])
t.test(prop ~ data, data = plt_dt[type == "neuronal"])

wilcox.test(prop ~ data, data = plt_dt[type == "non-neuronal"])
wilcox.test(prop ~ data, data = plt_dt[type == "neuronal"])


p_save_cols <- RColorBrewer::brewer.pal(3, "Dark2")[1:2]
names(p_save_cols) <- unique(plt_dt$type)

p_save <- ggplot(data = plt_dt, aes(x = data, y = prop, fill = type)) + 
  geom_boxplot() + 
  labs(y = "Proportion of expressed cells" ,x = "Reference dataset", size = rel(1.5)) +
  scale_fill_manual(values = p_save_cols, name = "Cell type") +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        # legend.position = "topright",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "prop_cell_gene.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")



# Plot --------------------------------------------------------------------

darmanis_gene_ave[, .(sc_ave_prop_n = mean(prop_exp_neu_n)/sd(prop_exp_neu_n),
                      sc_ave_prop_p = mean(prop_exp_neu_p)/sd(prop_exp_neu_p))]

nac_gene_ave[, .(sc_ave_prop_n = mean(prop_exp_neu_n)/sd(prop_exp_neu_n),
                 sc_ave_prop_p = mean(prop_exp_neu_p)/sd(prop_exp_neu_p))]

darmanis_gene_ave[, .(mean(prop_exp_neu_n), mean(prop_exp_neu_p))]

plt_dt <- rbind(nac_gene_ave_dt[, .(sc_ave_prop_n, sc_ave_prop_p, Iteration = 1:.N)],
                nac_gene_ave[, .(sc_ave_prop_n = mean(prop_exp_neu_n)/sd(prop_exp_neu_n),
                                 sc_ave_prop_p = mean(prop_exp_neu_p)/sd(prop_exp_neu_p),
                                 Iteration = "all")])

plt_dt <- rbind(plt_dt[, .(est = sc_ave_prop_n, Iteration, type = "non-neuronal") ],
                plt_dt[, .(est = sc_ave_prop_p, Iteration, type = "neuronal") ])

vline.dat <- plt_dt[Iteration == "all"] 

p_save_cols <- c(rep("gray", 1e3), "black")
names(p_save_cols) <- unique(plt_dt$Iteration)

p_save <- ggplot(data = plt_dt, aes(x = est, fill = type)) + 
  geom_histogram() + 
  labs(x = "Scaled proportion of expressed cells" ,x = "", size = rel(1.5)) +
  facet_grid( .~type) +
  geom_vline(aes(xintercept=est), data=vline.dat) + 
  scale_fill_discrete(name = "Cell type") +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = c(.8, .5),
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "prop_cell_gene_dt.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")


n_t_idx <- which(cell_mat$Neurons == T)
n_f_idx <- which(cell_mat$Neurons == F)


nac_cs <- apply(data.matrix(nac_gene_cells[, cell_mat$Cells, with = F]), 2, sum)

nac_cs <- data.table(cell_size_f = mean(nac_cs[n_f_idx]),
                  cell_size_t = mean(nac_cs[n_t_idx]))



plt_dt <- rbind(nac_cs_dt[, .(cell_size_f, cell_size_t, Iteration = 1:.N)],
                nac_cs[, .(cell_size_f, cell_size_t, Iteration = "all")])

p_save <- ggplot(data = plt_dt, aes(x = cell_size_t, y = cell_size_f, color = Iteration)) + 
  geom_point(size = 10) + 
  labs(y = "non-neuronal cells" ,x = "neuronal cells", size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "cell_size_dt.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")



nac_gene_ave_dt <- rbind(nac_gene_ave_dt[, .(ave_prop_n, ave_prop_p = ave_prop_)],
                         nac_gene_ave[, .(ave_prop_n = mean(prop_exp_neu_n), ave_prop_p = mean(prop_exp_neu_p))])
nac_cs_dt <- rbind(nac_cs_dt, nac_cs)
plt_dt <- cbind(mytable, nac_gene_ave_dt, nac_cs_dt)

# plt_dt <- rbind(plt_dt[, .(type, R_sqrd, RMSE, ave_prop_n, ave_prop_p, grp = "gen")],
#                 data.table(type ="Darmanis", R_sqrd = 0.5548, RMSE = 0.3653, 
#                            darmanis_gene_ave[, .(ave_prop_n = mean(prop_exp_neu_n), 
#                                                  ave_prop_p = mean(prop_exp_neu_p))],
#                            grp = "Darmanis"))
# 
# plt_dt$grp[plt_dt$type == "all"] <- "NAc all cells"
# 
# # p_save_cols <- c(rep("gray", 1e3), RColorBrewer::brewer.pal(3, "Dark2")[1:2])
# # names(p_save_cols) <- unique(plt_dt$type)
# # col_lbl <- c(rep(" ", 1e3), "NAc all", "Darmanis")
# 
# p_save_cols <- c("gray", RColorBrewer::brewer.pal(3, "Dark2")[1:2])
# names(p_save_cols) <- unique(plt_dt$grp)
# col_lbl <- c("Generated datasets", "NAc all cells", "Darmanis")
# 



p_save <- ggplot(data = plt_dt, aes(x = ave_prop_p, y = R_sqrd, color = type)) + 
  geom_point(size = 10) + 
  labs(x = "Proportion of expressed \n neuronal cells" , y = expression("R"^2), size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "n_v_r_sqrd.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")





p_save <- ggplot(data = plt_dt, aes(x = ave_prop_p, y = 1 - RMSE, color = type)) + 
  geom_point(size = 10) + 
  labs(x = "Proportion of expressed \n neuronal cells", y = expression("1 - RMSE"), size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "n_v_rmse.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")



p_save <- ggplot(data = plt_dt, aes(y = cell_size_t, x = R_sqrd, color = type)) + 
  geom_point(size = 10) + 
  labs(y = "Cell size" ,x = expression("RMSE"), size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "cs_v_r_sqrd.pdf"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")


summary(lm(scale(R_sqrd) ~ -1 + scale(ave_prop_p), data = plt_dt[-1001,]))
summary(lm(scale(1-RMSE) ~ -1 + scale(ave_prop_p), data = plt_dt[-1001,]))
summary(lm(scale(1-RMSE) ~ -1 + scale(ave_prop_p), data = plt_dt[-1001,]))$coefficients

fread("nac_gene_ave_dt.txt") -> nac_gene_ave_dt
dt <- cbind(mytable[-1001,], nac_gene_ave_dt)


summary(lm(scale(R_sqrd) ~ -1 + scale(sc_ave_prop_p), data = dt))
summary(lm(scale(RMSE) ~ -1 + scale(sc_ave_prop_p), data = dt))


plt_dt <- rbind(nac_gene_ave_dt[, .(sc_ave_prop_n, sc_ave_prop_p, Iteration = 1:.N)],
                nac_gene_ave[, .(sc_ave_prop_n = mean(prop_exp_neu_n)/sd(prop_exp_neu_n),
                                 sc_ave_prop_p = mean(prop_exp_neu_p)/sd(prop_exp_neu_p),
                                 Iteration = "all")])

plt_dt <- cbind(mytable, plt_dt)

p_save <- ggplot(data = plt_dt, aes(y = sc_ave_prop_p, x = R_sqrd, color = type)) + 
  geom_point(size = 10) + 
  labs(y = " Scaled proportion of expressed \n neuronal cells" ,x = expression("R"^2), size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "scn_v_r_sqrd.png"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")





p_save <- ggplot(data = plt_dt, aes(y = sc_ave_prop_p, x = RMSE, color = type)) + 
  geom_point(size = 10) + 
  labs(y = "Scaled proportion of expressed \n neuronal cells" ,x = expression("RMSE"), size = rel(1.5)) +
  scale_color_manual(values = p_save_cols, name = "") +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 20),
        axis.title = element_text(face="bold", size = 30),
        axis.text.y=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        axis.text.x=element_text(size = rel(1.3), face="bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(face="bold"), 
        legend.text = element_text(size = 12, face="bold"))
ggsave(file.path(".", "model", "manuscript", "scn_v_rmse.png"), plot = p_save, dpi = "retina", width = 20, height = 20, units = "cm")

