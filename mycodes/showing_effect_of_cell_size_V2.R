# Preamble ----------------------------------------------------------------

packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer", 
           "mixtools","matrixStats", "MuSiC", "convert", "xbioc", "ggplot2", "sva", "plotly")
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




# MuSiC + Darmanis -------------------------------------------------------------------


analysis_dt <- which(row.names(yExprs) %in% gg)
analysis_dt <- yExprs[analysis_dt, ]


sc_data <- rna_gene_x_cells[genes %in% row.names(analysis_dt),]
sc_data <- inner_join(data.table(genes = row.names(analysis_dt)), sc_data) %>% data.table  #sc_data[match(row.names(analysis_dt), genes),]
rn <- sc_data$genes
cn <- colnames(sc_data)
cn <- inner_join(data.table(cells = cell_type_data$cells), data.table(cells = cn)) %>% data.table #cn <- cn[cn %in% cell_type_data$cells]
sc_data <- sc_data[ , cn$cells, with = F]
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
sc_data <- sc_data[ , cn$cells, with = F]
sc_data <- data.matrix(sc_data)
row.names(sc_data) <- rn
sc_data <- ExpressionSet(assayData = sc_data)
bulk_data <- ExpressionSet(assayData = analysis_dt)



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
  
  par(mar=c(5,6,6,2), cex.axis=2,cex.lab=2))


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

# Based on > 146 genes
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


