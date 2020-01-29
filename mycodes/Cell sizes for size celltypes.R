
packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer", 
           "mixtools","matrixStats", "MuSiC", "convert", "xbioc", "ggplot2", "sva", "plotly",
           "doParallel", "parallel")
libs_loaded <- sapply(packs, library, character.only = T)

type_anal <- "Celltype" # or Celltype

scrna <- "single" # both

# Get common set of Genes -------------------------------------------------------------------

rna_gene_x_cells <- fread(file.path(".", "cell_data", "rna_gene_x_cells.txt"))
cell_type_data <- fread(file.path(".", "cell_data", "cell_type_info.txt"))
cell_type_data <- cell_type_data[ all %in% c("Neurons", "OPC",  "Astrocytes", "Endothelial",
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




# Get Cell size for 6 cell types ------------------------------------------


# Music
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
analysis_dt <- yExprs[analysis_dt, ]


# Mouse

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


osmfish_test_dt <- inner_join(my_cell_atts, cell_types_aj[, .(ClusterName, Celltype = Class, Subclass)]) %>% data.table 

osmfihs_ests <- my_cell_atts[, .(mean(CellArea), .N, var(CellArea)), by = .(ClusterName, valid)]

osmfihs_ests <- inner_join(osmfihs_ests, cell_types_aj[, .(ClusterName, Celltype = Class, Subclass)]) %>% data.table


my_ests <- osmfihs_ests[, .("w_mean" = weighted.mean(V1, N)), by = Celltype]
my_ests[, rel_size := w_mean/sum(w_mean)]



pd_data <- inner_join(my_ests[, .(Celltype, "osmFISH_est" = w_mean)], 
                      data.table("music_est" = my_s_k, Celltype = c("Astrocyte", "Endothelial", "Microglia", "Neuron", "Oligodendrocyte", "OPC"))) %>% data.table

pd_data_unique <- pd_data[, .(osmFISH_est, Celltype, music_est)] %>% unique




plot_dat_oshm <- my_cell_atts
plot_dat_oshm[, Celltype := NULL]
plot_dat_oshm <- inner_join(plot_dat_oshm, cell_types_aj[, .(ClusterName, Celltype = Class, Subclass)]) %>% data.table


rna_ests <- plot_dat_oshm[, .(mean(totalmolecules), .N, Celltype = unique(Celltype)), by = .(ClusterName, valid)]
rna_ests <- rna_ests[, .("w_mean" = weighted.mean(V1, N)), by = Celltype]
rna_ests[, rel_abun := w_mean/sum(w_mean)]

w_mole_count <- rna_ests$w_mean/my_ests$w_mean
w_mole_count/sum(w_mole_count)

out <- inner_join(rna_ests[, .(Celltype, "Rel_abun RNA in mouse" = rel_abun, total_mole = w_mean)], my_ests[, .(Celltype, "Rel_size in mouse" = rel_size, size = w_mean)]) %>% 
  inner_join(data.table("music_est" = my_s_k/sum(my_s_k), Celltype = c("Astrocyte", "Endothelial", "Microglia", "Neuron", "Oligodendrocyte", "OPC"))) %>% 
  inner_join(data.table(Celltype = rna_ests$Celltype, scaled_rna = w_mole_count/sum(w_mole_count))) %>% data.table



fwrite(out, "six_celltypes_sizes.txt")
