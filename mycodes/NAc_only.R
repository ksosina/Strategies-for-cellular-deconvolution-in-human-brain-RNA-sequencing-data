# Preamble ----------------------------------------------------------------

packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer", 
           "mixtools","matrixStats", "MuSiC", "convert", "xbioc", "ggplot2", "sva", "plotly")
libs_loaded <- sapply(packs, library, character.only = T)



# Get common set of Genes -------------------------------------------------------------------


load("./NAc_Nicotine_hg38_rseGene_rawCounts_postQCSamples_n223.rda")
load("./MatchedCellComp/singleCell_iPSC_quake_coefEsts_calibration_Zscale_adultOnly.rda")
yExprs <- log2(getRPKM(rse_gene, "Length")+1)





# Gene x Cell matrix
count_mat <- fread(file.path("cell_type_data", "countMatrix_n4169-NAc-nuclei.csv"))
cell_mat <- fread(file.path("cell_type_data", "pd-cellTypeAssignment_n4169.csv"))
cell_mat[, Neurons:= nucleusCellType == "Neuron"]
names(count_mat)[1] <- "Gene"
names(cell_mat)[1] <- "Cells"

my_out_cell_mat <- cell_mat[, .(size = mean(nCount_RNA)), by = Neurons]








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



# Missing genes
setdiff(count_mat$Gene, id_dt$Gene) %>% data.table

gg <- rownames(yExprs)
count_mat_id <- inner_join(id_dt[, .(gencodeID, Gene)], count_mat) %>%  data.table
count_mat_id <- inner_join(data.table( gencodeID = gg), count_mat_id) %>% data.table
count_mat_id[, Gene:=gencodeID]
count_mat_id[, `:=`(gencodeID = NULL)]



# Gene count will be different!
gg <- count_mat_id$Gene





# MuSiC + NAc -------------------------------------------------------------------


analysis_dt <- which(row.names(yExprs) %in% gg)
analysis_dt <- yExprs[analysis_dt, ]

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

# Note this is the same cell size estimated by MuSiC. Feel free to change this
cell_size_data <- data.frame(my_out_cell_mat[, .(Celltype = Neurons, w_mean = size )])

my_music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_mat$Neurons, samples = cell_mat$Cells,
                                 cell_size = cell_size_data)