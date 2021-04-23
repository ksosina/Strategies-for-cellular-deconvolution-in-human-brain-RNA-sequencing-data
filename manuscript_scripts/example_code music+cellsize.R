# Preamble ----------------------------------------------------------------

packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer", 
           "mixtools","matrixStats", "MuSiC", "convert", "xbioc", "ggplot2", "sva", "plotly",
           "doParallel", "parallel")
libs_loaded <- sapply(packs, library, character.only = T)


# Load data -------------------------------------------------------------------

# Bulk RNA-seq
bulk_rna <- fread(file.path(".", "toy_bulk.txt"))
# sn-RNA data. Rows are genes and columns are nuclei
ref_data <-  fread(file.path(".", "nac_ref.txt"))
# Cell sizes per cell type (Neuron vs glial)
cell_size <-  fread(file.path(".", "nac_size.txt"))
# Cell type  per cell in ref_data (Neuron vs glial)
cell_type_data <-  fread(file.path(".", "nac_cell_class.txt"))



# MuSiC + cell size -------------------------------------------------------------------

# Reference and Bulk datasets have to be matrices then converted to ExpressionSet classes
# Both datasets have to have the same row names

sc_data <- data.matrix(ref_data[, -1])
row.names(sc_data) <- ref_data[, 1] %>% unlist
sc_data <- ExpressionSet(assayData = sc_data)
bulk_data <- data.matrix(bulk_rna)
row.names(bulk_data) <- ref_data[, 1] %>% unlist
bulk_data <- ExpressionSet(assayData = bulk_data)

# Default
music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$Cells,
                        cell_size = NULL)

# Custom cell size
music_est_custom <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$Cells,
                        cell_size = cell_size)



# Results
music_est <- data.table(samples  = row.names(music_est$Est.prop.allgene), music_est$Est.prop.weighted)
music_est_custom <- data.table(samples  = row.names(music_est_custom$Est.prop.allgene), music_est_custom$Est.prop.weighted)