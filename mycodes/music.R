
# Preamble ----------------------------------------------------------------

packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer", 
           "mixtools","matrixStats", "MuSiC", "convert", "xbioc", "ggplot2", "sva")
libs_loaded <- sapply(packs, library, character.only = T)

# MuSic -------------------------------------------------------------------

rna_gene_x_cells <- fread(file.path(".", "cell_data", "rna_gene_x_cells.txt"))
cell_type_data <- fread(file.path(".", "cell_data", "cell_type_info.txt"))
cell_type_data <- cell_type_data[ all %in% c("Neurons", "OPC",  "Astrocytes",
                                             "Oligodendrocytes", "Microglia", "Endothelial"), ]


# Subset to informative genes
my_vals <- est_ref_params(rna_gene_x_cells, cell_type_data)
my_params <- my_vals$out
# my_params <- my_params[genes %in% pv[prop.v_e > 0.1, genes], ]
params_out_m <- my_vals$mean_rna
# params_out_m <- params_out_m[genes %in% pv[prop.v_e > 0.2, genes], ]
params_out_l <- my_vals$lambda
# params_out_l <- params_out_l[genes %in% pv[prop.v_e > 0.2, genes], ]

load("./NAc_Nicotine_hg38_rseGene_rawCounts_postQCSamples_n223.rda")
load("./MatchedCellComp/singleCell_iPSC_quake_coefEsts_calibration_Zscale_adultOnly.rda")
yExprs <- log2(getRPKM(rse_gene, "Length")+1)


# project
yExprs <- yExprs[rownames(coefEsts),]
yExprs_Z <- scale(yExprs[rownames(coefEsts),])


my_params <- my_params[genes %in% rownames(coefEsts),]

gg <- rownames(coefEsts)
gg2 <- gg[gg %in% params_out_m$genes]
gg <- gg[gg %in% my_params$genes]

my_params <- my_params[match(gg, genes),]
params_out_m <- params_out_m[match(gg2, genes),]
params_out_l <- params_out_l[match(gg2, genes),]


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

data.table::fwrite(data.table::data.table(results = my_s_k, cell_types = types_cell), "music_cell_size.txt")

# plot_test ---------------------------------------------------------------

{
  true_params <- est_ref_params(rna_gene_x_cells, cell_type_data)
  
  # Raw
  cell_means <- true_params$out
  cell_means <- inner_join(data.table(genes = row.names(yExprs)), cell_means) %>% data.table
  prop_info <- true_params$p_v
  names(prop_info)[2] <- "prop"
  prop_info <- inner_join(data.table(genes = row.names(yExprs)), prop_info) %>% data.table
  
  sc_1 <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 1e4, scaling = T, prop_data = prop_info, 
                          all_dt = cbind(0, yExprs), raw_dt = cell_means, eval = F, tol = 1e-8)
  
  ot_1 <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 0, scaling = T, prop_data = prop_info, 
                          all_dt = cbind(0, yExprs), raw_dt = cell_means, eval = F, tol = 1e-8)
  
  sc_1_np <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 1e4, scaling = T, prop_data = prop_info, 
                             all_dt = cbind(0, yExprs), raw_dt = cell_means, eval = F, tol = 1e-8, use.prior = F)
  
  # p > .2
  gg <- prop_info$genes[prop_info$prop > .2]
  prop_info2 <- prop_info[genes %in% gg,]
  cell_means_gg <- cell_means[genes %in% gg,]
  analysis_dt2 <- which(row.names(yExprs) %in% gg)
  analysis_dt2 <- yExprs[analysis_dt2, ]
  
  identical(row.names(analysis_dt2), cell_means_gg$genes)
  sc_1_02 <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 1e4, scaling = T, prop_data = prop_info2, 
                          all_dt = cbind(0, analysis_dt2), raw_dt = cell_means_gg, eval = F, tol = 1e-8)
  ot_1_02 <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 0, scaling = T, prop_data = prop_info2, 
                             all_dt = cbind(0, analysis_dt2), raw_dt = cell_means_gg, eval = F, tol = 1e-8)
  
  fwrite(data.table(genes = row.names(analysis_dt2)), "informative_genes.txt")
  
  
  
  # Model Averaged
  # With scaling
  iter <- 1e3
  n_ave <- 250
  B <- seq(0, .7, length.out = n_ave)
  
  # B <- quantile(true_params$p_v$prop.v_e, probs = B, type = 2)
  
  sample_size <- ncol(yExprs)
  
  sc_1_all <- lapply(B, my_sim_p_abbas3, use.sva = T, iter = iter, scaling = T, eval = F, prop_data = prop_info, all_dt = cbind(0, yExprs), raw_dt = cell_means)
  sc_1_all <- lapply(sc_1_all, function(df) {df[, -c(2), with = F]} )
  sc_1_all <- do.call(rbind, sc_1_all)
  sc_1_all[, p:=rep(B, each = sample_size) + 0.0050251]
  sc_1_all[, `:=`(F_1 = `FALSE`*p, `T_1` = `TRUE`*p)]
  sc_1_all_eval <- sc_1_all[, .(F_1 = sum(F_1), T_1 = sum(T_1), all_p = sum(p)), by = ID]
  sc_1_all_eval[, `:=`( F_1 = F_1/all_p, T_1 = T_1/all_p)]
  sc_1_all_eval[, c("all_p") := NULL]
  
  sc_1[, RNum := sapply(strsplit(ID, "_"), "[", 1)]
  ot_1[, RNum := sapply(strsplit(ID, "_"), "[", 1)]
  sc_1_np[, RNum := sapply(strsplit(ID, "_"), "[", 1)]
  sc_1_all_eval[, RNum := sapply(strsplit(ID, "_"), "[", 1)]
  sc_1_02[, RNum := sapply(strsplit(ID, "_"), "[", 1)]
  ot_1_02[, RNum := sapply(strsplit(ID, "_"), "[", 1)]
  load("./NAc_rse_gene_withCompEsts.rda")
  sc_1 <- sc_1[match(rse_gene$RNum, RNum),]
  ot_1 <- ot_1[match(rse_gene$RNum, RNum),]
  sc_1_np <- sc_1_np[match(rse_gene$RNum, RNum),]
  sc_1_all_eval <- sc_1_all_eval[match(rse_gene$RNum, RNum),]
  sc_1_02 <- sc_1_02[match(rse_gene$RNum, RNum),]
  ot_1_02 <- ot_1_02[match(rse_gene$RNum, RNum),]
}

palette(brewer.pal(4,"Dark2"))


r_2_size <- 1.2



r_2 <- cor(rse_gene$NeuN_neg_DNAm , sc_1_np$`FALSE`)^2
r_2 <- round(r_2, 2)
plot(rse_gene$NeuN_neg_DNAm , sc_1_np$`FALSE`,
     ylab = "LM", xlab = "DNAm-based", 
     main = "Cell Composition Estimation\n(Houseman RNA VS Houseman Methylation)",
     pch = 21, bg = 1, ylim = c(0,1),xlim=c(0,1))
abline(0,1,lty=2,lwd=2)
points(rse_gene$NeuN_pos_DNAm,  sc_1_np$`TRUE`,
       pch = 21, bg = 2)
legend("topleft", c("NeuN-", "NeuN+"), pch = 15, col=1:2, horiz = T, cex = .8, text.width = .15)
text(.9, .15, substitute(paste(R^2, nn), list(nn=paste("=", r_2))) ,
     cex = r_2_size)

png("./figs/comparisons.png", width=10, height=10, units="in", res=300)
par(mfrow = c(2,3), pty = "m")
r_2 <- cor(rse_gene$NeuN_neg_DNAm , music_est$`FALSE`)^2
r_2 <- round(r_2, 2)
plot(rse_gene$NeuN_neg_DNAm , music_est$`FALSE`,
     ylab = "MuSiC", xlab = "DNAm-based", 
     main = "Cell Composition Estimation\n(MuSic RNA VS Houseman Methylation)",
     pch = 21, bg = 1, ylim = c(0,1),xlim=c(0,1))
abline(0,1,lty=2,lwd=2)
points(rse_gene$NeuN_pos_DNAm,  music_est$`TRUE`,
       pch = 21, bg = 2)
legend("topleft", c("NeuN-", "NeuN+"), pch = 15, col=1:2, horiz = T, cex = .8, text.width = .15)
text(.9, .15, substitute(paste(R^2, nn), list(nn=paste("=", r_2))) ,
     cex = r_2_size)


r_2 <- cor(rse_gene$NeuN_neg_DNAm , sc_1$`FALSE`)^2
r_2 <- round(r_2, 2)
plot(rse_gene$NeuN_neg_DNAm , sc_1$`FALSE`,
     ylab = "2-Step", xlab = "DNAm-based", 
     main = "Cell Composition Estimation\n(2-Step;all genes VS Houseman Methylation)",
     pch = 21, bg = 1, ylim = c(0,1),xlim=c(0,1))
abline(0,1,lty=2,lwd=2)
points(rse_gene$NeuN_pos_DNAm,  sc_1$`TRUE`,
       pch = 21, bg = 2)
legend("topleft", c("NeuN-", "NeuN+"), pch = 15, col=1:2, horiz = T, cex = .8, text.width = .15)
text(.9, .15, substitute(paste(R^2, nn), list(nn=paste("=", r_2))) ,
     cex = r_2_size)

r_2 <- cor(rse_gene$NeuN_neg_DNAm , ot_1$`FALSE`)^2
r_2 <- round(r_2, 2)
plot(rse_gene$NeuN_neg_DNAm , ot_1$`FALSE`,
     ylab = "LM", xlab = "DNAm-based", 
     main = "Cell Composition Estimation\n(LM;all genes VS Houseman Methylation)",
     pch = 21, bg = 1, ylim = c(0,1),xlim=c(0,1))
abline(0,1,lty=2,lwd=2)
points(rse_gene$NeuN_pos_DNAm,  ot_1$`TRUE`,
       pch = 21, bg = 2)
legend("topleft", c("NeuN-", "NeuN+"), pch = 15, col=1:2, horiz = T, cex = .8, text.width = .15)
text(.9, .15, substitute(paste(R^2, nn), list(nn=paste("=", r_2))) ,
     cex = r_2_size)


r_2 <- cor(rse_gene$NeuN_neg_DNAm , sc_1_all_eval$F_1)^2
r_2 <- round(r_2, 2)
plot(rse_gene$NeuN_neg_DNAm , sc_1_all_eval$F_1,
     ylab = "Model averaged 2-Step", xlab = "DNAm-based", 
     main = "Cell Composition Estimation\n(2-Step VS Houseman Methylation)",
     pch = 21, bg = 1, ylim = c(0,1),xlim=c(0,1))
abline(0,1,lty=2,lwd=2)
points(rse_gene$NeuN_pos_DNAm,  sc_1_all_eval$T_1,
       pch = 21, bg = 2)
legend("topleft", c("NeuN-", "NeuN+"), pch = 15, col=1:2, horiz = T, cex = .8, text.width = .15)
text(.9, .15, substitute(paste(R^2, nn), list(nn=paste("=", r_2))) ,
     cex = r_2_size)

r_2 <- cor(rse_gene$NeuN_neg_DNAm , sc_1_02$`FALSE`)^2
r_2 <- round(r_2, 2)
plot(rse_gene$NeuN_neg_DNAm , sc_1_02$`FALSE`,
     ylab = "2-Step", xlab = "DNAm-based", 
     main = "Cell Composition Estimation\n(2-Step;p>0.2 VS Houseman Methylation)",
     pch = 21, bg = 1, ylim = c(0,1),xlim=c(0,1))
abline(0,1,lty=2,lwd=2)
points(rse_gene$NeuN_pos_DNAm,  sc_1_02$`TRUE`,
       pch = 21, bg = 2)
legend("topleft", c("NeuN-", "NeuN+"), pch = 15, col=1:2, horiz = T, cex = .8, text.width = .15)
text(.9, .15, substitute(paste(R^2, nn), list(nn=paste("=", r_2))) ,
     cex = r_2_size)

r_2 <- cor(rse_gene$NeuN_neg_DNAm , ot_1_02$`FALSE`)^2
r_2 <- round(r_2, 2)
plot(rse_gene$NeuN_neg_DNAm , ot_1_02$`FALSE`,
     ylab = "LM", xlab = "DNAm-based", 
     main = "Cell Composition Estimation\n(LM;p>0.2 VS Houseman Methylation)",
     pch = 21, bg = 1, ylim = c(0,1),xlim=c(0,1))
abline(0,1,lty=2,lwd=2)
points(rse_gene$NeuN_pos_DNAm,  ot_1_02$`TRUE`,
       pch = 21, bg = 2)
legend("topleft", c("NeuN-", "NeuN+"), pch = 15, col=1:2, horiz = T, cex = .8, text.width = .15)
text(.9, .15, substitute(paste(R^2, nn), list(nn=paste("=", r_2))) ,
     cex = r_2_size)
dev.off()


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
music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = rep(c( "Jurkat", "IM-9", "Raji", "THP-1"), each = 3), samples = cn, normalize = T)
music_est_nnls <- music_est$Est.prop.allgene
music_est <- music_est$Est.prop.weighted

eval_mix(music_est_nnls)
eval_mix(music_est)

music_est_nnls <- data.table("ID" = sapply(strsplit(row.names(music_est_nnls), " "), "[", 1),
                 "Mix" = sapply(strsplit(row.names(music_est_nnls), " "), "[", 2),
                 music_est_nnls)
music_est <- data.table("ID" = sapply(strsplit(row.names(music_est), " "), "[", 1),
                 "Mix" = sapply(strsplit(row.names(music_est), " "), "[", 2),
                 music_est)
# keep_idx <- which(p.adjust(cell_type_data_var$p_vals, "BH") <= 5e-2)
# 
# sample_ests_ts_info <- my_sim_p_abbas3(p = 0, use.sva = T, iter = 1e4, scaling = T, prop_data = cell_type_data_var[keep_idx,],
#                 all_dt = mixes_abbas[keep_idx,], raw_dt = cell_type_data[keep_idx,], eval = F)

# truth
MixA <- c("Jurkat" = 2.5, "IM-9" = 1.25, "Raji" =  2.5,"THP-1" = 3.75)/10
MixB <-  c("Jurkat" = 0.5, "IM-9" = 3.17, "Raji" = 4.75, "THP-1" = 1.58)/10
MixC <-  c("Jurkat" = 0.1, "IM-9" = 4.95, "Raji" = 1.65, "THP-1" = 3.3)/10
MixD  <- c("Jurkat" = 0.02, "IM-9" = 3.33, "Raji" = 3.33, "THP-1" = 3.33)/10

Mix <- rbind(MixA, MixB, MixC, MixD)
Mix <- data.table("Mix" = paste0("Mix", LETTERS[1:4]), Mix)

plot_dt_out <- rbind(inner_join(music_est[, .(Mix, Estimate = Jurkat)],  Mix[, .(Mix, Truth = Jurkat, type = "Jurkat")]),
                    inner_join(music_est[, .(Mix, Estimate = `IM-9`)],  Mix[, .(Mix, Truth = `IM-9`, type = "IM-9")]),
                    inner_join(music_est[, .(Mix, Estimate = Raji)],  Mix[, .(Mix, Truth = Raji, type = "Raji")]),
                    inner_join(music_est[, .(Mix, Estimate = `THP-1`)],  Mix[, .(Mix, Truth = `THP-1`, type = "THP-1")])) %>% data.table
with(plot_dt_out, cor(x = Truth, y = Estimate))
fwrite(plot_dt_out, file.path(".", "sim_results", "plot_dt_music_abbas_v2.txt"))


# Misspecification of the reference dataset -------------------------------

loop <- as.list(1:4)
g <- t(combn(1:4, 2))
g <- split(g, 1:nrow(g))
g <- c(loop, g)
rn <- raw_dt$`!Sample_title`
bulk_data <- ExpressionSet(assayData = music_dt)
music_dt <- data.matrix(all_dt[, -1])
row.names(music_dt) <- rn
all <- seq(1,12, by = 3)
sapply(g, function(k){
  
  k_old <- k
  if(length(k) == 1){
    k <- all[k]:(all[k]+2)
  }else{
    k <- c(all[k[1]]:(all[k[1]]+2), all[k[2]]:(all[k[2]]+2))
  }
  
  cov <- rep(c( "Jurkat", "IM-9", "Raji", "THP-1"), each = 3)
  cov <- cov[-k]
  

  print(cov)
  
  # Analyze
  sc_data <- raw_dt
  rn <- sc_data$`!Sample_title`
  cn <- colnames(sc_data)[-1]
  cn <- cn[-k]
  sc_data <- data.matrix(sc_data[, -1])[,-k]
  row.names(sc_data) <- rn
  sc_data <- ExpressionSet(assayData = sc_data)
  music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cov, samples = cn, normalize = T)
  music_est_nnls <- music_est$Est.prop.allgene
  music_est <- music_est$Est.prop.weighted
  
 
  
  music_est <- data.table("ID" = sapply(strsplit(row.names(music_est), " "), "[", 1),
                          "Mix" = sapply(strsplit(row.names(music_est), " "), "[", 2),
                          music_est)
  
  # truth
  MixA <- c("Jurkat" = 2.5, "IM-9" = 1.25, "Raji" =  2.5,"THP-1" = 3.75)/10
  MixB <-  c("Jurkat" = 0.5, "IM-9" = 3.17, "Raji" = 4.75, "THP-1" = 1.58)/10
  MixC <-  c("Jurkat" = 0.1, "IM-9" = 4.95, "Raji" = 1.65, "THP-1" = 3.3)/10
  MixD  <- c("Jurkat" = 0.02, "IM-9" = 3.33, "Raji" = 3.33, "THP-1" = 3.33)/10
  
  k <- k_old
  Mix <- rbind(MixA, MixB, MixC, MixD)
  Mix <- Mix[, -(k)]
  Mix <- data.table("Mix" = paste0("Mix", LETTERS[1:4]), Mix)
  
  r <- names(Mix)[-1]
  plot_dt_ts <- lapply(r, function(x){
    dt <- music_est[, c("Mix", x), with = F]
    dt_m <- Mix[,  c("Mix", x), with = F]
    dt_m$type <- x
    names(dt) <- c("Mix", "Estimate")
    names(dt_m) <- c("Mix", "Truth", "type")
    out <- inner_join(dt,  dt_m) %>% data.table
    return(out)
  })
  plot_dt_ts <- do.call(rbind, plot_dt_ts)
  with(plot_dt_ts, round(cor(x = Truth, y = Estimate), 2))
})

sapply(g, function(k){
  # truth
  MixA <- c("Jurkat" = 2.5, "IM-9" = 1.25, "Raji" =  2.5,"THP-1" = 3.75)/10
  MixB <-  c("Jurkat" = 0.5, "IM-9" = 3.17, "Raji" = 4.75, "THP-1" = 1.58)/10
  MixC <-  c("Jurkat" = 0.1, "IM-9" = 4.95, "Raji" = 1.65, "THP-1" = 3.3)/10
  MixD  <- c("Jurkat" = 0.02, "IM-9" = 3.33, "Raji" = 3.33, "THP-1" = 3.33)/10
  
  Mix <- rbind(MixA, MixB, MixC, MixD)
  colnames(Mix)[(k - 1)]
}) %>% c

# Mixed/ Unknown content -------------------------------------------------------------------

all_objs <- ls()
all_objs <- all_objs[!all_objs %in% c("mixes_abbas", "cell_type_data_var", "cell_type_data",
                                      "eval_mix", "my_sim_p_abbas3")]
rm(list = all_objs)
mixes_abbas_hct <- fread(file.path(".", "cell_data", "abbas_et_al_data", "FinalMixtureMatrix.Abbas.HCT116.Res30.subset.txt"))
# ciber_sig <- fread(file.path(".", "cell_data", "abbas_et_al_data", "GSE11103_matrix_classes.GSE11103_matrix_pure.bm.K999.0.txt"))
# ciber_sig[, V6:=NULL]
keep_idx <- which(mixes_abbas$`!Sample_title` %in% mixes_abbas_hct$Probeset)

ciber_sig <- cell_type_data[keep_idx,]


names(ciber_sig)[1] <- names(cell_type_data)[1]
all <- seq(1, 64, by = 4) + 1

# Normalize
{
  yExprs_normfactors <- edgeR::calcNormFactors(cbind(data.matrix(mixes_abbas_hct[, -1]), ciber_sig[, -c(1), with = F]), method = "TMM" )
  a <- ncol(mixes_abbas_hct)  - 1 + 1
  b <- a + ncol(ciber_sig[, -c(1), with = F]) - 1
  q_out <- yExprs_normfactors[a:b]
  q_out <- replicate(nrow(mixes_abbas_hct)  , q_out)
  q_out <- t(q_out)
  raw_dt <- cbind(ciber_sig[, 1], ciber_sig[, -c(1), with = F]*q_out)
  
  yExprs_normfactors <- yExprs_normfactors[-c(a:b)]
  yExprs_normfactors <- replicate(nrow(mixes_abbas_hct), yExprs_normfactors)
  yExprs_normfactors <- t(yExprs_normfactors)
  all_dt <- cbind(mixes_abbas_hct[, 1], mixes_abbas_hct[,-1]*yExprs_normfactors)
}


# Analyze


sc_data <- raw_dt
rn <- sc_data$`!Sample_title`
cn <- colnames(sc_data)[-1]
sc_data <- data.matrix(sc_data[, -1])
row.names(sc_data) <- rn
music_dt <- data.matrix(all_dt[, -1])
row.names(music_dt) <- rn
sc_data <- ExpressionSet(assayData = sc_data)
bulk_data <- ExpressionSet(assayData = music_dt)
music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = rep(c( "Jurkat", "IM-9", "Raji", "THP-1"), each = 3), samples = cn, normalize = T)
music_est_nnls <- music_est$Est.prop.allgene
music_est <- music_est$Est.prop.weighted

eval_mix(music_est_nnls)
eval_mix(music_est)

music_est_nnls <- data.table("ID" = sapply(strsplit(row.names(music_est_nnls), " "), "[", 1),
                             "Mix" = sapply(strsplit(row.names(music_est_nnls), " "), "[", 2),
                             music_est_nnls)
music_est <- data.table("ID" = sapply(strsplit(row.names(music_est), " "), "[", 1),
                        "Mix" = sapply(strsplit(row.names(music_est), " "), "[", 2),
                        music_est)







# out <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 1e4, tol = 1e-10, eval = F, use.sum = NULL, scaling = T, 
#                        prop_data = cell_type_data_var[keep_idx, ], all_dt = analysis_dt, raw_dt = cell_type_data[keep_idx,])
out <- music_est
t_info <- sapply(strsplit(colnames(mixes_abbas_hct)[-1], "\\("), "[", 2)
t_info <- sapply(strsplit(t_info, "\\)"), "[", 1)
s_info <-  sapply(strsplit(t_info, ","), "[", 2)
t_info <-  sapply(strsplit(t_info, ","), "[", 1)
t_info <- trimws(t_info);s_info <- trimws(s_info)

out[, `:=`(noise_info = sapply(strsplit(s_info, " "), "[", 1),
           tumor_info = sapply(strsplit(t_info, " "), "[", 1),
           Mix = paste0(ID, Mix))]

# truth
MixA <- c("Jurkat" = 2.5, "IM-9" = 1.25, "Raji" =  2.5,"THP-1" = 3.75)/10
MixB <-  c("Jurkat" = 0.5, "IM-9" = 3.17, "Raji" = 4.75, "THP-1" = 1.58)/10
MixC <-  c("Jurkat" = 0.1, "IM-9" = 4.95, "Raji" = 1.65, "THP-1" = 3.3)/10
MixD  <- c("Jurkat" = 0.02, "IM-9" = 3.33, "Raji" = 3.33, "THP-1" = 3.33)/10



Mix <- rbind(MixA, MixB, MixC, MixD)
Mix <- data.table("Mix" = paste0("Mix", LETTERS[1:4]), Mix)


# Evaluate

cor_out <- apply(out, 1, function(x){
  dt <- as.numeric(x[3:6])
  m <- x[2]
  t <- Mix[Mix == m, -1, with = F]
  t <- unlist(t)
  if(m != "MixE"){
    r <- cor(dt, t)
    d <- dt-t
    d <- crossprod(d)/length(d)
    out <- c("pearson" = r, "rmse" = d) 
  }else{
    out <- c("pearson" = NA, "rmse" = NA) 
  }
  
  return(out)
})
cor_out <- t(cor_out)
out <- cbind(out, cor_out)
fwrite(out, file.path(".", "sim_results", "mixes_out_music.txt"))
