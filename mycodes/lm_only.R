
# Preamble ----------------------------------------------------------------

# Load libraries
packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer",  "mixtools","matrixStats")
loaded_libs <- sapply(packs, library, character.only = T)

rna_gene_x_cells <- fread(file.path(".", "cell_data", "rna_gene_x_cells.txt"))
cell_type_data <- fread(file.path(".", "cell_data", "cell_type_info.txt"))
cell_type_data <- cell_type_data[ all %in% c("Neurons", "OPC",  "Astrocytes",
                                             "Oligodendrocytes", "Microglia", "Endothelial"), ]



# Write functions ---------------------------------------------------------

# This script contain all the functions needed to estimate lambda, and mu from the reference data
# pi_rna, pi_cells, and mu require the genes in the rows and the cells in the columns
# the function, est_ref_params, calculates lambda*mu and returns it.
# The function, gene_datasets, generates test validation and train according to specified proportions
# Update 14th August 2018, changed the calculation of mu



gene_datasets <- function(test = NULL, train = NULL, validate = NULL, seed = NULL, cell_data = NULL){
  p <- c(test, train, validate)/sum(c(test, train, validate))
  test <- p[1];train <- p[2]; validate <- p[3]
  if (test + train + validate != 1)
    stop("The probabilites need to sum to one")
  
  if(!is.null(seed))
    set.seed(seed)
  
  
  my_list <- split(cell_data$cells, cell_data[,2])
  split_data <- lapply(my_list, function(x){
    splitSample <- sample(1:3, size=length(x), prob=c(validate,train,test), replace = TRUE)
    train.x <- x[splitSample==2]
    valid.x <- x[splitSample==1]
    test.x <- x[splitSample==3]
    list(train = train.x, valid = valid.x, test = test.x)
  })
  k <- length(split_data)
  k <- rep(1:3, each = k)
  
  type <- c("train", "valid", "test")
  
  # str(split_data)
  # str(do.call(rbind, split_data))
  
  split_data <- do.call(rbind, split_data)
  
  re <- lapply(1:3, function(i){
    sub <- which(k == i)
    re <-  unlist(split_data[sub])
    re <- unname(re)
    return(re)
  })
  names(re) <- type
  
  
  return(re)
  
}


est_ref_params <- function(cell_rna_data, cell_type_info){
  require(data.table)
  # QC
  # Drop genes with total sum == 0
  ctc <- cell_type_info$cells
  ctc <- intersect(colnames(cell_rna_data), ctc)
  cell_rna_data <- cell_rna_data[, c("genes", ctc), with = F]
  total_rna <- apply(cell_rna_data[,-1], 1, sum)
  keep_gene <- which(total_rna > 0)
  #subset to individuals in cell_rna_data
  names(cell_type_info)[1] <- "cells"
  cell_type_info <- cell_type_info[cells %in% ctc,]
  
  gene_x_cells <- cell_rna_data[keep_gene, -1]
  keep_gene <- cell_rna_data[keep_gene, genes]
  rlist <- split(cell_type_info$cells, cell_type_info[,2])
  
  
  # pi cells
  r <- nrow(gene_x_cells);m <- ncol(gene_x_cells)
  n_c <- sapply(rlist, length)
  pi_c <- n_c/m
  pi_c <- apply(t(pi_c), 2, rep, r)
  
  # pi rna
  total_rna <- apply(gene_x_cells, 1, sum)
  denom <- 1/total_rna
  
  
  numr <- sapply(rlist, function(x){
    c(apply(gene_x_cells[,x,with = F], 1, sum))
  })
  s_n <- colnames(numr)
  pi_r <- as.matrix(numr) * denom
  
  # mu rna
  # denom <- apply(t(n_c), 2, rep, r)
  # mu <- numr/denom
  
  # mu rna
  #update 14th August 2018 ---> recalulate numerator
  # numr <- apply(gene_x_cells, 1, sum)
  # numr <- apply(t(numr), 2, rep, length(s_n))
  # numr <- t(numr)
  # colnames(numr) <- s_n
  # denom <- apply(t(n_c), 2, rep, r)
  # mu <- numr/denom
  
  numr <- apply(gene_x_cells, 1, mean)
  numr <- apply(t(numr), 2, rep, length(s_n))
  numr <- t(numr)
  colnames(numr) <- s_n
  mu <- numr
  
  # Lambda
  l <- pi_r/pi_c
  out <- mu*l
  mu_c <- out
  
  out_old <- data.table(genes = keep_gene, data.matrix(out))
  lambda <- data.table(genes = keep_gene, data.matrix( pi_r/pi_c))
  
  #update August 15th 2018
  mean_rna <- apply(gene_x_cells, 1, mean)
  lm_data <- data.table(y = mean_rna, data.matrix(out))
  
  # LM
  n <- ncol(lm_data)
  
  # Obtain all columns except the last
  mysubset <- lapply(names(lm_data)[1:(n-1)], function(x) {c(x, names(lm_data)[n])})
  lm_x <- lm_data[, lapply(mysubset, function(x) (get(x[1])) - (get(x[2])) )]
  fit_x <- lm(V1 ~ -1 + ., data = lm_x)
  beta <- coef(fit_x)
  beta <- c(beta)
  
  beta <- c(beta, 1 - sum(beta))
  names(beta) <- s_n
  l <- beta/unique(pi_c)
  # l <- 1/l
  l <- apply(l, 2, rep, r)
  
  out <- mu_c*l
  out <- data.table(genes = keep_gene, data.matrix(out))
  
  mu <- data.table(genes = keep_gene, data.matrix(mu))
  mean_rna <- data.table(genes = keep_gene, data.matrix(mean_rna))
  
  # Variances
  v_e <- sapply(rlist, function(x){
    e <- c(apply(gene_x_cells[,x,with = F], 1, mean))
    e
  })
  v_e <- apply(v_e, 1, var)
  
  e_v <- sapply(rlist, function(x){
    v <- c(apply(gene_x_cells[,x,with = F], 1, var))
    v
  })
  e_v <- apply(e_v, 1, mean)
  
  t_v <- cbind(v_e,v_e + e_v)
  
  p_v <- apply(t_v, 1, function(y) { ifelse(y[2] > 0, y[1]/y[2], 0) })
  
  p_v <- data.table(genes = keep_gene, prop.v_e = p_v)
  
  
  out <- list(out = out, pi_c = unique(pi_c), pi_r = data.table(pi_r), p_v = p_v, mu_c = mu, out_old = out_old, mean_rna = mean_rna, lambda = lambda)
  return(out)
}



est_pi_c <- function(my_data_all = rna_gene_x_cells, cell_type_info = cell_type_data, ppv = .1, lower = T, test,...){
  
  # Generate and use training data------
  
  data_split <- gene_datasets(cell_data = cell_type_info, test = test, ...)
  my_data_train <- my_data_all[, c(names(my_data_all)[1], data_split$train), with = F]
  params <- est_ref_params(my_data_train, cell_type_info)
  pv <- params$p_v
  my_vals <- params$out
  
  # Subset to informative genes---------
  if(lower){
    my_vals <- my_vals[genes %in% pv[prop.v_e <= ppv, genes], ]
  }else(
    my_vals <- my_vals[genes %in% pv[prop.v_e > ppv, genes], ]
  )
  
  # Load libs------
  packs <- c("data.table", "dplyr")
  libs_loaded <- sapply(packs, require, character.only = T)
  
  # Prep validation data-----
  my_val_data <- my_data_all[, c(names(my_data_all)[1], data_split$valid), with = F]
  my_val_data <- my_val_data[genes %in% my_vals$genes]
  
  
  # get data ---------
  y <- apply(my_val_data[,-1], 1, mean, na.rm = T)
  
  # LM---------
  n <- ncol(my_vals)
  
  # Obtain all columns except the last
  mysubset <- lapply(names(my_vals)[2:(n-1)], function(x) {c(x, names(my_vals)[n])})
  
  x <- my_vals[, lapply(mysubset, function(x) get(x[1]) - get(x[2])  )]
  y_star <- (y - unlist(my_vals[,n, with = F]))
  
  x <- scale(x, scale = F)
  x <- cbind(1, x)
  y_star <- scale(y_star, scale = F)
  
  if(!all(is.na(y)) & !all(is.na(x))){
    ft <- lm(y_star~ x[,-1])
    beta <- coef(ft)
    beta <- c(beta)
    beta <- c(beta[-1], 1 - sum(beta[-1]))
  }else{
    beta <- rep(NA, (n-1))
  }
  
  names(beta) <- names(my_vals)[2:(n)]
  
  
  # LM2---------
  
  params_out_m <- params$mean_rna
  params_out_l <- params$lambda
  # Subset to informative genes
  if(lower){
    params_out_m <- params_out_m[genes %in% pv[prop.v_e <= ppv, genes], ]
    params_out_l <- params_out_l[genes %in% pv[prop.v_e <= ppv, genes], ]
  }else{
    params_out_m <- params_out_m[genes %in% pv[prop.v_e > ppv, genes], ]
    params_out_l <- params_out_l[genes %in% pv[prop.v_e > ppv, genes], ]
  }
  
  y_n <- y/params_out_m$V1
  my_dt <- data.table(y = y_n, params_out_l[, -1])
  my_dt <- data.frame(my_dt)
  my_dt[is.na(my_dt)] <- NA
  my_dt <- data.table(na.omit(my_dt))
  
  # y_n <- y
  # my_dt <- data.table(y = y_n, params_out_l[, -1]*y_n)
  
  n <- ncol(my_dt)
  # Obtain all columns except the last
  mysubset <- lapply(names(my_dt)[1:(n-1)], function(x) {c(x, names(my_dt)[n])})
  my_dt <- my_dt[, lapply(mysubset, function(x) scale(get(x[1]), scale = F) - scale(get(x[2]), scale = F)  )]
  
  test_my_dt <- apply(my_dt, 2, function(x){
    all(is.na(x))
  })
  test_my_dt <- sum(test_my_dt)
  test_my_dt <- test_my_dt == ncol(my_dt)
  if(!test_my_dt){
    ft <- lm(V1 ~ -1 + ., data = my_dt)
    b <- coef(ft)
    b <- c(b, 1-sum(b))
  }else{
    b <- rep(NA, (n-1))
  }
  
  names(b) <- names(params_out_l)[-1]
  
  
  
  # Truth ---------
  # all_pi_c <- unique(est_ref_params(my_data_all, cell_type_info)$pi_c)
  all_pi_c <- unique(est_ref_params(my_data_all[, c(names(my_data_all)[1], data_split$valid), with = F], cell_type_info)$pi_c) 
  
  # Output 
  a <- names(my_vals)[-1]
  beta <- beta[match(a,names(beta))]
  b <- b[match(a,names(b))]
  
  re <- list(LM = beta, LM2 = b,truth = all_pi_c)
  re <- do.call(rbind, re)
  re <- data.table(method = c("LM", "LM2", "Truth"), re)
  
  vt_data <- apply(my_val_data[, -1], 1, mean)
  vt_data <- data.table(genes = my_val_data$genes, ave = vt_data)
  out <- list(re = re, val_dt = vt_data)
  return(out)
  
  
}




# All cell types testing ----------------------------------------------------------
# Neurons----
cell_type_data <- fread(file.path(".", "cell_data", "cell_type_info.txt"))
cell_type_data <- cell_type_data[ all %in% c("Neurons", "OPC",  "Astrocytes",
                                             "Oligodendrocytes", "Microglia", "Endothelial"), ]

p_t <- .0
p_tr <- .5
p_tv <- .5
cell_type_data[, Neurons:= ifelse(Neurons, "N+", "N-")]
rna_gene_sub <- rna_gene_x_cells
rna_gene_sub <- rna_gene_sub[, c("genes", cell_type_data$cells), with = F]
# neu <- est_pi_c(my_data_all = rna_gene_sub, cell_type_info = cell_type_data, test = p_t, train = p_tr, validate = p_tv, seed= 368, ppv = 0.05)
# neu
# 
# # use upper
# neu2 <- est_pi_c(my_data_all = rna_gene_sub, cell_type_info = cell_type_data, test = p_t, train = p_tr, validate = p_tv, seed= 368, ppv = 0.5, lower = F)
# neu2
# 
# 
# neu_all_genes <- est_pi_c(my_data_all = rna_gene_sub, cell_type_info = cell_type_data, test = p_t, train = p_tr, validate = p_tv, seed= 368, ppv = -1, lower = F)
# neu_all_genes

# sims----
# show error over all 1e3 iterations using all genes
packs <- c("parallel")
loaded_libs <- sapply(packs, library, character.only = T)
ncpu <- detectCores() - 1
cl <- makeCluster(ncpu)
clusterExport(cl, c("rna_gene_sub", "cell_type_data", "p_t", "p_tr", "p_tv",
                    "gene_datasets","est_ref_params", "est_pi_c"))
sims1 <- parSapply(cl, 1:1e3, function(i){
  require("data.table")
  re <- est_pi_c(my_data_all = rna_gene_sub, cell_type_info = cell_type_data, 
                            test = p_t, train = p_tr, validate = p_tv, ppv = -1, lower = F)
  
  neu_all_genes <- re$re
  sse <- lapply(seq_len(nrow(neu_all_genes) - 1), function(x){
    neu <- neu_all_genes
    a <- unlist(neu[x, -1])
    b <- unlist(neu[nrow(neu), -1])
    d <- a - b
    c(sum(d*d))
    
  })
  names(sse) <- c("LM", "LM2")
  return(sse)
})

sims1 <- t(sims1)


# Check the effect of ppv
sims2 <- parSapply(cl, seq(0, .7, by = .1), function(p){
  require("data.table")
  out <- sapply(1:1e3, function(i){
    re <- est_pi_c(my_data_all = rna_gene_sub, cell_type_info = cell_type_data, test = p_t, train = p_tr, validate = p_tv, ppv = p, lower = F)
    neu_all_genes <- re$re
    sse <- sapply(seq_len(nrow(neu_all_genes) - 1), function(x){
      neu <- neu_all_genes
      a <- unlist(neu[x, -1])
      b <- unlist(neu[nrow(neu), -1])
      d <- a - b
      c(sum(d*d))
      
    })
    names(sse) <- c("LM", "LM2")
    return(sse)
  })
  out <- t(out)
  out <- apply(out, 2, mean)
  return(out)
})
sims2 <- t(sims2)
sims2 <- cbind(seq(0, .7, by = .1), sims2)

fwrite(data.table(sims1), file.path(".","sim_results", "neurons_all.txt"))
fwrite(data.table(sims2), file.path(".","sim_results", "neurons_ppv.txt"))

# ALL cell types----
cell_type_data[, Neurons:= NULL]
rna_gene_sub <- rna_gene_x_cells
rna_gene_sub <- rna_gene_sub[, c("genes", cell_type_data$cells), with = F]

# all <- est_pi_c(my_data_all = rna_gene_sub, cell_type_info = cell_type_data, test = p_t, train = p_tr, validate = p_tv, seed= 3434, ppv = 1e-1)
# all
# 
# all_upper <- est_pi_c(my_data_all = rna_gene_sub, cell_type_info = cell_type_data, test = p_t, train = p_tr, validate = p_tv, seed= 3434, ppv = 5e-1, lower = F)
# all_upper
# 
# all_genes <- est_pi_c(my_data_all = rna_gene_sub, cell_type_info = cell_type_data, test = p_t, train = p_tr, validate = p_tv, seed= 3434, ppv = -1, lower = F)
# all_genes

clusterExport(cl, c("rna_gene_sub", "cell_type_data", "p_t", "p_tr", "p_tv",
                    "gene_datasets","est_ref_params", "est_pi_c"))
sims1 <- parSapply(cl, 1:1e3, function(i){
  require("data.table")
  re <- est_pi_c(my_data_all = rna_gene_sub, cell_type_info = cell_type_data, test = p_t, train = p_tr, validate = p_tv, ppv = -1, lower = F)
  neu_all_genes <- re$re
  sse <- sapply(seq_len(nrow(neu_all_genes) - 1), function(x){
    neu <- neu_all_genes
    a <- unlist(neu[x, -1])
    b <- unlist(neu[nrow(neu), -1])
    d <- a - b
    c(sum(d*d))
    
  })
  names(sse) <- c("LM", "LM2")
  return(sse)
})

sims1 <- t(sims1)

# Check the effect of ppv
sims2 <- parSapply(cl, seq(0, .7, by = .1), function(p){
  require("data.table")
  out <- sapply(1:5e2, function(i){
    re <- est_pi_c(my_data_all = rna_gene_sub, cell_type_info = cell_type_data, test = p_t, train = p_tr, validate = p_tv, ppv = p, lower = F)
    neu_all_genes <- re$re
    sse <- sapply(seq_len(nrow(neu_all_genes) - 1), function(x){
      neu <- neu_all_genes
      a <- unlist(neu[x, -1])
      b <- unlist(neu[nrow(neu), -1])
      d <- a - b
      c(sum(d*d))
      
    })
    names(sse) <- c("LM", "LM2")
    return(sse)
  })
  out <- t(out)
  out <- apply(out, 2, mean, na.rm = T)
  return(out)
})
sims2 <- t(sims2)

fwrite(data.table(sims1), file.path(".","sim_results", "all_types.txt"))
fwrite(data.table(sims2), file.path(".","sim_results", "all_types_ppv.txt"))


stopCluster(cl)