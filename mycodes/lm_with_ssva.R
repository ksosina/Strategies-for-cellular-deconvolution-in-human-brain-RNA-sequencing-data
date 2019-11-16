
# Preamble ----------------------------------------------------------------

# Load libraries
packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer",  
           "mixtools","matrixStats", "ggplot2", "sva", "ggplot2")
loaded_libs <- sapply(packs, library, character.only = T, quietly = T)


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



est_pi_c <- function(my_data_all = rna_gene_x_cells, yes_plot = F, scaling = F, cell_type_info = cell_type_data, ppv = .1, lower = T, test, bulk = F,...){
  
  
  {
    
    # Generate and use training data
    data_split <- gene_datasets(cell_data = cell_type_info, test = test, ...)
    my_data_train <- my_data_all[, c(names(my_data_all)[1], data_split$train), with = F]
    params <- est_ref_params(my_data_train, cell_type_info)
    pv <- params$p_v
    my_vals <- params$out
    
    
    
    # Subset to informative genes
    my_vals_em <-  my_vals
    if(lower){
      my_vals <- my_vals[genes %in% pv[prop.v_e < ppv, genes], ]
    }else(
      my_vals <- my_vals[genes %in% pv[prop.v_e > ppv, genes], ]
    )
    
    # Load libs
    packs <- c("data.table", "dplyr", "mixtools", "SummarizedExperiment", "penalized")
    libs_loaded <- sapply(packs, require, character.only = T)
    
    # Prep validation data
    my_val_data <- my_data_all[, c(names(my_data_all)[1], data_split$valid), with = F]
    my_val_data_em <- my_val_data[genes %in% my_vals_em$genes]
    my_val_data <- my_val_data[genes %in% my_vals$genes]
    
    
    
    # Calculate lamba using validation dataset instead
    if(bulk){
      my_vals_n <- params$mu_c
      my_vals_n <- my_vals_n[genes %in% my_vals$genes,]
      
      # mu total rna
      mean_rna <- apply(my_val_data[,-1], 1, mean)
      denom <- 1/mean_rna
      
      my_vals <- data.table(genes = my_vals_n$genes,
                            my_vals_n[, -1] * denom)
      
      # Removing nan genes
      nan_genes <- my_vals_n$genes[is.nan(denom)]
      nan_idx <- apply(my_vals[,-1], 2, function(x) {which(is.nan(x))})
      nan_idx <- unique(do.call(c, nan_idx))
      nan_genes <- c(nan_genes, my_vals_n$genes[nan_idx])
      if(length(nan_genes) > 0){
        my_vals <- my_vals[!genes %in% nan_genes, ]
        my_val_data <- my_val_data[!genes %in% nan_genes,]
      }
      
    }
    
    # get data for EM
    y <- apply(my_val_data[,-1], 1, mean)
    my_val_data_em <- apply(my_val_data_em[,-1], 2, mean)
    
    
    # EM approach
    means <- apply(my_vals_em[,-1], 2, mean)
    sds <- apply(my_vals_em[,-1], 2, sd)
    wait1 <- normalmixEM(my_val_data_em, lambda = .2, mu = means, sigma =  sds, arbvar = T, maxrestarts = 1e3, maxit = 5e3)
    
    if(yes_plot){
      plot(wait1, whichplots = 2, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8, lwd2 =6,
           main2="Density of gene expression levels", xlab2="Mean expression levels")
    }
    
    # summary(wait1)
    
    wait1 <- wait1[c("lambda", "mu", "sigma")]
    
    # LM
    n <- ncol(my_vals)
    
    # Obtain all columns except the last
    mysubset <- lapply(names(my_vals)[2:(n-1)], function(x) {c(x, names(my_vals)[n])})
    
    if(scaling){
      
      s <- apply(my_vals[, -1], 2, sd)
      s.max <- s[s == max(s)]
      # n.max <- which(names(s) == names(s.max))
      n.max <- which(names(my_vals) == names(s.max))
      # print( names(s.max))
      r <- 2:n
      r <- r[!r %in% n.max]
      
      
      
      # Obtain all columns except the max
      mysubset <- lapply(names(my_vals)[r], function(x) {c(x, names(my_vals)[n.max])})
      
      # print(c(r, n, n.max, names(my_vals), mysubset[[1]]))
      
      
      # x <- my_vals[, lapply(mysubset, function(x) { ( get(x[1])/sd(get(x[1])) ) - ( get(x[2])/sd(get(x[2])) ) }  )]
      # x <- my_vals[, lapply(mysubset, function(x) get(x[1]) - get(x[2])  )]
      x <- my_vals[, lapply(mysubset, function(x) scale(get(x[1])) - scale(get(x[2])) )]
      
      
      # y_star <- {y/sd(y)} -  { unlist(my_vals[,n.max, with = F])/sd(unlist(my_vals[,n.max, with = F])) }
      # y_star <- y - unlist(my_vals[,n.max, with = F]) 
      y_star <- scale(y) - scale( unlist(my_vals[,n.max, with = F]))
      
    }else{
      
      x <- my_vals[, lapply(mysubset, function(x) get(x[1]) - get(x[2])  )]
      y_star <- (y - unlist(my_vals[,n, with = F]))
      
      
    }
    
    # x <- data.matrix(cbind(1,x))
    x <- scale(x, scale = F)
    x <- cbind(1, x)
    y_star <- scale(y_star, scale = F)
    
    
    
    
    
    ft <- lm(y_star~ x[,-1])
    # beta <- solve(crossprod(x)) %*% crossprod(x, y_star)
    beta <- coef(ft)
    beta <- c(beta)
    
    beta <- c(beta[-1], 1 - sum(beta[-1]))
    names(beta) <- names(my_vals)[2:(n)]
    
    if (scaling){
      names(beta) <- names(my_vals)[c(r, n.max)]
    }
    
    
    
    
    x <- my_vals[, lapply(mysubset, function(x) get(x[1]) - get(x[2])  )]
    if(scaling){
      y_star <- (y - unlist(my_vals[,n.max, with = F]) )
    }else{
      y_star <- (y - unlist(my_vals[,n, with = F]) )
    }
    x <- data.matrix(cbind(1,x))
    
    
    
    
    mod1 <- penalized(scale(y_star, scale = F), scale(x[, -1], scale = F), ~0,lambda1=.9, positive = T)
    const.beta <- coef(mod1)
    const.beta <- c(const.beta, 1-sum(const.beta))
    
    if (length(const.beta) != (n-1)){
      mod1 <- lm(scale(y_star, scale = F) ~ scale(x[, -1], scale = F))
      const.beta <- coef(mod1)
      const.beta[const.beta < 0 ] <- 0
      const.beta <- c(const.beta[-1], 1-sum(const.beta))
    }
    
    
    names(const.beta) <- names(my_vals)[2:(n)]
    if (scaling){
      names(const.beta) <- names(my_vals)[c(r, n.max)]
    }
    
    # LM2
    
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
    # y_n <- y
    # my_dt <- data.table(y = y_n, params_out_l[, -1]*y_n)
    
    n <- ncol(my_dt)
    # Obtain all columns except the last
    mysubset <- lapply(names(my_dt)[1:(n-1)], function(x) {c(x, names(my_dt)[n])})
    my_dt <- my_dt[, lapply(mysubset, function(x) scale(get(x[1]), scale = F) - scale(get(x[2]), scale = F)  )]
    
    ft <- lm(V1 ~ -1 + ., data = my_dt)
    b <- coef(ft)
    b <- c(b, 1-sum(b))
    names(b) <- names(params_out_l)[-1]
    
    
    
    
    
    
    # Truth 
    # all_pi_c <- unique(est_ref_params(my_data_all, cell_type_info)$pi_c)
    all_pi_c <- unique(est_ref_params(my_data_all[, c(names(my_data_all)[1], data_split$valid), with = F], cell_type_info)$pi_c) 
    wait_out <- wait1$`lambda`
    names(wait_out) <- names(wait1$`mu`)
    
    a <- names(my_vals)[-1]
    # bn <- names(beta)
    # print(c(match(a,bn), a, bn))
    # print(a)
    beta <- beta[match(a,names(beta))]
    b <- b[match(a,names(b))]
    
    # if(class(ft2) != "try-error"){
    #   beta2 <- beta2[match(a,names(beta2))]
    # }
    
    const.beta <- const.beta[match(a,names(const.beta))]
    
    
    
    # if (!b | class(ft2) != "try-error"){
    #   re <- list(EM = wait_out, LM = beta, const.LM = const.beta, truth = all_pi_c)
    #   re <- do.call(rbind, re)
    #   re <- data.table(method = c("EM", "LM", "constrd.LM", "Truth"), re)
    # }else{
    #   re <- list(EM = wait_out, LM = beta, WLS = beta2, const.LM = const.beta ,truth = all_pi_c)
    #   re <- do.call(rbind, re)
    #   re <- data.table(method = c("EM", "LM", "WLS", "constrd.LM", "Truth"), re)
    # }
    
    re <- list(EM = wait_out, LM = beta, const.LM = const.beta, LM2 = b,truth = all_pi_c)
    re <- do.call(rbind, re)
    re <- data.table(method = c("EM", "LM", "constrd.LM", "LM2", "Truth"), re)
  }
  
  # Testing
  if(test > 0){
    
    # EM
    {
      
      # test_dt <- my_data_all[, c(names(my_data_all)[1], data_split$test), with = F]
      # test_dt <- test_dt[genes %in% my_vals_em$genes]
      # cell_type_data_test <- cell_type_info[cells %in% colnames(test_dt)[-1],]
      # cell_type_data_test <- cell_type_data_test[match(colnames(test_dt)[-1], cells),]
      # names(cell_type_data_test)[2] <- "Neurons"
      # cell_type_data_test2 <- cell_type_data_test
      # if (ncol(cell_type_data_test2) > 2){
      #   out_col <- 3:ncol(cell_type_data_test2)
      #   out_col <- names(cell_type_data_test2)[-out_col]
      #   cell_type_data_test <- reshape2::dcast(cell_type_data_test[,out_col, with = F], cells~Neurons)
      # }else{
      #   cell_type_data_test <- reshape2::dcast(cell_type_data_test, cells~Neurons)
      # }
      # 
      # 
      # cell_type_data_test <- data.table(cells = cell_type_data_test[,1], 
      #                                   apply(cell_type_data_test[,-1], 2, function(x){
      #                                     x <- as.numeric(!is.na(x))
      #                                   }))
      # 
      # 
      # y_bar <- apply(test_dt[, -1], 1, mean)
      # 
      # 
      # y_em <- apply(test_dt[, -1], 1, function(x){
      #   
      #   r <- x * (data.matrix(cell_type_data_test[, -1]) %*% wait_out)
      #   r <- unlist(r)
      #   a <- tapply(r, cell_type_data_test2$Neurons, mean)
      #   # a <- data.frame(Neurons = cell_type_data_test2$Neurons, "k" = r) %>% data.table
      #   # a[, k:= r]
      #   # a[, .(s = mean(k)), by = Neurons]
      #   # r <- sum(a$s)
      #   r <- sum(a)
      #   return(r)
      #   
      #   
      # })
      # y_em <- y_bar - y_em
      # y_em <- t(y_em) %*% y_em
      # y_em <- sqrt(y_em)
      
      y_em <- NA
    }
    
    # LM
    {
      test_dt <- my_data_all[, c(names(my_data_all)[1], data_split$test), with = F]
      test_dt <- test_dt[genes %in% my_vals$genes]
      cell_type_data_test <- cell_type_info[cells %in% colnames(test_dt)[-1],]
      cell_type_data_test <- cell_type_data_test[match(colnames(test_dt)[-1], cells),]
      names(cell_type_data_test)[2] <- "Neurons"
      cell_type_data_test2 <- cell_type_data_test
      
      if (ncol(cell_type_data_test2) > 2){
        out_col <- 3:ncol(cell_type_data_test2)
        out_col <- names(cell_type_data_test2)[-out_col]
        cell_type_data_test <- reshape2::dcast(cell_type_data_test[,out_col, with = F], cells~Neurons)
      }else{
        cell_type_data_test <- reshape2::dcast(cell_type_data_test, cells~Neurons)
      }
      
      cell_type_data_test <- data.table(cells = cell_type_data_test[,1], 
                                        apply(cell_type_data_test[,-1], 2, function(x){
                                          x <- as.numeric(!is.na(x))
                                        }))
      y_bar <- apply(test_dt[, -1], 1, mean)
      
      nt <- ncol(cell_type_data_test[, -1])
      if(nt == length(beta)){
        y_lm <- apply(test_dt[, -1], 1, function(x){
          
          r <- x * (data.matrix(cell_type_data_test[, -1]) %*% beta)
          r <- unlist(r)
          a <- tapply(r, cell_type_data_test2$Neurons, mean)
          r <- sum(a)
          return(r)
          
        })
        y_lm <- y_bar - y_lm
        y_lm <- t(y_lm) %*% y_lm
        y_lm <- sqrt(y_lm)
        
        y_lm2 <- apply(test_dt[, -1], 1, function(x){
          
          r <- x * (data.matrix(cell_type_data_test[, -1]) %*% b)
          r <- unlist(r)
          a <- tapply(r, cell_type_data_test2$Neurons, mean)
          r <- sum(a)
          return(r)
          
        })
        y_lm2 <- y_bar - y_lm2
        y_lm2 <- t(y_lm2) %*% y_lm2
        y_lm2 <- sqrt(y_lm2)
        
        y_clm <- apply(test_dt[, -1], 1, function(x){
          
          r <- x * (data.matrix(cell_type_data_test[, -1]) %*% const.beta)
          r <- unlist(r)
          a <- tapply(r, cell_type_data_test2$Neurons, mean)
          r <- sum(a)
          return(r)
          
        })
        
        y_clm <- y_bar - y_clm
        y_clm <- t(y_clm) %*% y_clm
        y_clm <- sqrt(y_clm)
        
      }else{
        y_lm <- y_clm <- y_lm2 <- NA
      }
      
    }
    
    re_2 <- c("EM" = y_em, "LM" = y_lm, "Const.LM" = y_clm, "LM2" = y_lm2)
    
    re <- list("estimates" = re, "Sqrt residual error sum of squares" = re_2)
  }
  
  
  
  return(re)
}





# Two step approach --------------------------------------------------------
# Generate/load data
rna_gene_x_cells <- fread(file.path(".", "cell_data", "rna_gene_x_cells.txt"))
cell_type_data <- fread(file.path(".", "cell_data", "cell_type_info.txt"))
cell_type_data <- cell_type_data[ all %in% c("Neurons", "OPC",  "Astrocytes",
                                             "Oligodendrocytes", "Microglia", "Endothelial"), ]

my_vals <- est_ref_params(rna_gene_x_cells, cell_type_data)
pv <- my_vals$p_v

my_sim_p <- function(p, use.sva = T, iter = 1, tol = 1e-6, scaling = T, eval = T){
  # Subset to informative genes
  my_params <- my_vals$out
  
  # my_params <- my_params[genes %in% pv[prop.v_e > 0.1, genes], ]
  if(p == 0){
    params_out_m <- my_vals$mean_rna
    params_out_m <- params_out_m[genes %in% pv[prop.v_e >= p, genes], ]
    params_out_l <- my_vals$lambda
    params_out_l <- params_out_l[genes %in% pv[prop.v_e >= p, genes], ]
    # my_params <- my_params[genes %in% pv[prop.v_e >= p, genes], ]
    
  }else{
    params_out_m <- my_vals$mean_rna
    params_out_m <- params_out_m[genes %in% pv[prop.v_e > p, genes], ]
    params_out_l <- my_vals$lambda
    params_out_l <- params_out_l[genes %in% pv[prop.v_e > p, genes], ]
    # my_params <- my_params[genes %in% pv[prop.v_e > p, genes], ]
  }
  
  e <- new.env()
  load("./NAc_Nicotine_hg38_rseGene_rawCounts_postQCSamples_n223.rda", envir = e)
  load("./MatchedCellComp/singleCell_iPSC_quake_coefEsts_calibration_Zscale_adultOnly.rda", envir = e)
  attach(e, warn.conflicts = F)
  yExprs <- log2(getRPKM(rse_gene, "Length")+1)
  
  
  # project
  yExprs <- yExprs[rownames(coefEsts),]
  yExprs_Z <- scale(yExprs[rownames(coefEsts),])
  
  
  # my_params <- my_params[genes %in% rownames(coefEsts),]
  
  gg <- rownames(coefEsts)
  gg2 <- gg[gg %in% params_out_m$genes]
  # gg <- gg[gg %in% my_params$genes]
  
  # my_params <- my_params[match(gg, genes),]
  params_out_m <- params_out_m[match(gg2, genes),]
  params_out_l <- params_out_l[match(gg2, genes),]
  
  
  # analysis_dt <- which(row.names(yExprs) %in% gg)
  # analysis_dt <- yExprs[analysis_dt, ]
  analysis_dt2 <- which(row.names(yExprs) %in% gg2)
  analysis_dt2 <- yExprs[analysis_dt2, ]
  
  # Get initial estimates of pi
  sample_ests_LM2 <- apply(analysis_dt2, 2, function(y){
    
    # LM
    n <- ncol(my_params)
    
    
    y_n <- y/params_out_m$V1
    my_dt <- data.table(y = y_n, params_out_l[, -1])
    n <- ncol(my_dt)
    # Obtain all columns except the last
    mysubset <- lapply(names(my_dt)[1:(n-1)], function(x) {c(x, names(my_dt)[n])})
    my_dt <- my_dt[, lapply(mysubset, function(x) scale(get(x[1]), scale = F) - scale(get(x[2]), scale = F)  )]
    
    ft <- lm(V1 ~ -1 + ., data = my_dt)
    b <- coef(ft)
    b <- c(b, 1-sum(b))
    names(b) <- names(params_out_l)[-1]
    if(any(b < 0)){
      b[b <0] <- 1e-3
      b[b > 1]<- 1
    }
    b <- b/sum( b )
    
    
    return(b)
    
  })
  sample_ests_LM2 <- t(sample_ests_LM2)
  sample_ests_LM2_old <- sample_ests_LM2
  # head(sample_ests_LM2)
  # summary(sample_ests_LM2)
  
  if(scaling){
    analysis_dt2 <- t(scale(t(analysis_dt2)))
  }
  
  # SSVA to eliminate noise
  x <- sample_ests_LM2[, -1]
  x <- cbind(1, x)
  if(use.sva){
    mod1 <-  x
    mod0 <- cbind(mod1[,1])
    analysis_dt2_sva <- apply(analysis_dt2, 2, function(x) {exp(x) })
    svseq <- svaseq(analysis_dt2_sva,mod1,mod0, constant = 0)$sv
  }
  
  
  
  gene_ests_LM2 <- apply(analysis_dt2, 1, function(y){
    
    # LM
    n <- ncol(sample_ests_LM2)
    if(use.sva){
      ft <- lm(y ~ -1 + x + svseq)
    }else{
      ft <- lm(y ~ -1 + x )
    }
    
    beta <- coef(ft)
    beta <- beta[1:n]
    beta_1 <- sapply(beta[-1], function(b) {b <- b + beta[1]; return(b)} )
    b <- c(beta[1], beta_1)
    names(b) <- colnames(sample_ests_LM2)
    return(b)
    
  })
  gene_ests_LM2 <- t(gene_ests_LM2)
  # head(gene_ests_LM2)
  # summary(gene_ests_LM2)
  
  # Get new estimates of pi
  sample_ests_LM2 <- apply(analysis_dt2, 2, function(y){
    
    # LM
    y_n <- y
    my_dt <- data.table(y = y_n, gene_ests_LM2)
    n <- ncol(my_dt)
    # Obtain all columns except the last
    mysubset <- lapply(names(my_dt)[1:(n-1)], function(x) {c(x, names(my_dt)[n])})
    my_dt <- my_dt[, lapply(mysubset, function(x) scale(get(x[1]), scale = F) - scale(get(x[2]), scale = F)  )]
    
    ft <- lm(V1 ~ -1 + ., data = my_dt)
    b <- coef(ft)
    b <- c(b, 1-sum(b))
    names(b) <- names(params_out_l)[-1]
    if(any(b < 0)){
      b[b <0] <- 1e-3
      b[b > 1]<- 1
    }
    b <- b/sum( b )
    
    return(b)
    
  })
  sample_ests_LM2 <- t(sample_ests_LM2)
  # head(sample_ests_LM2)
  # summary(sample_ests_LM2)
  
  if(iter > 1){
    test_est <- sample_ests_LM2_old - sample_ests_LM2
    sample_ests_LM2_out <- sample_ests_LM2
    test_est <- apply(test_est, 2, function(x){
      x <- crossprod(x)/length(x)
      return(x)
    })
    err <- all(test_est < tol)
    while (!err) {
      sample_ests_LM2_old <- sample_ests_LM2
      x <- sample_ests_LM2[, -1]
      x <- cbind(1, x)
      if(use.sva){
        mod1 <-  x
        mod0 <- cbind(mod1[,1])
        analysis_dt2_sva <- apply(analysis_dt2, 2, function(x) {exp(x) })
        svseq <- svaseq(analysis_dt2_sva,mod1,mod0, constant = 1e-5)$sv
      }
      gene_ests_LM2 <- apply(analysis_dt2, 1, function(y){
        
        # LM
        n <- ncol(sample_ests_LM2)
        if(use.sva){
          ft <- lm(y ~ -1 + x + svseq)
        }else{
          ft <- lm(y ~ -1 + x )
        }
        ft <- lm(y ~ -1 + x )
        beta <- coef(ft)
        beta <- beta[1:n]
        beta_1 <- sapply(beta[-1], function(b) {b <- b + beta[1]; return(b)} )
        b <- c(beta[1], beta_1)
        names(b) <- colnames(sample_ests_LM2)
        return(b)
        
      })
      gene_ests_LM2 <- t(gene_ests_LM2)
      sample_ests_LM2 <- apply(analysis_dt2, 2, function(y){
        
        # LM
        y_n <- y
        my_dt <- data.table(y = y_n, gene_ests_LM2)
        n <- ncol(my_dt)
        # Obtain all columns except the last
        mysubset <- lapply(names(my_dt)[1:(n-1)], function(x) {c(x, names(my_dt)[n])})
        my_dt <- my_dt[, lapply(mysubset, function(x) scale(get(x[1]), scale = F) - scale(get(x[2]), scale = F)  )]
        
        ft <- lm(V1 ~ -1 + ., data = my_dt)
        b <- coef(ft)
        b <- c(b, 1-sum(b))
        names(b) <- names(params_out_l)[-1]
        if(any(b < 0)){
          b[b <0] <- 1e-3
          b[b > 1]<- 1
        }
        b <- b/sum( b )
        
        return(b)
        
      })
      sample_ests_LM2 <- t(sample_ests_LM2)
      test_est <- sample_ests_LM2_old - sample_ests_LM2
      test_est <- apply(test_est, 2, function(x){
        x <- crossprod(x)/length(x)
        return(x)
      })
      err <- all(test_est < tol)
      if(err){
        sample_ests_LM2_out <<- sample_ests_LM2
      }
      # print(test_est)
    }
    sample_ests_LM2 <<- sample_ests_LM2_out
    
  }
  
  # sample_ests_LM2 <- apply(sample_ests_LM2, 2, round, 4)
  sample_ests_LM2 <- data.table(samples = row.names(sample_ests_LM2), sample_ests_LM2)
  sample_ests_LM2[, RNum := sapply(strsplit(samples, "_"), "[", 1)]
  

  
  if(!eval){
    out <- sample_ests_LM2
  }else{
    load("NAc_rse_gene_withCompEsts.rda", envir = e)
    attach(e, warn.conflicts = F)
    sample_ests_LM2 <- sample_ests_LM2[match(rse_gene$RNum, RNum),]
    
    rmse_dnam <- rse_gene$NeuN_neg_DNAm-sample_ests_LM2$`FALSE`
    rmse_dnam <- crossprod(rmse_dnam)/length(rmse_dnam)
    rmse_dnam <- sqrt(rmse_dnam)
    
    rmse_rna <- rse_gene$NeuN_neg_RNA-sample_ests_LM2$`FALSE`
    rmse_rna <- crossprod(rmse_rna)/length(rmse_rna)
    rmse_rna <- sqrt(rmse_rna)
    
    
    out <- rbind(data.table(p, ngenes = nrow(params_out_l), pearson = cor(rse_gene$NeuN_neg_RNA , sample_ests_LM2$`FALSE`), "rmse" = rmse_rna, type = "RNA"),
                 data.table(p, ngenes = nrow(params_out_l), pearson = cor(rse_gene$NeuN_neg_DNAm , sample_ests_LM2$`FALSE`), "rmse" = rmse_dnam, type = "DNAm"))
  }
  return(out)
  
}

# With scaling 
pv_sims <- lapply(seq(0,.7, by = 1e-2), my_sim_p, iter = 2)
pv_sims <- do.call(rbind, pv_sims)
pv_sims[, `:=`("scale" = T, "SVA" = T)]
fwrite(pv_sims, file.path(".", "sim_results", "two_stage_sva_scale.txt"))

pv_sims_no_sva <- lapply(seq(0,.7, by = 1e-2), my_sim_p, use.sva = F, iter = 2)
pv_sims_no_sva <- do.call(rbind, pv_sims_no_sva)
pv_sims_no_sva[, `:=`("scale" = T, "SVA" = F)]
fwrite(pv_sims_no_sva, file.path(".", "sim_results", "two_stage_nsva_scale.txt"))

# No scaling
pv_sims_no_scale <- lapply(rev(seq(0,.7, by = 1e-2)), my_sim_p, iter = 2, scaling = F)
pv_sims_no_scale <- do.call(rbind, pv_sims_no_scale)
pv_sims_no_scale[, `:=`("scale" = F, "SVA" = T)]
fwrite(pv_sims_no_scale, file.path(".", "sim_results", "two_stage_sva_nscale.txt"))

pv_sims_no_sva_no_scale <- lapply(seq(0,.7, by = 1e-2), my_sim_p, use.sva = F, iter = 2, scaling = F)
pv_sims_no_sva_no_scale <- do.call(rbind, pv_sims_no_sva_no_scale)
pv_sims_no_sva_no_scale[, `:=`("scale" = F, "SVA" = F)]
fwrite(pv_sims_no_sva_no_scale, file.path(".", "sim_results", "two_stage_nsva_nscale.txt"))

pv_all <- rbind(pv_sims, pv_sims_no_sva, pv_sims_no_scale, pv_sims_no_sva_no_scale)
fwrite(pv_all, file.path(".", "sim_results", "two_stage_sva_sims_lm_all.txt"))
fwrite(pv_sims, file.path(".", "sim_results", "two_stage_sva_sims_lm.txt"))
fwrite(pv_sims_no_sva, file.path(".", "sim_results", "two_stage_sims_lm.txt"))


pv_all <- fread( file.path(".", "sim_results", "two_stage_sva_sims_lm_all.txt"))
pv_sims <- fread( file.path(".", "sim_results", "two_stage_sva_sims_lm.txt"))
pv_sims_no_sva <- fread(file.path(".", "sim_results", "two_stage_sims_lm.txt"))
# plot --------------------------------------------------------------------


transparent_legend =  theme(
  legend.background = element_rect(fill ="transparent"),
  legend.key = element_rect(fill = "transparent",
                            color = "transparent")
)

remove_grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))

no_x_axis_label <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

cols <- RColorBrewer::brewer.pal(length(unique(pv_all$type)), "Dark2")
cols <- cols[1:2]
names(cols) <- unique(pv_all$type)

ggplot(data = pv_all[scale==T & SVA== T,], aes(x = p, y = pearson, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  ylab("Pearson correlation coefficient") + xlab("Cutoff") +
  ggtitle("Relationship between cutoff and model fit\n(Two Stage approach: Scaled and SVA)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.25,0.35),
        legend.title = element_blank())
ggsave(file.path(".", "figs","effect_of_cutoff_two_step.png"), dpi = "retina")

plot_dt <- pv_all
plot_dt[, `:=`(scale = ifelse(scale, "Scaled", "Unscaled"), 
               SVA = ifelse(SVA, "SVA", "No SVA"))]
cols <- RColorBrewer::brewer.pal(4, "Dark2")
names(cols) <- unique(with(plot_dt, interaction(scale, SVA)))

min(plot_dt$rmse.V1*100)

ggplot(data = plot_dt, aes(x = p, y = rmse.V1*100, color = interaction(scale, SVA), linetype = interaction(scale, SVA))) + 
  geom_point() +
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  facet_grid(.~type, scales = "free_x") +
  ylab(expression(RMSE%*%100)) + xlab("Cutoff") +
  # ylim(4,20) +
  ggtitle("Relationship between cutoff and model fit\n(Two Stage approach)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave(file.path(".", "figs","rmse_vs_cutoff.png"), dpi = "retina")


ggplot(data = plot_dt, aes(x = p, y = rmse.V1*100, color = interaction(scale, SVA), linetype = interaction(scale, SVA))) + 
  # geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  facet_grid(.~type, scales = "free_x") +
  ylab(expression(RMSE%*%100)) + xlab("Cutoff") +
  ylim(4,20) +
  ggtitle("Relationship between cutoff and model fit\n(Two Stage approach)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave(file.path(".", "figs","rmse_vs_cutoff_trim.png"), dpi = "retina")


sample_ests_LM2 <- my_sim_p(p = .2, use.sva = T, scaling = T, iter = 2, eval = F)
load("NAc_rse_gene_withCompEsts.rda")

png(file.path(".", "figs", "ts_houseman.png"), width=10, height=10, units="in", res=300)
palette(brewer.pal(4,"Dark2"))
par(mfrow = c(1,2), pty = "m")
r_2_size <- 1.2
r_2 <- cor(rse_gene$NeuN_neg_RNA, sample_ests_LM2$`FALSE`)^2
r_2 <- round(r_2, 2)
plot(rse_gene$NeuN_neg_RNA , sample_ests_LM2$`FALSE`,
     xlab = "RNA-based Houseman", ylab = "RNA-based TS",
     main = "Houseman RNA VS TS",
     pch = 21, bg = 1, ylim = c(0,1),xlim=c(0,1))
abline(0,1,lty=2,lwd=2)
points(rse_gene$Neurons_RNA , sample_ests_LM2$`TRUE`,
       pch = 21, bg = 2)
legend("topleft", c("NeuN-", "NeuN+"), pch = 15, col=1:2, horiz = T, cex = .8, text.width = .15)
text(.8, .5, substitute(paste(R^2, nn), list(nn=paste("=", r_2))) ,
     cex = r_2_size)

r_2 <- cor(rse_gene$NeuN_neg_DNAm , sample_ests_LM2$`FALSE`)^2
r_2 <- round(r_2, 2)
plot(rse_gene$NeuN_neg_DNAm , sample_ests_LM2$`FALSE`,
     ylab = "RNA-based TS", xlab = "DNAm-based", 
     main = "Houseman Methylation VS TS",
     pch = 21, bg = 1, ylim = c(0,1),xlim=c(0,1))
abline(0,1,lty=2,lwd=2)
points(rse_gene$NeuN_pos_DNAm,  sample_ests_LM2$`TRUE`,
       pch = 21, bg = 2)
legend("topleft", c("NeuN-", "NeuN+"), pch = 15, col=1:2, horiz = T, cex = .8, text.width = .15)
text(.8, .5, substitute(paste(R^2, nn), list(nn=paste("=", r_2))) ,
     cex = r_2_size)
dev.off()