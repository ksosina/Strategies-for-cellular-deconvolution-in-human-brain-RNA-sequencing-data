
# Functions + testing -----------------------------------------------------


# Application of a weigthed scheme
# Change the convergence criterion to the mean RNA level instead of the cell type proportion
all_objs <- ls()
all_objs <- all_objs[!all_objs %in% c("mixes_abbas", "cell_type_data_var", "cell_type_data")]
rm(list = all_objs)
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
my_sim_p_abbas3 <- function(p, use.sva = F, iter = 1e3, tol = 1e-12, eval = T, use.sum = NULL, scaling = T, prop_data = cell_type_data_var, 
                            all_dt = mixes_abbas, raw_dt = cell_type_data, norm = T, method = "TMM", cov_data = NULL){
  # require(nnls)
  require(sva)
  require(edgeR)
  require(data.table)
  require(dplyr)
  # Subset to informative genes from reference
  
  # TODO: Clean this part up
  if(!is.null(use.sum)){
    if (use.sum == "pval") {
      use.pval = T; use.sva = T
    }else if (use.sum == "r2"){
      use.r2 = T; use.sva = T
    }
    
  }else{
    use.pval = F; use.r2 = F
  }
  
  # raw_dt[, V2:= NULL]
  
  if(norm){
    
    
    
    yExprs_normfactors <- calcNormFactors(cbind(data.matrix(all_dt[, -1]), raw_dt[, -c(1), with = F]), method = method )
    a <- ncol(all_dt)  - 1 + 1
    b <- a + ncol(raw_dt[, -c(1), with = F]) - 1
    q_out <- yExprs_normfactors[a:b]
    q_out <- replicate(nrow(all_dt)  , q_out)
    q_out <- t(q_out)
    raw_dt <- cbind(raw_dt[, 1], raw_dt[, -c(1), with = F]*q_out)
    
    yExprs_normfactors <- yExprs_normfactors[-c(a:b)]
    yExprs_normfactors <- replicate(nrow(all_dt)  , yExprs_normfactors)
    yExprs_normfactors <- t(yExprs_normfactors)
    all_dt <- cbind(all_dt[, 1], all_dt[,-1]*yExprs_normfactors)
  }
  
  
  
  # Estimate
  cutoff <- quantile(prop_data$prop, probs = c(1-p,p))
  cutoff[2] <- ifelse(p == 0, 0, cutoff[2])
  info_idx <- which(prop_data$prop >= cutoff[2] )
  
  if(iter == 0){
    sample_ests_LM2 <- apply(all_dt[, -1], 2, function(y){
      
      # LM
      y <- y[info_idx]
      my_params <- raw_dt[info_idx, -c(1), with = F]
      n <- ncol(my_params)
      
      # Obtain all columns except the last
      mysubset <- lapply(names(my_params)[1:(n-1)], function(x) {c(x, names(my_params)[n])})
      x <- my_params[, lapply(mysubset, function(x) get(x[1])-get(x[2]))]
      x <- data.matrix(cbind(1,x))
      
      y_star <- (y - unlist(my_params[,n, with = F]))
      
      
      
      beta <- solve(crossprod(x)) %*% crossprod(x, y_star)
      beta <- c(beta)
      
      beta <- c(beta[-1], 1 - sum(beta[-1]))
      names(beta) <- names(my_params)[1:(n)]
      if(any(beta < 0)){
        beta[beta < 0] <- 1e-5
        # beta[beta > 0] <- beta[beta >0] /sum( beta[beta >0] )
        beta <- beta /sum( beta )
      }
      
      
      return(beta)
      
    })
    sample_ests_LM2 <- t(sample_ests_LM2)
  }else{
    
    
    
    # Get initial estimates of pi
    sample_ests_LM2 <- apply(all_dt[, -1], 2, function(y){
      
      # LM
      y <- y[info_idx]
      my_params <- raw_dt[info_idx, -c(1), with = F]
      n <- ncol(my_params)
      
      # Obtain all columns except the last
      mysubset <- lapply(names(my_params)[1:(n-1)], function(x) {c(x, names(my_params)[n])})
      x <- my_params[, lapply(mysubset, function(x) get(x[1])-get(x[2]))]
      x <- data.matrix(cbind(1,x))
      
      y_star <- (y - unlist(my_params[,n, with = F]))
      
      
      
      beta <- solve(crossprod(x)) %*% crossprod(x, y_star)
      beta <- c(beta)
      
      beta <- c(beta[-1], 1 - sum(beta[-1]))
      names(beta) <- names(my_params)[1:(n)]
      if(any(beta < 0)){
        beta[beta < 0] <- 1e-5
        # beta[beta > 0] <- beta[beta >0] /sum( beta[beta >0] )
        beta <- beta /sum( beta )
      }
      
      
      return(beta)
      
    })
    sample_ests_LM2 <- t(sample_ests_LM2)
    sample_ests_LM2_old <- sample_ests_LM2
    
    # Scale
    if(scaling){
      w_mixes_abbas <- t(scale(t(data.matrix(all_dt[,-1]))))
    }
    
    
    
    # SSVA to eliminate noise
    x <- sample_ests_LM2[, -1]
    x <- cbind(1, x)
    if(!is.null(cov_data)){
      x <- cbind(x, cov_data)
    }
    if(use.sva){
      mod1 <-  x
      mod0 <- cbind(mod1[,1])
      if (scaling){
        analysis_dt2_sva <- apply(w_mixes_abbas, 2, function(x) {exp(x) })
        svseq <- svaseq(analysis_dt2_sva,mod1,mod0, constant = 0)$sv
      }else{
        svseq <- svaseq(data.matrix(all_dt[,-1]),mod1,mod0, constant = 1)$sv
      }
      if( all(svseq == 0) ){
        svseq <- rep(0, nrow(x))
      }
      
    }
    
    if (scaling){
      if(use.sva){
        gene_ests_LM2 <- lm(data.matrix(t(w_mixes_abbas)) ~ 0 + x + svseq)
        gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
        gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
        colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
      }else{
        
        gene_ests_LM2 <- lm(data.matrix(t(w_mixes_abbas)) ~ 0 + x)
        gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
        gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
        colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
        
      }
      
      
    }else{
      if(use.sva){
        gene_ests_LM2 <- lm(data.matrix(t(all_dt[,-1])) ~ 0 + x + svseq)
        gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
        gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
        gene_ests_LM2 <- apply(gene_ests_LM2, 2, function(b){ b[b <= 0] <- 0; return(b) })
        colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
      }else{
        gene_ests_LM2 <- lm(data.matrix(t(all_dt[,-1])) ~ 0 + x)
        gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
        gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
        colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
        gene_ests_LM2 <- apply(gene_ests_LM2, 2, function(b){ b[b <= 0] <- 0; return(b)})
      }
      
      
    }
    
    
    
    
    
    # Get new estimates of pi
    if(scaling){
      sample_ests_LM2 <- apply(w_mixes_abbas, 2, function(y){
        
        # LM
        y <- y[info_idx]
        my_params <- data.table(gene_ests_LM2[info_idx, ]) #cell_type_data_mu[info_idx, -1]
        n <- ncol(my_params)
        
        # Obtain all columns except the last
        mysubset <- lapply(names(my_params)[1:(n-1)], function(x) {c(x, names(my_params)[n])})
        x <- my_params[, lapply(mysubset, function(x) get(x[1])-get(x[2]))]
        x <- data.matrix(cbind(1,x))
        
        y_star <- (y - unlist(my_params[,n, with = F]))
        
        # ft <- lm(y_star ~ -1 + x, weights = cell_type_data_var$w_2[info_idx])
        ft <- lm(y_star ~ -1 + x)
        # beta <- solve(crossprod(x)) %*% crossprod(x, y_star)
        beta <- coef(ft)
        beta <- c(beta)
        # if(scaling){
        #   beta <- c(beta, 1 - sum(beta))
        #   names(beta) <- names(my_vals)[2:(n)]
        # }
        beta <- c(beta[-1], 1 - sum(beta[-1]))
        names(beta) <- names(my_params)[1:(n)]
        if(any(beta < 0)){
          # beta[beta < 0] <- 0
          beta[beta < 0] <- 1e-5
          # beta[beta > 0] <- beta[beta >0] /sum( beta[beta >0] )
          beta <- beta /sum( beta )
        }
        
        
        return(beta)
        
      })
    }else{
      sample_ests_LM2 <- apply(all_dt[,-1], 2, function(y){
        
        # LM
        y <- y[info_idx]
        my_params <- data.table(gene_ests_LM2[info_idx, ]) #cell_type_data_mu[info_idx, -1]
        n <- ncol(my_params)
        
        # Obtain all columns except the last
        mysubset <- lapply(names(my_params)[1:(n-1)], function(x) {c(x, names(my_params)[n])})
        x <- my_params[, lapply(mysubset, function(x) get(x[1])-get(x[2]))]
        x <- data.matrix(cbind(1,x))
        
        y_star <- (y - unlist(my_params[,n, with = F]))
        
        
        # if(scaling){
        #   y_star <- scale(y_star)
        #   x <- scale(x[, -1])
        # }
        
        
        # ft <- lm(y_star ~ -1 + x, weights = cell_type_data_var$w_2[info_idx])
        ft <- lm(y_star ~ -1 + x)
        # beta <- solve(crossprod(x)) %*% crossprod(x, y_star)
        beta <- coef(ft)
        beta <- c(beta)
        # if(scaling){
        #   beta <- c(beta, 1 - sum(beta))
        #   names(beta) <- names(my_vals)[2:(n)]
        # }
        beta <- c(beta[-1], 1 - sum(beta[-1]))
        names(beta) <- names(my_params)[1:(n)]
        if(any(beta < 0)){
          # beta[beta < 0] <- 0
          beta[beta < 0] <- 1e-5
          # beta[beta > 0] <- beta[beta >0] /sum( beta[beta >0] )
          beta <- beta /sum( beta )
        }
        
        
        
        return(beta)
        
      })
      
    }
    
    sample_ests_LM2 <- t(sample_ests_LM2)
    
    # Iterate until convergence
    if(iter > 1){
      gene_ests_LM2_old <- raw_dt[info_idx, colnames(sample_ests_LM2), with = F] 
      gene_ests_LM2_new <- gene_ests_LM2[info_idx, ]
      test_est_gene <- gene_ests_LM2_old - gene_ests_LM2_new
      test_est_gene <- apply(test_est_gene, 1, function(x){
        x <- crossprod(x)/length(x)
        return(x)
      })
      
      c <- 0
      test_est <- sample_ests_LM2_old - sample_ests_LM2
      sample_ests_LM2_out <- sample_ests_LM2
      test_est <- apply(test_est, 1, function(x){
        x <- crossprod(x)/length(x)
        return(x)
      })
      err <- max(test_est_gene) < tol
      while (!err ) {
        sample_ests_LM2_old <- sample_ests_LM2
        gene_ests_LM2_old <- gene_ests_LM2[info_idx, ]
        
        x <- sample_ests_LM2[, -1]
        x <- cbind(1, x)
        
        # if(use.sva){
        #   mod1 <-  x
        #   mod0 <- cbind(mod1[,1])
        #   if (scaling){
        #     analysis_dt2_sva <- apply(w_mixes_abbas, 2, function(x) {exp(x) })
        #     svseq <- svaseq(analysis_dt2_sva,mod1,mod0, constant = 1e-5)$sv
        #   }else{
        #     svseq <- svaseq(data.matrix(all_dt[,-1]),mod1,mod0, constant = 1)$sv
        #   }
        #   
        # }
        
        if(!is.null(cov_data)){
          x <- cbind(x, cov_data)
        }
        
        if (scaling){
          if(use.sva){
            gene_ests_LM2 <- lm(data.matrix(t(w_mixes_abbas)) ~ 0 + x + svseq)
            gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
            gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
            colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
          }else{
            gene_ests_LM2 <- lm(data.matrix(t(w_mixes_abbas)) ~ 0 + x)
            gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
            gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
            colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
            
          }
          
          
        }else{
          if(use.sva){
            gene_ests_LM2 <- lm(data.matrix(t(all_dt[,-1])) ~ 0 + x + svseq)
            gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
            gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
            gene_ests_LM2 <- apply(gene_ests_LM2, 2, function(b){ b[b <= 0] <- 0; return(b) })
            colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
          }else{
            gene_ests_LM2 <- lm(data.matrix(t(all_dt[,-1])) ~ 0 + x)
            gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
            gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
            colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
            gene_ests_LM2 <- apply(gene_ests_LM2, 2, function(b){ b[b <= 0] <- 0; return(b)})
          }
          
          
        }
        
        if(scaling){
          sample_ests_LM2 <- apply(w_mixes_abbas, 2, function(y){
            
            # LM
            y <- y[info_idx]
            my_params <- data.table(gene_ests_LM2[info_idx, ]) #cell_type_data_mu[info_idx, -1]
            n <- ncol(my_params)
            
            # Obtain all columns except the last
            mysubset <- lapply(names(my_params)[1:(n-1)], function(x) {c(x, names(my_params)[n])})
            x <- my_params[, lapply(mysubset, function(x) get(x[1])-get(x[2]))]
            x <- data.matrix(cbind(1,x))
            
            y_star <- (y - unlist(my_params[,n, with = F]))
            
            
            
            # ft <- lm(y_star ~ -1 + x, weights = cell_type_data_var$w_2[info_idx])
            ft <- lm(y_star ~ -1 + x)
            # beta <- solve(crossprod(x)) %*% crossprod(x, y_star)
            beta <- coef(ft)
            beta <- c(beta)
            
            beta <- c(beta[-1], 1 - sum(beta[-1]))
            names(beta) <- names(my_params)[1:(n)]
            if(any(beta < 0)){
              # beta[beta < 0] <- 0
              beta[beta < 0] <- 1e-5
              # beta[beta > 0] <- beta[beta >0] /sum( beta[beta >0] )
              beta <- beta /sum( beta )
            }
            
            return(beta)
            
          })
        }else{
          sample_ests_LM2 <- apply(all_dt[,-1], 2, function(y){
            
            # LM
            y <- y[info_idx]
            my_params <- data.table(gene_ests_LM2[info_idx, ]) #cell_type_data_mu[info_idx, -1]
            n <- ncol(my_params)
            
            # Obtain all columns except the last
            mysubset <- lapply(names(my_params)[1:(n-1)], function(x) {c(x, names(my_params)[n])})
            x <- my_params[, lapply(mysubset, function(x) get(x[1])-get(x[2]))]
            x <- data.matrix(cbind(1,x))
            
            y_star <- (y - unlist(my_params[,n, with = F]))
            
            
            
            
            beta <- solve(crossprod(x)) %*% crossprod(x, y_star)
            beta <- c(beta)
            # if(scaling){
            #   beta <- c(beta, 1 - sum(beta))
            #   names(beta) <- names(my_vals)[2:(n)]
            # }
            beta <- c(beta[-1], 1 - sum(beta[-1]))
            names(beta) <- names(my_params)[1:(n)]
            if(any(beta < 0)){
              # beta[beta < 0] <- 0
              beta[beta < 0] <- 1e-5
              # beta[beta > 0] <- beta[beta >0] /sum( beta[beta >0] )
              beta <- beta /sum( beta )
            }
            
            return(beta)
            
          })
        }
        sample_ests_LM2 <- t(sample_ests_LM2)
        gene_ests_LM2_new <- gene_ests_LM2[info_idx, ]
        
        c <- c + 1
        
        test_est <- sample_ests_LM2_old - sample_ests_LM2
        test_est <- apply(test_est, 1, function(x){
          x <- crossprod(x)/length(x)
          return(x)
        })
        
        test_est_gene <- gene_ests_LM2_old - gene_ests_LM2_new
        test_est_gene <- apply(test_est_gene, 1, function(x){
          x <- crossprod(x)/length(x)
          return(x)
        })
        
        # err <- all(test_est < tol)
        err <- max(test_est_gene) < tol
        if(err){
          sample_ests_LM2_out <<- sample_ests_LM2
        }
        if(c > iter){
          stop("convergence not achieved, increase tol")
        }
        # print(max(test_est_gene))
      }
      sample_ests_LM2 <<- sample_ests_LM2_out
      
    }
    
    
  }
  
  
  sample_ests_LM2 <- apply(sample_ests_LM2, 2, round, 4)
  sample_ests_LM2 <- data.table("ID" = sapply(strsplit(row.names(sample_ests_LM2), " "), "[", 1),
                                "Mix" = sapply(strsplit(row.names(sample_ests_LM2), " "), "[", 2),
                                sample_ests_LM2)
  
  
  
  if(eval){
    # Evaluate
    pear <- eval_mix(sample_ests_LM2)
    pear <- data.table(pear, "p" = p, ngenes = length(info_idx))
    return(pear)
  }else{
    return(sample_ests_LM2)
  }
  
  
}

# Subsampling
my_sim_p_abbas4 <- function(p = 0.1, use.sva = F, iter = 1e3, tol = ifelse(scaling, 1e-12, 1), eval = T, use.sum = NULL, scaling = T, prop_data = cell_type_data_var, 
                            all_dt = mixes_abbas, raw_dt = cell_type_data, norm = T, method = "TMM", cov_data = NULL, perc.sample = .5, r = 1e3){
  
  #   Improve estimation
  #   1.) Filter out non-informative genes I.e require r^2_g > 0.1
  
  # Estimate
  my_cutoff <- quantile(prop_data$prop, probs = c(1-p,p))
  my_cutoff[2] <- ifelse(p == 0, 0, my_cutoff[2])
  my_info_idx <- which(prop_data$prop >= my_cutoff[2] )
  
  
  #   2.) For this new set, randomly sample x% of genes and also calculate the average r^2 for the set
  n_set <- length(my_info_idx)*perc.sample
  n_set <- round(n_set)
  r <- min(r, choose(length(my_info_idx), n_set))
  r <- seq_len(r)
  
  #   3.) Do this for r iterations, and obtain the estimate of pi as the weighted average over all iterations. Weigh with the average r^2 calculated
  est_lm <- lapply(r, function(i){
    my_info_idx <- sample(my_info_idx, size = n_set)
    ave_r2 <- mean(prop_data$prop[my_info_idx])
    out <- my_sim_p_abbas3(p = 0, use.sva = use.sva, iter = iter, tol = tol, eval = F, use.sum = use.sum, scaling = scaling, prop_data = prop_data[my_info_idx,], 
                              all_dt = all_dt[my_info_idx,], raw_dt = raw_dt[my_info_idx, ], norm = norm, method = method, cov_data = cov_data)
    out <- data.table(out, ave_rsqrd = ave_r2)
  })
  est_lm <- do.call(rbind, est_lm)
  est_lm[, denom := sum(ave_rsqrd), by = .(ID, Mix)]
  # denom <- unique(est_lm$ave_rsqrd)
  # denom <- sum(denom)
  
  b <- est_lm[, .( ID,Mix, Jurkat =  (Jurkat*ave_rsqrd)/denom,
                               `IM-9` = (`IM-9`*ave_rsqrd)/denom,
                               Raji = (Raji*ave_rsqrd)/denom,
                               `THP-1` = (`THP-1`*ave_rsqrd)/denom)]
  b <- b[, lapply(.SD, sum, na.rm=TRUE), by=.(ID, Mix), .SDcols=c("Jurkat","IM-9","Raji","THP-1") ] 
  out <- b
  if(eval){
    out <- eval_mix(b)
    out <- data.table(out, B = n_set)
  }
  
  
  return(out)
  
}


old <- my_sim_p_abbas3(p = .95, use.sva = T, iter = 1e4, tol = 1e-12)
my_sim_p_abbas3(p = .95, use.sva = F, iter = 1e3)
my_sim_p_abbas4(p = .1, use.sva = T, iter = 1e4, tol = 1e-12, perc.sample = 1e-2, r = 10)
my_sim_p_abbas4(p = .1, use.sva = F, iter = 1e3)

keep_idx <- which(p.adjust(cell_type_data_var$f_test, "BH") <= 5e-2)

my_sim_p_abbas3(p = 0.2, use.sva = F, iter = 1e4, scaling = T, prop_data = cell_type_data_var[keep_idx,], 
                all_dt = mixes_abbas[keep_idx,], raw_dt = cell_type_data[keep_idx,])


# Simulations -------------------------------------------------------------

# Two-step Model avergaed version 1
p_sims_ts_with_sva_bagging <- my_sim_p_abbas4(p = .1, use.sva = T, iter = 1e4, tol = 1e-12, perc.sample = 1e-2, r = 1e3)
p_sims_ts_with_sva_bagging2 <- my_sim_p_abbas4(p = .1, use.sva = T, iter = 1e4, tol = 1e-12, perc.sample = 1e-2, r = 1e3)
# Two-step with and w/o SVA

# With scaling
end <- .999
iter <- 1e3
jump <- 1e-3


p_sims_ts <- lapply(seq(0, end, by = jump), my_sim_p_abbas3, use.sva = F, iter = iter, scaling = T)
p_sims_ts <- do.call(rbind, p_sims_ts)
p_sims_ts[, `:=`("scale" = T, "SVA" = F)]

fwrite(p_sims_ts, file.path(".", "sim_results", "p_sims_fullTS_with_filter_abbas.txt"))


p_sims_ts_with_sva <- lapply(seq(0,end, by = jump), my_sim_p_abbas3, use.sva = T, iter = iter, scaling = T)
p_sims_ts_with_sva <- do.call(rbind, p_sims_ts_with_sva)
p_sims_ts_with_sva[, `:=`("scale" = T, "SVA" = T)]
fwrite(p_sims_ts_with_sva, file.path(".", "sim_results", "fullTS_v2.txt"))


# One step Lm with/without sva
#sim across all values of p
end <- .999
iter <- 0
jump <- 1e-3

p_sims_LM1 <- lapply(seq(0, end, by = jump), my_sim_p_abbas3, use.sva = F, iter = iter, scaling = T)
p_sims_LM1 <- do.call(rbind, p_sims_LM1)

fwrite(p_sims_LM1, file.path(".", "sim_results", "p_sims_LM1_abbas.txt"))


# 2-step results at a pre-specified threshold Abbas et al

sample_ests_ts_info <- my_sim_p_abbas3(p = 0.75, use.sva = T, iter = 1e3, eval = F)
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

plot_dt_ts <- rbind(inner_join(sample_ests_ts_info[, .(Mix, Estimate = Jurkat)],  Mix[, .(Mix, Truth = Jurkat, type = "Jurkat")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = `IM-9`)],  Mix[, .(Mix, Truth = `IM-9`, type = "IM-9")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = Raji)],  Mix[, .(Mix, Truth = Raji, type = "Raji")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = `THP-1`)],  Mix[, .(Mix, Truth = `THP-1`, type = "THP-1")])) %>% data.table
with(plot_dt_ts, cor(x = Truth, y = Estimate))
fwrite(plot_dt_ts, file.path(".", "sim_results", "plot_dt_TS_abbas_v2.txt"))

# LM results at a pre-specified threshold Abbas et al

sample_ests_ts_info <- my_sim_p_abbas3(p = .75, use.sva = T, iter = 0, eval = F)
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

plot_dt_ts <- rbind(inner_join(sample_ests_ts_info[, .(Mix, Estimate = Jurkat)],  Mix[, .(Mix, Truth = Jurkat, type = "Jurkat")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = `IM-9`)],  Mix[, .(Mix, Truth = `IM-9`, type = "IM-9")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = Raji)],  Mix[, .(Mix, Truth = Raji, type = "Raji")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = `THP-1`)],  Mix[, .(Mix, Truth = `THP-1`, type = "THP-1")])) %>% data.table
with(plot_dt_ts, cor(x = Truth, y = Estimate))
fwrite(plot_dt_ts, file.path(".", "sim_results", "plot_dt_LM1_abbas.txt"))


# 2-step model averaged Abbas et al
sample_ests_ts_info <- fread(file.path(".", "sim_results", "p_sims_ts_with_sva_all.txt"))
ns <- seq(0,1, by = 1/300)[-1]
ns <- quantile(sample_ests_ts_info$p, probs = ns, type = 1)

sample_ests_ts_info <- sample_ests_ts_info[p %in% ns, ]

nl <- sample_ests_ts_info$p[sample_ests_ts_info$p %in% ns]
nl <- unique(nl)
denom <-  sum(nl)

sample_ests_ts_info <- sample_ests_ts_info[, .( ID,Mix, Jurkat =  (Jurkat*p)/denom,
                             `IM-9` = (`IM-9`*p)/denom,
                             Raji = (Raji*p)/denom,
                             `THP-1` = (`THP-1`*p)/denom)]
sample_ests_ts_info <- sample_ests_ts_info[, lapply(.SD, sum, na.rm=TRUE), by=.(ID, Mix), .SDcols=c("Jurkat","IM-9","Raji","THP-1") ] 

# truth
MixA <- c("Jurkat" = 2.5, "IM-9" = 1.25, "Raji" =  2.5,"THP-1" = 3.75)/10
MixB <-  c("Jurkat" = 0.5, "IM-9" = 3.17, "Raji" = 4.75, "THP-1" = 1.58)/10
MixC <-  c("Jurkat" = 0.1, "IM-9" = 4.95, "Raji" = 1.65, "THP-1" = 3.3)/10
MixD  <- c("Jurkat" = 0.02, "IM-9" = 3.33, "Raji" = 3.33, "THP-1" = 3.33)/10

Mix <- rbind(MixA, MixB, MixC, MixD)
Mix <- data.table("Mix" = paste0("Mix", LETTERS[1:4]), Mix)

plot_dt_ts <- rbind(inner_join(sample_ests_ts_info[, .(Mix, Estimate = Jurkat)],  Mix[, .(Mix, Truth = Jurkat, type = "Jurkat")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = `IM-9`)],  Mix[, .(Mix, Truth = `IM-9`, type = "IM-9")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = Raji)],  Mix[, .(Mix, Truth = Raji, type = "Raji")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = `THP-1`)],  Mix[, .(Mix, Truth = `THP-1`, type = "THP-1")])) %>% data.table
with(plot_dt_ts, cor(x = Truth, y = Estimate))
fwrite(plot_dt_ts, file.path(".", "sim_results", "plot_dt_TS_ave_abbas.txt"))


keep_idx <- which(p.adjust(cell_type_data_var$p_vals, "BH") <= 5e-2)

sample_ests_ts_info <- my_sim_p_abbas3(p = 0, use.sva = T, iter = 1e4, scaling = T, prop_data = cell_type_data_var[keep_idx,],
                                       all_dt = mixes_abbas[keep_idx,], raw_dt = cell_type_data[keep_idx,], eval = F)


plot_dt_ts <- rbind(inner_join(sample_ests_ts_info[, .(Mix, Estimate = Jurkat)],  Mix[, .(Mix, Truth = Jurkat, type = "Jurkat")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = `IM-9`)],  Mix[, .(Mix, Truth = `IM-9`, type = "IM-9")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = Raji)],  Mix[, .(Mix, Truth = Raji, type = "Raji")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = `THP-1`)],  Mix[, .(Mix, Truth = `THP-1`, type = "THP-1")])) %>% data.table
with(plot_dt_ts, cor(x = Truth, y = Estimate))

# Misspecification of the reference dataset -------------------------------

loop <- as.list(2:5)
g <- t(combn(2:5, 2))
g <- split(g, 1:nrow(g))
g <- c(loop, g)
p <- .75
ts_misspec_nsva <- sapply(g, function(k){
  sample_ests_ts_info <- my_sim_p_abbas3(p = p, use.sva = F, iter = 1e3, eval = F, raw_dt = cell_type_data[, -k, with = F])
  
  # truth
  MixA <- c("Jurkat" = 2.5, "IM-9" = 1.25, "Raji" =  2.5,"THP-1" = 3.75)/10
  MixB <-  c("Jurkat" = 0.5, "IM-9" = 3.17, "Raji" = 4.75, "THP-1" = 1.58)/10
  MixC <-  c("Jurkat" = 0.1, "IM-9" = 4.95, "Raji" = 1.65, "THP-1" = 3.3)/10
  MixD  <- c("Jurkat" = 0.02, "IM-9" = 3.33, "Raji" = 3.33, "THP-1" = 3.33)/10
  
  Mix <- rbind(MixA, MixB, MixC, MixD)
  Mix <- Mix[, -(k-1)]
  Mix <- data.table("Mix" = paste0("Mix", LETTERS[1:4]), Mix)
  
  r <- names(Mix)[-1]
  plot_dt_ts <- lapply(r, function(x){
    dt <- sample_ests_ts_info[, c("Mix", x), with = F]
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

ts_misspec <- sapply(g, function(k){
  sample_ests_ts_info <- my_sim_p_abbas3(p = p, use.sva = T, iter = 1e3, eval = F, raw_dt = cell_type_data[, -k, with = F])
  
  # truth
  MixA <- c("Jurkat" = 2.5, "IM-9" = 1.25, "Raji" =  2.5,"THP-1" = 3.75)/10
  MixB <-  c("Jurkat" = 0.5, "IM-9" = 3.17, "Raji" = 4.75, "THP-1" = 1.58)/10
  MixC <-  c("Jurkat" = 0.1, "IM-9" = 4.95, "Raji" = 1.65, "THP-1" = 3.3)/10
  MixD  <- c("Jurkat" = 0.02, "IM-9" = 3.33, "Raji" = 3.33, "THP-1" = 3.33)/10
  
  Mix <- rbind(MixA, MixB, MixC, MixD)
  Mix <- Mix[, -(k-1)]
  Mix <- data.table("Mix" = paste0("Mix", LETTERS[1:4]), Mix)
  
  r <- names(Mix)[-1]
  plot_dt_ts <- lapply(r, function(x){
    dt <- sample_ests_ts_info[, c("Mix", x), with = F]
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

lm_misspec <- sapply(g, function(k){
  sample_ests_ts_info <- my_sim_p_abbas3(p = p, use.sva = F, iter = 0, eval = F, raw_dt = cell_type_data[, -k, with = F])
  
  # truth
  MixA <- c("Jurkat" = 2.5, "IM-9" = 1.25, "Raji" =  2.5,"THP-1" = 3.75)/10
  MixB <-  c("Jurkat" = 0.5, "IM-9" = 3.17, "Raji" = 4.75, "THP-1" = 1.58)/10
  MixC <-  c("Jurkat" = 0.1, "IM-9" = 4.95, "Raji" = 1.65, "THP-1" = 3.3)/10
  MixD  <- c("Jurkat" = 0.02, "IM-9" = 3.33, "Raji" = 3.33, "THP-1" = 3.33)/10
  
  Mix <- rbind(MixA, MixB, MixC, MixD)
  Mix <- Mix[, -(k-1)]
  Mix <- data.table("Mix" = paste0("Mix", LETTERS[1:4]), Mix)
  
  r <- names(Mix)[-1]
  plot_dt_ts <- lapply(r, function(x){
    dt <- sample_ests_ts_info[, c("Mix", x), with = F]
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



# Model averaging ---------------------------------------------------------

# With scaling
end <- .999
start <- 1e-3
iter <- 1e3

p_sims_ts_with_sva_all <- lapply(seq(start,end, by = 1e-3), my_sim_p_abbas3, use.sva = T, iter = iter, scaling = T, eval = F)
p_sims_ts_with_sva_all <- do.call(rbind, p_sims_ts_with_sva_all)
p_sims_ts_with_sva_all[, p:=rep(seq(start,end, length.out = iter), each = 12)]
fwrite(p_sims_ts_with_sva_all, file.path(".", "sim_results", "p_sims_ts_with_sva_all.txt"))

p_sims_ts_with_sva_ave <- lapply(seq(10, 1e3, by = 10), function(k){
  
  ns <- seq(0,1, by = 1/k)[-1]
  ns <- quantile(p_sims_ts_with_sva_all$p, probs = ns, type = 1)
  
  p_sims_ts_with_sva <- p_sims_ts_with_sva_all[p %in% ns, ]
  
  nl <- p_sims_ts_with_sva_all$p[p_sims_ts_with_sva_all$p %in% ns]
  nl <- unique(nl)
  denom <-  sum(nl)
  b <- p_sims_ts_with_sva[, .( ID,Mix, Jurkat =  (Jurkat*p)/denom,
                               `IM-9` = (`IM-9`*p)/denom,
                               Raji = (Raji*p)/denom,
                               `THP-1` = (`THP-1`*p)/denom)]
  b <- b[, lapply(.SD, sum, na.rm=TRUE), by=.(ID, Mix), .SDcols=c("Jurkat","IM-9","Raji","THP-1") ] 
  out <- eval_mix(b)
  out <- data.table(out, B = length(nl))
  return(out)
})
p_sims_ts_with_sva_ave <- do.call(rbind, p_sims_ts_with_sva_ave)

fwrite(p_sims_ts_with_sva_ave, file.path(".", "sim_results", "fullTS_v2_ave.txt"))

