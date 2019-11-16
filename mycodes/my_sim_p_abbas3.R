
# Functions + testing -----------------------------------------------------


# Application of a weigthed scheme
# Change the convergence criterion to the mean RNA level instead of the cell type proportion
all_objs <- ls()
all_objs <- all_objs[!all_objs %in% c("mixes_abbas", "cell_type_data_var", "cell_type_data")]
rm(list = all_objs)
row_std_data <- function(x){
  require(Rfast)
  J <- rep(1, ncol(x))
  rs <-  Rfast::rowmeans(x)
  kk <- tcrossprod(rs, J)
  sd_y <- Rfast::rowVars(x, std = T)
  d <- (1/sd_y) * (x - kk)
  return(d)
}

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


my_sim_p_abbas3 <- function(p, use.sva = F, iter = 1e3, tol = 1e-12, eval = T, use.sum = NULL, scaling = T, prop_data = cell_type_data_var, pi_c = true_params$pi_c,
                            all_dt = mixes_abbas, raw_dt = cell_type_data, norm = T, method = "TMM", cov_data = NULL, use.super = F, digits = 4, use.prior = T, output.gene_est = F){
  # require(nnls)
  require(sva)
  require(edgeR)
  require(data.table)
  require(dplyr)
  require(Rfast)
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
    my_params <- raw_dt[info_idx, -c(1), with = F]
    n <- ncol(my_params)
    # Obtain all columns except the last
    mysubset <- lapply(names(my_params)[1:(n-1)], function(x) {c(x, names(my_params)[n])})
    X <- my_params[, lapply(mysubset, function(x) get(x[1])-get(x[2]))]
    X <- data.matrix(cbind(1,X))
    
    # sample_ests_LM2 <- lm(data.matrix(all_dt[info_idx,-1]) ~ 0 + X, offset = unlist(my_params[,n, with = F]))
    # sample_ests_LM2 <- t(sample_ests_LM2$coefficients)
    y_star <- sweep(data.matrix(all_dt[info_idx,-1]), MARGIN = 1, unlist(my_params[,n, with = F]) , `-`)
    sample_ests_LM2 <- Rfast::lmfit(X, y_star)
    sample_ests_LM2 <- t(sample_ests_LM2$be)
    colnames(sample_ests_LM2) <- names(my_params)[1:(n)]
    sample_ests_LM2 <- apply(sample_ests_LM2, 1, function(beta){
      beta <- c(beta[-1], 1 - sum(beta[-1]))
      names(beta) <- names(my_params)[1:(n)]
      if(any(beta < 0)){
        beta[beta < 0] <- 1e-5
        beta <- beta /sum( beta )
      }
      return(beta)
    })
    sample_ests_LM2 <- t(sample_ests_LM2)
    
  }else{
    
    
    
    # Get initial estimates of pi
    my_params <- raw_dt[info_idx, -c(1), with = F]
    n <- ncol(my_params)
    # Obtain all columns except the last
    mysubset <- lapply(names(my_params)[1:(n-1)], function(x) {c(x, names(my_params)[n])})
    X <- my_params[, lapply(mysubset, function(x) get(x[1])-get(x[2]))]
    X <- data.matrix(cbind(1,X))
    
    # sample_ests_LM2 <- lm(data.matrix(all_dt[info_idx,-1]) ~ 0 + X, offset = unlist(my_params[,n, with = F]))
    # sample_ests_LM2 <- t(sample_ests_LM2$coefficients)
    y_star <- sweep(data.matrix(all_dt[info_idx,-1]), MARGIN = 1, unlist(my_params[,n, with = F]) , `-`)
    sample_ests_LM2 <- Rfast::lmfit(X, y_star)
    sample_ests_LM2 <- t(sample_ests_LM2$be)
    
    
    colnames(sample_ests_LM2) <- names(my_params)[1:(n)]
    sample_ests_LM2 <- apply(sample_ests_LM2, 1, function(beta){
      beta <- c(beta[-1], 1 - sum(beta[-1]))
      names(beta) <- names(my_params)[1:(n)]
      if(any(beta < 0)){
        beta[beta < 0] <- 1e-5
        beta <- beta /sum( beta )
      }
      
      # Add noninformative prior?
      if(use.prior){
        beta <- beta + pi_c
        # beta <- beta + seq_len(length(beta))/length(beta)
        beta <- beta/2
        names(beta) <- names(my_params)[1:(n)]
      }
     
      
      
      return(beta)
    })
    sample_ests_LM2 <- t(sample_ests_LM2)
    sample_ests_LM2_old <- sample_ests_LM2
    
    # Scale
    if(scaling){
      w_mixes_abbas <-  row_std_data(data.matrix(all_dt[,-1]))
      # w_mixes_abbas <- t(standardise(t(data.matrix(all_dt[,-1]))))
      
      # w_mixes_abbas <- t(scale(t(data.matrix(all_dt[,-1]))))
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
        
        cnt <- NULL
        if(use.super){
          cnt <- 1 - prop_data$prop
        }
        svseq <- svaseq(analysis_dt2_sva,mod1,mod0, constant = 0, controls = cnt)$sv
        # svseq <- svaseq(analysis_dt2_sva,mod1,mod0, constant = 0)$sv
      }else{
        svseq <- svaseq(data.matrix(all_dt[,-1]),mod1,mod0, constant = 1)$sv
      }
      if( all(svseq == 0) ){
        svseq <- rep(0, nrow(x))
      }
      
    }
    
    if (scaling){
      if(use.sva){
        # gene_ests_LM2 <- lm(data.matrix(t(w_mixes_abbas)) ~ 0 + x + svseq)
        # gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
        X_n <- cbind(x, svseq)
        gene_ests_LM2 <- Rfast::lmfit(X_n, data.matrix(t(w_mixes_abbas)))
        gene_ests_LM2 <- t(gene_ests_LM2$be)[,1:ncol(sample_ests_LM2)]
        gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
        colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
      }else{
        
        # gene_ests_LM2 <- lm(data.matrix(t(w_mixes_abbas)) ~ 0 + x)
        # gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
        gene_ests_LM2 <- Rfast::lmfit(x, data.matrix(t(w_mixes_abbas)))
        gene_ests_LM2 <- t(gene_ests_LM2$be)[,1:ncol(sample_ests_LM2)]
        gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
        colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
        
      }
      
      
    }else{
      if(use.sva){
        # gene_ests_LM2 <- lm(data.matrix(t(all_dt[,-1])) ~ 0 + x + svseq)
        # gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
        X_n <- cbind(x, svseq)
        gene_ests_LM2 <- Rfast::lmfit(X_n, data.matrix(t(all_dt[,-1])))
        gene_ests_LM2 <- t(gene_ests_LM2$be)[,1:ncol(sample_ests_LM2)]
        gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
        gene_ests_LM2 <- apply(gene_ests_LM2, 2, function(b){ b[b <= 0] <- 0; return(b) })
        colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
      }else{
        # gene_ests_LM2 <- lm(data.matrix(t(all_dt[,-1])) ~ 0 + x)
        # gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
        gene_ests_LM2 <- Rfast::lmfit(x, data.matrix(t(all_dt[,-1])))
        gene_ests_LM2 <- t(gene_ests_LM2$be)[,1:ncol(sample_ests_LM2)]
        gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
        colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
        gene_ests_LM2 <- apply(gene_ests_LM2, 2, function(b){ b[b <= 0] <- 0; return(b)})
      }
      
      
    }
    
    
    
    
    
    # Get new estimates of pi
    if(scaling){
      my_params <- data.table(gene_ests_LM2[info_idx, ])
      n <- ncol(my_params)
      # Obtain all columns except the last
      mysubset <- lapply(names(my_params)[1:(n-1)], function(x) {c(x, names(my_params)[n])})
      X <- my_params[, lapply(mysubset, function(x) get(x[1])-get(x[2]))]
      X <- data.matrix(cbind(1,X))
      
      sample_ests_LM2 <- lm(data.matrix(w_mixes_abbas[info_idx, ]) ~ 0 + X, offset = unlist(my_params[,n, with = F]))
      sample_ests_LM2 <- t(sample_ests_LM2$coefficients)
      
      # y_star <- sweep(data.matrix(w_mixes_abbas[info_idx, ]) , MARGIN = 1, unlist(my_params[,n, with = F]) , `-`)
      # sample_ests_LM2 <- Rfast::lmfit(X, y_star)
      # sample_ests_LM2 <- t(sample_ests_LM2$be)
      
      colnames(sample_ests_LM2) <- names(my_params)[1:(n)]
      sample_ests_LM2 <- apply(sample_ests_LM2, 1, function(beta){
        beta <- c(beta[-1], 1 - sum(beta[-1]))
        names(beta) <- names(my_params)[1:(n)]
        if(any(beta < 0)){
          beta[beta < 0] <- 1e-5
          beta <- beta /sum( beta )
        }
        return(beta)
      })
      
    }else{
      
      my_params <- data.table(gene_ests_LM2[info_idx, ])
      n <- ncol(my_params)
      # Obtain all columns except the last
      mysubset <- lapply(names(my_params)[1:(n-1)], function(x) {c(x, names(my_params)[n])})
      X <- my_params[, lapply(mysubset, function(x) get(x[1])-get(x[2]))]
      X <- data.matrix(cbind(1,X))
      
      # sample_ests_LM2 <- lm(data.matrix(all_dt[info_idx, -1]) ~ 0 + X, offset = unlist(my_params[,n, with = F]))
      # sample_ests_LM2 <- t(sample_ests_LM2$coefficients)
      
      y_star <- sweep(data.matrix(all_dt[info_idx, -1]) , MARGIN = 1, unlist(my_params[,n, with = F]) , `-`)
      sample_ests_LM2 <- Rfast::lmfit(X, y_star)
      sample_ests_LM2 <- t(sample_ests_LM2$be)
      
      colnames(sample_ests_LM2) <- names(my_params)[1:(n)]
      sample_ests_LM2 <- apply(sample_ests_LM2, 1, function(beta){
        beta <- c(beta[-1], 1 - sum(beta[-1]))
        names(beta) <- names(my_params)[1:(n)]
        if(any(beta < 0)){
          beta[beta < 0] <- 1e-5
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
            # gene_ests_LM2 <- lm(data.matrix(t(w_mixes_abbas)) ~ 0 + x + svseq)
            # gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
            X_n <- cbind(x, svseq)
            gene_ests_LM2 <- Rfast::lmfit(X_n, data.matrix(t(w_mixes_abbas)))
            gene_ests_LM2 <- t(gene_ests_LM2$be)[,1:ncol(sample_ests_LM2)]
            
            gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
            colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
          }else{
            # gene_ests_LM2 <- lm(data.matrix(t(w_mixes_abbas)) ~ 0 + x)
            # gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
            gene_ests_LM2 <- Rfast::lmfit(x, data.matrix(t(w_mixes_abbas)))
            gene_ests_LM2 <- t(gene_ests_LM2$be)[,1:ncol(sample_ests_LM2)]
            
            gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
            colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
            
          }
          
          
        }else{
          if(use.sva){
            # gene_ests_LM2 <- lm(data.matrix(t(all_dt[,-1])) ~ 0 + x + svseq)
            # gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
            X_n <- cbind(x, svseq)
            gene_ests_LM2 <- Rfast::lmfit(X_n, data.matrix(t(all_dt[,-1])))
            gene_ests_LM2 <- t(gene_ests_LM2$be)[,1:ncol(sample_ests_LM2)]
            
            gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
            gene_ests_LM2 <- apply(gene_ests_LM2, 2, function(b){ b[b <= 0] <- 0; return(b) })
            colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
          }else{
            # gene_ests_LM2 <- lm(data.matrix(t(all_dt[,-1])) ~ 0 + x)
            # gene_ests_LM2 <- t(gene_ests_LM2$coefficients)[,1:ncol(sample_ests_LM2)]
            gene_ests_LM2 <- Rfast::lmfit(x, data.matrix(t(all_dt[,-1])))
            gene_ests_LM2 <- t(gene_ests_LM2$be)[,1:ncol(sample_ests_LM2)]
            
            gene_ests_LM2 <- cbind(gene_ests_LM2[,1], (gene_ests_LM2[, -1] + gene_ests_LM2[,1]))
            colnames(gene_ests_LM2) <- colnames(sample_ests_LM2)
            gene_ests_LM2 <- apply(gene_ests_LM2, 2, function(b){ b[b <= 0] <- 0; return(b)})
          }
          
          
        }
        
        if(scaling){
          my_params <- data.table(gene_ests_LM2[info_idx, ])
          n <- ncol(my_params)
          # Obtain all columns except the last
          mysubset <- lapply(names(my_params)[1:(n-1)], function(x) {c(x, names(my_params)[n])})
          X <- my_params[, lapply(mysubset, function(x) get(x[1])-get(x[2]))]
          X <- data.matrix(cbind(1,X))
          
          sample_ests_LM2 <- lm(data.matrix(w_mixes_abbas[info_idx, ]) ~ 0 + X, offset = unlist(my_params[,n, with = F]))
          sample_ests_LM2 <- t(sample_ests_LM2$coefficients)
          
          # y_star <- sweep(data.matrix(w_mixes_abbas[info_idx, ]), MARGIN = 1, unlist(my_params[,n, with = F]) , `-`)
          # sample_ests_LM2 <- Rfast::lmfit(X, y_star)
          # sample_ests_LM2 <- t(sample_ests_LM2$be)
          
          
          colnames(sample_ests_LM2) <- names(my_params)[1:(n)]
          sample_ests_LM2 <- apply(sample_ests_LM2, 1, function(beta){
            beta <- c(beta[-1], 1 - sum(beta[-1]))
            names(beta) <- names(my_params)[1:(n)]
            if(any(beta < 0)){
              beta[beta < 0] <- 1e-5
              beta <- beta /sum( beta )
            }
            return(beta)
          })
        }else{
          
          my_params <- data.table(gene_ests_LM2[info_idx, ])
          n <- ncol(my_params)
          # Obtain all columns except the last
          mysubset <- lapply(names(my_params)[1:(n-1)], function(x) {c(x, names(my_params)[n])})
          X <- my_params[, lapply(mysubset, function(x) get(x[1])-get(x[2]))]
          X <- data.matrix(cbind(1,X))
          
          # sample_ests_LM2 <- lm(data.matrix(all_dt[info_idx, -1]) ~ 0 + X, offset = unlist(my_params[,n, with = F]))
          # sample_ests_LM2 <- t(sample_ests_LM2$coefficients)
          y_star <- sweep(data.matrix(all_dt[info_idx,-1]), MARGIN = 1, unlist(my_params[,n, with = F]) , `-`)
          sample_ests_LM2 <- Rfast::lmfit(X, y_star)
          sample_ests_LM2 <- t(sample_ests_LM2$be)
          
          colnames(sample_ests_LM2) <- names(my_params)[1:(n)]
          sample_ests_LM2 <- apply(sample_ests_LM2, 1, function(beta){
            beta <- c(beta[-1], 1 - sum(beta[-1]))
            names(beta) <- names(my_params)[1:(n)]
            if(any(beta < 0)){
              beta[beta < 0] <- 1e-5
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
          gene_ests_LM2_old <<- gene_ests_LM2[info_idx, ]
        }
        if(c > iter){
          stop("convergence not achieved, increase tol")
        }
        # print(max(test_est_gene))
      }
      sample_ests_LM2 <<- sample_ests_LM2_out
      gene_ests_LM2_old <<- gene_ests_LM2[info_idx, ]
      
    }
    
    
  }
  
  
  
  sample_ests_LM2 <- apply(sample_ests_LM2, 2, round, digits = digits)
  sample_ests_LM2 <- data.table("ID" = sapply(strsplit(row.names(sample_ests_LM2), " "), "[", 1),
                                "Mix" = sapply(strsplit(row.names(sample_ests_LM2), " "), "[", 2),
                                sample_ests_LM2)
  
  
  
  if(eval){
    # Evaluate
    pear <- eval_mix(sample_ests_LM2)
    pear <- data.table(pear, "p" = p, ngenes = length(info_idx))
    return(pear)
  }else{
    if(!output.gene_est){
      return(sample_ests_LM2)
    }else{
      gene_ests_LM2_out <- data.table(all_dt[info_idx,1],gene_ests_LM2_old)
      return(gene_ests_LM2_out)
    }
    
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
p_sims_ts_with_sva_bagging <- my_sim_p_abbas4(p = .1, use.sva = T, iter = 1e4, tol = 1e-8, perc.sample = 1e-2, r = 250)
p_sims_ts_with_sva_bagging2 <- my_sim_p_abbas4(p = .1, use.sva = T, iter = 1e4, tol = 1e-8, perc.sample = 2e-2, r = 250)
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



# In scRNA data -----------------------------------------------------------

# Functions 
{
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
    
    
    out <- list(out = out, pi_c = unique(pi_c), pi_r = data.table(pi_r), p_v = p_v, mu = mu, out_old = out_old, mean_rna = mean_rna, lambda = lambda)
    return(out)
  }
  
  
  eval_sc <- function(est, truth, digits = 4){
    r <- ncol(est)
    out <- sapply(1:r, function(i){
      e <- est[, i];p <- truth[, i]
      e <- round(e, digits = digits);round(p, digits = digits)
      mean_r <- mean(e/p)
      d <- e-p
      d <- crossprod(d)
      d <- d/length(e)
      rmse <- sqrt(d)
      pear <- cor(e, p)
     return( c(pear, rmse, mean_r))
    })
    out <- t(out)
    colnames(out) <- c("pearson", "rmse", "mean_ratio")
    rownames(out) <- colnames(est)
    return(out)
  }
}

rna_gene_x_cells <- data.table::fread(file.path(".", "cell_data", "rna_gene_x_cells.txt"))
cell_type_data <- data.table::fread(file.path(".", "cell_data", "cell_type_info.txt"))
cell_type_data <- cell_type_data[ all %in% c("Neurons", "OPC",  "Astrocytes",
                                             "Oligodendrocytes", "Microglia", "Endothelial"), ]



# Generate bulk data 
# TODO: Add noise
# After adding noise account for genes with negative expression
# 1. Fix lambda and generate pi from beta for 1e3 subjects
# 2. Do for two cell types
true_params <- est_ref_params(rna_gene_x_cells, cell_type_data)
priors <- cell_type_data[, .N/285, by = Neurons]
sample_size <- 1e2
true_props <- MCMCpack::rdirichlet(sample_size, alpha = round(priors$V1, 2))
true_props <- data.table(true_props)
# true_props$noise <- rnorm(sample_size)
# cor_data <- cor(true_params$out[, -1])
# sim_bulk_data <- rmultinom(sample_size, 1e2, prob =  round(priors$V1, 2))
# t(sim_bulk_data)

rna_gene_x_cells_edited <- rna_gene_x_cells[, c("genes", cell_type_data$cells), with = F]

tt_1 <- apply(rna_gene_x_cells_edited[,-1], 1, function(x){
  cor(as.numeric(cell_type_data$Neurons),x)
})

tt <- apply(true_params$out[, -1], 1, function(x){
  cor(x,c(true_params$pi_c))
})

tt_c <- cor(t(data.matrix(true_params$out[, -1])), t(data.matrix(true_params$pi_c)))
set.seed(49348394)
edited_cell_means <- data.matrix(true_params$out[, c(-1), with = F])
nr <- nrow(edited_cell_means)
perc <- .8
s_nr <- sample(seq_len(nr), size = perc*nr)
edited_cell_means[s_nr,] <- edited_cell_means[s_nr,c(2,1)]
# edited_cell_means[s_nr,] <- t(apply(edited_cell_means[s_nr,], 1, function(x){
#   if(x[1] <= x[2]){
#     x <- x + c(10, 1e2)
#   }else{
#     x <- x + c(1e2, 1e1)
#   }
# }))

# edited_cell_means[s_nr,] <- edited_cell_means[s_nr,] + cbind(rep(10, length(s_nr)), rep(100, length(s_nr)))


tt2 <- apply(edited_cell_means, 1, function(x){
  cor(c(true_params$pi_c),x)
})

sim_bulk_data <- tcrossprod(edited_cell_means, data.matrix(true_props))
sim_bulk_data <- apply(sim_bulk_data, 1, function(x){
  noise <- rnorm(sample_size, sd = sd(x))
  x <- x + abs(noise)
  # x <- ifelse(x < 0 , 0, x)
  return(x)
})
sim_bulk_data <- t(sim_bulk_data)
sim_bulk_data <- data.table(genes = true_params$out$genes , sim_bulk_data)
# 
# sim_bulk_data_cor <- scale(t(sim_bulk_data[, -1]))
# loop <- combn(ncol(sim_bulk_data_cor),2)
# sim_bulk_data_cor <- apply(sim_bulk_data_cor, 1, function(x){
#   out <- 
# })


# drop_gene <- apply(sim_bulk_data, 1, function(x){any(x< 0)})
# sum(drop_gene)

prop_info <- true_params$p_v
names(prop_info)[2] <- "prop"

cell_means <- true_params$out

sc_1 <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 1e4, scaling = T, prop_data = prop_info, 
                        all_dt = sim_bulk_data, raw_dt = cell_means, eval = F, tol = 1e-8)

sc_1_np <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 1e4, scaling = T, prop_data = prop_info, 
                        all_dt = sim_bulk_data, raw_dt = cell_means, eval = F, tol = 1e-8, use.prior = F)

sc_2 <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 1e4, scaling = T, prop_data = prop_info, 
                        all_dt = sim_bulk_data, raw_dt = cell_means, eval = F, tol = 1e-8)



eval_sc(est = data.matrix(sc_1[, -c(1:2), with = F]), truth = data.matrix(true_props))
eval_sc(est = data.matrix(sc_1_np[, -c(1:2), with = F]), truth = data.matrix(true_props))
eval_sc(est = data.matrix(sc_2[, -c(1:2), with = F]), truth = data.matrix(true_props))

# Test MuSiC
{
  packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer", 
             "mixtools","matrixStats", "MuSiC", "convert", "xbioc", "ggplot2", "sva")
  libs_loaded <- sapply(packs, library, character.only = T)
  sc_data <- rna_gene_x_cells_edited[genes %in% sim_bulk_data$genes,]
  sc_data <- sc_data[match(sim_bulk_data$genes, genes), ]
  rn <- sc_data$genes
  cn <- colnames(sc_data)
  cn <- cn[cn %in% cell_type_data$cells]
  sc_data <- sc_data[ , cn, with = F]
  sc_data <- data.matrix(sc_data)
  row.names(sc_data) <- rn
  sc_data <- ExpressionSet(assayData = sc_data)
  bulk_data <- data.matrix(sim_bulk_data[, -1])
  row.names(bulk_data) <- rn
  bulk_data <- ExpressionSet(assayData = bulk_data)
  music_est <- music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters = cell_type_data$Neurons, samples = cell_type_data$cells)
  music_est_nnls <- data.table(samples  = row.names(music_est$Est.prop.allgene), music_est$Est.prop.allgene)
  music_est <- data.table(samples  = row.names(music_est$Est.prop.allgene),music_est$Est.prop.weighted) 
  
  eval_sc(est = data.matrix(music_est_nnls[, -c(1), with = F]), truth = data.matrix(true_props))
  eval_sc(est = data.matrix(music_est[, -c(1), with = F]), truth = data.matrix(true_props))
}

# Check where the estimated genes converge to

{
  # Answer: It depends. If the effect is wide spead then most genes have ranks similar to that of the reference which is problematic. 
  edited_cell_means <- data.matrix(true_params$out[, c(-1), with = F])
  nr <- nrow(edited_cell_means)
  perc <- .65
  s_nr <- sample(seq_len(nr), size = perc*nr)
  edited_cell_means[s_nr,] <- edited_cell_means[s_nr,c(2,1)]
  
  sim_bulk_data <- tcrossprod(edited_cell_means, data.matrix(true_props))
  sim_bulk_data <- apply(sim_bulk_data, 1, function(x){
    noise <- rnorm(sample_size, sd = sd(x))
    x <- x + abs(noise)
    # x <- ifelse(x < 0 , 0, x)
    return(x)
  })
  sim_bulk_data <- t(sim_bulk_data)
  sim_bulk_data <- data.table(genes = true_params$out$genes , sim_bulk_data)
  
  prop_info <- true_params$p_v
  names(prop_info)[2] <- "prop"
  
  cell_means <- true_params$out
  
  sc_1_gene <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 1e4, scaling = T, prop_data = prop_info, 
                               all_dt = sim_bulk_data, raw_dt = cell_means, eval = F, tol = 1e-8, output.gene_est = T)
  
  check_cell_means <- cell_means[genes %in% sc_1_gene$genes,]
  
  check_cell_means_ranks <- sapply( seq_len(nrow(sc_1_gene[,-1])),  function(i){
    identical(rank(sc_1_gene[i,-1]), rank(check_cell_means[i,-1]))
  })
  
  check_cell_means_ranks_original <- sapply( seq_len(nrow(sc_1_gene[,-1])),  function(i){
    identical(rank(sc_1_gene[i,-1]), rank(edited_cell_means[i,]))
  })
  
  mean(!check_cell_means_ranks)
  mean(check_cell_means_ranks_original)
  mean(which(!check_cell_means_ranks) %in% s_nr)
}
  
cor(data.matrix(sc_1[, -c(1:2), with = F]), y = data.matrix(true_props))
tt_est <- cor(t(data.matrix(sc_1[, -c(1:2), with = F])), y = t(data.matrix(edited_cell_means)))
apply(tt_est, 1, function(x) {mean(x > 0)})
mean(tt_c > 0)

apply(tt_est[-82,], 1, function(x) {wilcox.test(x, tt_c, paired = T)$`p.value`})

plot_dt <- data.table(est = unlist(c(sc_1[, -(1:2), with = F])),
                      truth = unlist(c(true_props)),
                      cell_type = rep(c("N-", "N+"), each = sample_size))

ggplot(plot_dt, aes(x = truth, y = est, color = cell_type)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)

# Model Averaged
# With scaling
# end <- .999
# start <- 1e-3
iter <- 1e3



B <- seq(0, 1, length.out = n_ave)
# B <- quantile(true_params$p_v$prop.v_e, probs = B, type = 2)



sc_1_all <- lapply(B[-200], my_sim_p_abbas3, use.sva = F, iter = iter, scaling = T, eval = F, prop_data = prop_info, all_dt = sim_bulk_data, raw_dt = cell_means)
sc_1_all <- lapply(sc_1_all, function(df) {df[, -c(2), with = F]} )
sc_1_all <- do.call(rbind, sc_1_all)
sc_1_all[, p:=rep(B[-200], each = sample_size) + 0.0050251]
sc_1_all[, `:=`(F_1 = `FALSE`*p, `T_1` = `TRUE`*p)]
sc_1_all_eval <- sc_1_all[, .(F_1 = sum(F_1), T_1 = sum(T_1), all_p = sum(p)), by = ID]
sc_1_all_eval[, `:=`( F_1 = F_1/all_p, T_1 = T_1/all_p)]
sc_1_all_eval[, c("ID", "all_p") := NULL]


eval_sc(est = data.matrix(sc_1_all_eval), truth = data.matrix(true_props)) 

plot_dt <- data.table(est = unlist(c(sc_1_all_eval)),
                      truth = unlist(c(true_props)),
                      cell_type = rep(c("N-", "N+"), each = sample_size))

ggplot(plot_dt, aes(x = truth, y = est, color = cell_type)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)


# Simulations 100
sim_code <- function(sample_size = 1e2, prior_dat = priors, cell_mean_data = true_params$out,iter = 1e3, seed, flip,
                     true_params_pv = true_params$p_v, n_ave = 200, diff_means = F, pred = F, use.sva = F, perc = ifelse(diff_means, .7, NULL)){
  require(MCMCpack)
  require(data.table)
  
  # set.seed(49348394)
  if(!missing(seed)){
    set.seed(seed)
  }
  
  true_props <- MCMCpack::rdirichlet(sample_size, alpha = round(prior_dat$V1, 2))
  true_props <- data.table(true_props)
  
 
  
  # true_props$noise <- rnorm(sample_size)
  # cor_data <- cor(true_params$out[, -1])
  # sim_bulk_data <- rmultinom(sample_size, 1e2, prob =  round(priors$V1, 2))
  # t(sim_bulk_data)
  
  
 
  permute_col <- NA
  if(diff_means){
    nr <- nrow(cell_mean_data)
    nc <- ncol(cell_mean_data[, -1])
    
    edited_cell_means <- data.matrix(cell_mean_data[, -1])
    s_nr <- sample(seq_len(nr), size = perc*nr)
    permute_col <- rev(seq_len(nc)) #sample(seq_len(nc))
    if(!missing(flip)){
      if(flip){
        edited_cell_means[s_nr,] <- edited_cell_means[s_nr, permute_col]
        edited_cell_means[s_nr,] <- t(apply(edited_cell_means[s_nr,], 1, function(x){
          x <- x + exp(rank(x)/max(rank(x)))
        }))
      }else{
        edited_cell_means[s_nr,] <- t(apply(edited_cell_means[s_nr,], 1, function(x){
          # if(x[1] <= x[2]){
          #   x <- x + c(10, 1e2)
          # }else{
          #   x <- x + c(1e2, 1e1)
          # }
          x <- x + exp(rank(x)/max(rank(x)))
        }))
        # edited_cell_means[s_nr,] <- edited_cell_means[s_nr,]*exp(fold_change)
        permute_col <- seq_len(nc)
      }
    }else{
      edited_cell_means[s_nr,] <- edited_cell_means[s_nr, permute_col]
    }
    
    
    sim_bulk_data <- tcrossprod(edited_cell_means, data.matrix(true_props))
    
    # scale_pi <- matrix(1, nrow = nr, ncol = nc) %*% diag(rev(1:nc) * scale_factor)
    # scale_pi <- scale_pi * cell_mean_data[, -1]

    
  }else{
    sim_bulk_data <- tcrossprod(data.matrix(cell_mean_data[, -1]), data.matrix(true_props))
  }
  
  
  sim_bulk_data <- apply(sim_bulk_data, 1, function(x){
    noise <- rnorm(sample_size, sd = sd(x))
    x <- x + abs(noise) # <- Half normal noise
    # x <- ifelse(x < 0 , 0, x)
    return(x)
  })
  sim_bulk_data <- t(sim_bulk_data)
  sim_bulk_data <- data.table(genes = cell_mean_data$genes , sim_bulk_data)
  
  # drop_gene <- apply(sim_bulk_data, 1, function(x){any(x< 0)})
  # sum(drop_gene)
  
  prop_info <- true_params_pv
  names(prop_info)[2] <- "prop"
  
  cell_means <- cell_mean_data
  
 
  
  
  B <- seq(0, 1, length.out = n_ave)
  if(n_ave > 1){
    B <- B[-n_ave]
  }
  # B <- quantile(true_params$p_v$prop.v_e, probs = B, type = 2)
  
  
  
  sc_1_all <- lapply(B, my_sim_p_abbas3, use.sva = use.sva, iter = iter, scaling = T, eval = F, prop_data = prop_info, all_dt = sim_bulk_data, raw_dt = cell_means)
  sc_1_all <- lapply(sc_1_all, function(df) {df[, -c(2), with = F]} )
  sc_1_all <- do.call(rbind, sc_1_all)
  
  
  if(n_ave < 1){
    sc_1_all_eval <- sc_1_all
  }else{
    p_dat <- rep(B, each = sample_size) + 1e-6
    p_dat <- p_dat/sum(B + 1e-6)
    sc_1_all[, p:=p_dat]
    n <- ncol(sc_1_all)
    # all_p <- sum(B + 1e-6)
   
    
    mysubset <- lapply(names(sc_1_all)[2:(n-1)], function(x) {c(x, names(sc_1_all)[n])})
    X <- sc_1_all[, lapply(mysubset, function(x) get(x[1])*get(x[2]) )]
    
    
    
    sc_1_all <- data.table(sc_1_all[,1], X)
    
    sc_1_all_eval <- sc_1_all[, lapply(.SD, sum, na.rm=TRUE), by=.(ID), .SDcols=2:(n-1)] 
  }
  
  # p_dat <- p_dat/sum(p_dat)
  
  
  
  
  sc_1_all_eval[, c("ID") := NULL]
  
  # sc_1_all[, `:=`(F_1 = `FALSE`*p, `T_1` = `TRUE`*p)]
  # sc_1_all_eval <- sc_1_all[, .(F_1 = sum(F_1), T_1 = sum(T_1), all_p = sum(p)), by = ID]
  # sc_1_all_eval[, `:=`( F_1 = F_1/all_p, T_1 = T_1/all_p)]
  # sc_1_all_eval[, c("ID", "all_p") := NULL]
  out <- eval_sc(est = data.matrix(sc_1_all_eval), truth = data.matrix(true_props)) 
  if(pred){
    
    ref_alpha <- cor(t(data.matrix(cell_means[, -1])), (data.matrix(prior_dat$V1)))
    
    if(diff_means){
      est_alpha <- cor(t(data.matrix(sc_1_all_eval)), y = t(data.matrix(edited_cell_means)))
    }else{
      est_alpha <- cor(t(data.matrix(sc_1_all_eval)), y = t(data.matrix(cell_mean_data[, -1])))
    }
    
    
    # test_alpha <- apply(est_alpha, 1, function(x) {
    #   p_b <- mean(x > 0)
    #   p_null <- mean(ref_alpha>0)
    #   z <- p_b - p_null
    #   z <- z/sqrt(p_b*(1-p_b) + p_null*(1-p_null))
    #   p_val <- 2*pnorm(z, lower.tail = F)
    #   return(p_val)
    # })
    
    
    # test_alpha <- apply(est_alpha, 1, function(x) {wilcox.test(x, ref_alpha, paired = T)$`p.value`})
    test_alpha <- apply(est_alpha, 1, function(x) {t.test(x, ref_alpha, paired = T)$`p.value`})
    
    sc_1_all_eval$p_val <- test_alpha
    
    
    return(list(est = sc_1_all_eval, truth = true_props))
  }else{
    return(list(eval = out, permuted_cols = permute_col))
  }
  
}

# Scenario 1
# Mean increases by a certain percentage, no flip

sim_1_means_change <- lapply(seq(0,1, by = .1), function(p){
  sim_code(sample_size = 1e2, cell_mean_data = true_params$out, n_ave = 2e2, iter = 1e3, seed = 49348395, diff_means = T, pred = F, flip = F, perc = p)
})

sim_1_means_change <- lapply(sim_1_means_change, "[[", 1)
sim_1_means_change <- do.call(rbind, sim_1_means_change)
sim_1_means_change <- data.table(sim_1_means_change)
sim_1_means_change[, perc := rep(seq(0,1, by = .1), each = 6)]


sim_1_means_change_raw <- lapply(seq(0,1, by = .1), function(p){
  sim_code(sample_size = 1e2, cell_mean_data = true_params$out, n_ave = 2e2, iter = 0, seed = 49348395, diff_means = T, pred = F, flip = F, perc = p)
})

sim_1_means_change_raw <- lapply(sim_1_means_change_raw, "[[", 1)
sim_1_means_change_raw <- do.call(rbind, sim_1_means_change_raw)
sim_1_means_change_raw <- data.table(sim_1_means_change_raw)
sim_1_means_change_raw[, perc := rep(seq(0,1, by = .1), each = 2)]


# Scenario 2 ----> 2 and 3 are the same!
# Means flip only, no increase
sim_2_means_flip <- lapply(seq(0,1, by = .1), function(p){
  sim_code(sample_size = 1e2, cell_mean_data = true_params$out, n_ave = 2e2, iter = 1e3, seed = 49348395, diff_means = T, pred = F, perc = p)
})

sim_2_means_flip_raw <- lapply(seq(0,1, by = .1), function(p){
  sim_code(sample_size = 1e2, cell_mean_data = true_params$out, n_ave = 2e2, iter = 0, seed = 49348395, diff_means = T, pred = F, perc = p)
})

sim_2_means_flip_raw <- lapply(sim_2_means_flip_raw, "[[", 1)
sim_2_means_flip_raw <- do.call(rbind, sim_2_means_flip_raw)
sim_2_means_flip_raw <- data.table(sim_2_means_flip_raw)
sim_2_means_flip_raw[, perc := rep(seq(0,1, by = .1), each = 2)]

# Scenario 3
# Means flip and change
sim_3_means_flip_change <- lapply(seq(0,1, by = .1), function(p){
  sim_code(sample_size = 1e2, cell_mean_data = true_params$out, n_ave = 2e2, iter = 1e3, seed = 49348395, diff_means = T, pred = F, flip = T, perc = p)
})

sim_3_means_flip_change_raw <- lapply(seq(0,1, by = .1), function(p){
  sim_code(sample_size = 1e2, cell_mean_data = true_params$out, n_ave = 2e2, iter = 0, seed = 49348395, diff_means = T, pred = F, flip = T, perc = p)
})

sim_3_means_flip_change_raw <- lapply(sim_3_means_flip_change_raw, "[[", 1)
sim_3_means_flip_change_raw <- do.call(rbind, sim_3_means_flip_change_raw)
sim_3_means_flip_change_raw <- data.table(sim_3_means_flip_change_raw)
sim_3_means_flip_change_raw[, perc := rep(seq(0,1, by = .1), each = 2)]

sim_1 <- replicate(5e2, sim_code(sample_size = 1e1))
sim_2 <- replicate(5e2, sim_code(sample_size = 5e1))
sim_3 <- replicate(5e2, sim_code(sample_size = 1e2))

sim_code(sample_size = 1e1, cell_mean_data = true_params$out, n_ave = 1, iter = 0, seed = 49348395, diff_means = T, pred = F)

# 3. All six cell types
true_params <- est_ref_params(rna_gene_x_cells, cell_type_data[, -2])
priors <- cell_type_data[, .N/285, by = all]


# Sims

{
  
  # Scenario 1
  # Mean increases by a certain percentage, no flip
  
  sim_1_means_sva <- sim_code(sample_size = 1e2, cell_mean_data = true_params$out, n_ave = 2e2, iter = 1e3, seed = 49348395, diff_means = F, pred = F, flip = F, use.sva = T)
  
  
  sim_1_means_change <- lapply(seq(0,1, by = .1), function(p){
    sim_code(sample_size = 1e2, cell_mean_data = true_params$out, n_ave = 2e2, iter = 1e3, seed = 49348395, diff_means = T, pred = F, flip = F, perc = p)
  })
  
  sim_1_means_change <- lapply(sim_1_means_change, "[[", 1)
  sim_1_means_change <- do.call(rbind, sim_1_means_change)
  sim_1_means_change <- data.table(sim_1_means_change)
  sim_1_means_change[, perc := rep(seq(0,1, by = .1), each = 6)]
  sim_1_means_change[, cell_type := rep(priors$all, 11)]
  
  
  sim_1_means_change_raw <- lapply(seq(0,1, by = .1), function(p){
    sim_code(sample_size = 1e2, cell_mean_data = true_params$out, n_ave = 2e2, iter = 0, seed = 49348395, diff_means = T, pred = F, flip = F, perc = p)
  })
  
  sim_1_means_change_raw <- lapply(sim_1_means_change_raw, "[[", 1)
  sim_1_means_change_raw <- do.call(rbind, sim_1_means_change_raw)
  sim_1_means_change_raw <- data.table(sim_1_means_change_raw)
  sim_1_means_change_raw[, perc := rep(seq(0,1, by = .1), each = 6)]
  sim_1_means_change_raw[, cell_type := rep(priors$all, 11)]
  
  
  # Scenario 2 ----> 2 and 3 are the same!
  # Means flip only, no increase
  sim_2_means_flip <- lapply(seq(0,1, by = .1), function(p){
    sim_code(sample_size = 1e2, cell_mean_data = true_params$out, n_ave = 2e2, iter = 1e3, seed = 49348395, diff_means = T, pred = F, perc = p)
  })
  
  sim_2_means_flip <- lapply(sim_2_means_flip, "[[", 1)
  sim_2_means_flip <- do.call(rbind, sim_2_means_flip)
  sim_2_means_flip <- data.table(sim_2_means_flip)
  sim_2_means_flip[, perc := rep(seq(0,1, by = .1), each = 6)]
  sim_2_means_flip[, cell_type := rep(priors$all, 11)]
  
  
  sim_2_means_flip_raw <- lapply(seq(0,1, by = .1), function(p){
    sim_code(sample_size = 1e2, cell_mean_data = true_params$out, n_ave = 2e2, iter = 0, seed = 49348395, diff_means = T, pred = F, perc = p)
  })
  
  sim_2_means_flip_raw <- lapply(sim_2_means_flip_raw, "[[", 1)
  sim_2_means_flip_raw <- do.call(rbind, sim_2_means_flip_raw)
  sim_2_means_flip_raw <- data.table(sim_2_means_flip_raw)
  sim_2_means_flip_raw[, perc := rep(seq(0,1, by = .1), each = 6)]
  sim_2_means_flip_raw[, cell_type := rep(priors$all, 11)]
  
  # Save
  fwrite(sim_1_means_change, file.path(".", "sim_results", "sim_1_means_change_6.txt"))
  fwrite(sim_1_means_change_raw, file.path(".", "sim_results", "sim_1_means_change_raw_6.txt"))
  fwrite(sim_2_means_flip, file.path(".", "sim_results", "sim_2_means_flip_6.txt"))
  fwrite(sim_2_means_flip_raw, file.path(".", "sim_results", "sim_2_means_flip_raw_6.txt"))
  
}

# Check threshold
{
  # Threshold effect
  sample_size <- 1e2
  set.seed(49348395)
  true_props <- MCMCpack::rdirichlet(sample_size, alpha = round(priors$V1, 2))
  true_props <- data.table(true_props)
  edited_cell_means <- data.matrix(true_params$out[, c(-1), with = F])
  sim_bulk_data <- tcrossprod(edited_cell_means, data.matrix(true_props))
  sim_bulk_data <- apply(sim_bulk_data, 1, function(x){
    noise <- rnorm(sample_size, sd = sd(x))
    x <- x + abs(noise)
    # x <- ifelse(x < 0 , 0, x)
    return(x)
  })
  sim_bulk_data <- t(sim_bulk_data)
  sim_bulk_data <- data.table(genes = true_params$out$genes , sim_bulk_data)
  prop_info <- true_params$p_v
  names(prop_info)[2] <- "prop"
  cell_means <- true_params$out
  
  
  sims_1_threshs <- lapply(seq(0,.99, by = 1e-1), function(p){
    print(p)
    out <- my_sim_p_abbas3(p = p, use.sva = F, iter = 1e4, scaling = T, prop_data = prop_info, 
                           all_dt = sim_bulk_data, raw_dt = cell_means, eval = F, tol = 1e-8)
    out <- eval_sc(est = data.matrix(out[, -c(1:2), with = F]), truth = data.matrix(true_props))
    return(out)
  })
  sims_1_threshs <- do.call(rbind, sims_1_threshs)
  sims_1_threshs <- data.table(sims_1_threshs)
  sims_1_threshs[, type := rep(paste0("Type", 1:6),10)]
  fwrite(sims_1_threshs, file.path(".", "sim_results", "sims_1_thresh_6.txt"))
  
  
}

# Plot sims
{
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
  
  library(ggplot2)
  transparent_legend =  theme(
    legend.background = element_rect(fill ="transparent"),
    legend.key = element_rect(fill = "transparent",
                              color = "transparent")
  )
  
  remove_grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(), axis.line = element_line(colour = "black"))
  cols <- RColorBrewer::brewer.pal(length(unique(priors$all)), "Dark2")
  names(cols) <- unique(priors$all)
  lbls <- paste("Type", 1:6)
  names(lbls) <- unique(priors$all)
  shp <- 1:6
  names(shp) <- priors$all
  
  plot_dt <- rbind(sim_1_means_change[, .(pearson, rmse, mean_ratio, perc, cell_type, model = "2-Step")],
                   sim_1_means_change_raw[, .(pearson, rmse, mean_ratio, perc, cell_type, model = "LM")])
  p_1 <- ggplot(data = sim_1_means_change, aes(x = perc, y = pearson, color = cell_type, shape = cell_type)) + 
    geom_point() + 
    geom_smooth(se = F) +
    scale_color_manual(values = cols, labels = lbls) +
    scale_shape_manual(labels = lbls, values = shp) +
    transparent_legend + remove_grid +
    # scale_y_continuous(breaks = seq(.8, 1, by = .05)) + 
    scale_y_continuous(breaks = seq(0, 1, by = .05), limits = c(0,1)) +
    ylab(expression(rho)) + xlab("Fraction of genes affected") +
    ggtitle("1st scenario") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          text = element_text(size = 12),
          axis.title = element_text(face="bold", size = 9),
          axis.text.y=element_text(size = 8, face="bold"),
          axis.text.x=element_text(size = 8, face="bold"),
          legend.position = c(0.5,0.40),
          legend.direction = "horizontal",
          legend.title = element_blank())
  
  p_2 <- ggplot(data = sim_1_means_change, aes(x = perc, y = rmse, color = cell_type, shape = cell_type)) + 
    geom_point() + 
    geom_smooth(se = F) +
    scale_color_manual(values = cols, labels = lbls) +
    scale_shape_manual(labels = lbls, values = shp) +
    transparent_legend + remove_grid +
    # scale_y_continuous(breaks = seq(.8, 1, by = .05)) + 
    scale_y_continuous(breaks = seq(0, .5, by = .05), limits = c(0,.4)) +
    ylab("RMSE") + xlab("Fraction of genes affected") +
    ggtitle("1st scenario") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          text = element_text(size = 12),
          axis.title = element_text(face="bold", size = 9),
          axis.text.y=element_text(size = 8, face="bold"),
          axis.text.x=element_text(size = 8, face="bold"),
          legend.position = "none", #c(0.8,0.20),
          legend.title = element_blank())
  
  p_3 <- ggplot(data = plot_dt, aes(x = perc, y = log(mean_ratio), color = cell_type, shape = model)) + 
    geom_point(size = 2) + 
    geom_smooth(se = F) +
    scale_color_manual(values = cols, labels = lbls) +
    scale_shape_manual(labels = lbls, values = shp) +
    transparent_legend + remove_grid +
    # scale_y_continuous(breaks = seq(.8, 1, by = .05)) + 
    ylab("log(Mean ratio)") + xlab("Fraction of genes affected") +
    ggtitle("Impact of differences in cell-type expression between bulk and reference on estimation \n 2-Step") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          text = element_text(size = 12),
          axis.title = element_text(face="bold", size = 9),
          axis.text.y=element_text(size = 8, face="bold"),
          axis.text.x=element_text(size = 8, face="bold"),
          legend.position = "none",#c(0.8,0.20),
          legend.title = element_blank())
  
  
  plot_dt2 <- rbind(sim_2_means_flip[, .(pearson, rmse, mean_ratio, perc, cell_type, model = "2-Step")],
                    sim_2_means_flip_raw[, .(pearson, rmse, mean_ratio, perc, cell_type, model = "LM")])
  
  o_1 <- ggplot(data = sim_2_means_flip, aes(x = perc, y = pearson, color = cell_type, shape = cell_type)) + 
    geom_point(size = 2) + 
    geom_smooth(se = F) +
    scale_color_manual(values = cols, labels = lbls) +
    scale_shape_manual(labels = lbls, values = shp) +
    transparent_legend + remove_grid +
    scale_y_continuous(breaks = seq(0, 1, by = .05), limits = c(0,1)) +
    ylab(expression(rho)) + xlab("Fraction of genes affected") +
    ggtitle("2nd scenario") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          text = element_text(size = 12),
          axis.title = element_text(face="bold", size = 9),
          axis.text.y=element_text(size = 8, face="bold"),
          axis.text.x=element_text(size = 8, face="bold"),
          legend.position = "none",#c(0.2,0.20),
          legend.title = element_blank())
  
  o_2 <- ggplot(data = sim_2_means_flip, aes(x = perc, y = rmse, color = cell_type, shape = cell_type)) + 
    geom_point(size = 2) + 
    geom_smooth(se = F) +
    scale_color_manual(values = cols, labels = lbls) +
    scale_shape_manual(labels = lbls, values = shp) +
    transparent_legend + remove_grid +
    # scale_y_continuous(breaks = seq(.8, 1, by = .05)) + 
    scale_y_continuous(breaks = seq(0, .5, by = .05), limits = c(0,.4)) +
    ylab("RMSE") + xlab("Fraction of genes affected") +
    ggtitle("2nd scenario") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          text = element_text(size = 12),
          axis.title = element_text(face="bold", size = 9),
          axis.text.y=element_text(size = 8, face="bold"),
          axis.text.x=element_text(size = 8, face="bold"),
          legend.position = "none",#c(0.8,0.20),
          legend.title = element_blank())
  
  o_3 <- ggplot(data = plot_dt2, aes(x = perc, y = log(as.numeric(mean_ratio)), color = cell_type, shape = model)) + 
    geom_point(size = 2) + 
    geom_smooth(se = F) +
    scale_color_manual(values = cols, labels = lbls) +
    scale_shape_manual(labels = lbls, values = shp) +
    transparent_legend + remove_grid +
    # scale_y_continuous(breaks = seq(.8, 1, by = .05)) + 
    ylab("log(Mean ratio)") + xlab("Fraction of genes affected") +
    ggtitle("2nd scenario") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          text = element_text(size = 12),
          axis.title = element_text(face="bold", size = 9),
          axis.text.y=element_text(size = 8, face="bold"),
          axis.text.x=element_text(size = 8, face="bold"),
          legend.position = c(0.8,0.20),
          legend.title = element_blank())
  
  png(file.path(".", "figs", "scenes.png"), res = 200, width = 1920, height = 1080)
  multiplot(p_1, o_1, p_2, o_2, cols = 2)
  dev.off()
}

sim_bulk_data <- tcrossprod(data.matrix(true_params$out[, -1]), data.matrix(true_props))
sim_bulk_data <- apply(sim_bulk_data, 1, function(x){
  noise <- rnorm(sample_size, sd = sd(x))
  x <- x + noise
  x <- ifelse(x < 0 , 0, x)
  return(x)
})
sim_bulk_data <- t(sim_bulk_data)
sim_bulk_data <- data.table(genes = true_params$out$genes , sim_bulk_data)
prop_info <- true_params$p_v
names(prop_info)[2] <- "prop"

cell_means <- true_params$out

sc_1 <- my_sim_p_abbas3(p = .2, use.sva = T, iter = 1e4, scaling = T, prop_data = prop_info, 
                        all_dt = sim_bulk_data, raw_dt = cell_means, eval = F, tol = 1e-6)

eval_sc(est = data.matrix(sc_1[, -c(1:2), with = F]), truth = data.matrix(true_props))

plot_dt <- data.table(est = unlist(c(sc_1[, -(1:2), with = F])),
                      truth = unlist(c(true_props)),
                      cell_type = rep(unique(cell_type_data$all), each = sample_size))

ggplot(plot_dt, aes(x = truth, y = est, color = cell_type)) +
  geom_point(size = 2) + 
  geom_smooth(se = F) +
  geom_abline(slope = 1, intercept = 0)


# 4. Fix pi and total mu, change lambda
priors <- cell_type_data[, .N/285, by = all]
# true_props <- MCMCpack::rdirichlet(sample_size, alpha = round(priors$V1, 2))
true_props <- tcrossprod(rep(1,sample_size),priors$V1)

s <- cov(data.matrix(true_params$out[, -1]))
# diag(s) <- rep(1e-3, ncol(s))
lamd <- mvtnorm::rmvnorm(sample_size, mean = seq(10, 60, by = 10), sigma = s)



sim_bulk_data <- tcrossprod(data.matrix(true_params$out[, -1]), data.matrix(true_props))
sim_bulk_data <- apply(sim_bulk_data, 1, function(x){
  noise <- rnorm(sample_size, sd = sd(x))
  x <- x + noise
  x <- ifelse(x < 0 , 0, x)
  return(x)
})
sim_bulk_data <- t(sim_bulk_data)
sim_bulk_data <- data.table(genes = true_params$out$genes , sim_bulk_data)
prop_info <- true_params$p_v
names(prop_info)[2] <- "prop"

cell_means <- true_params$out

sc_1 <- my_sim_p_abbas3(p = .2, use.sva = T, iter = 1e4, scaling = T, prop_data = prop_info, 
                        all_dt = sim_bulk_data, raw_dt = cell_means, eval = F, tol = 1e-6)

eval_sc(est = data.matrix(sc_1[, -c(1:2), with = F]), truth = data.matrix(true_props))

plot_dt <- data.table(est = unlist(c(sc_1[, -(1:2), with = F])),
                      truth = unlist(c(true_props)),
                      cell_type = rep(unique(cell_type_data$all), each = sample_size))

ggplot(plot_dt, aes(x = truth, y = est, color = cell_type)) +
  geom_point(size = 2) + 
  geom_smooth(se = F) +
  geom_abline(slope = 1, intercept = 0)


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
data.table::fwrite(data.table::data.table(results = my_s_k, cell_types = types_cell), "music_cell_size.txt")

# Toy plots ---------------------------------------------------------------

ps <- cell_type_data[, .N/285, by = all]

toy_cell <- rmvnorm(1e3, mu = seq(10,60,by = 10), sigma = diag(1e-2,6,6))

# samps <- round(ps$V1, 3)*200
k <- 3
samps <- rep(5,k)
mu_s <- seq(10,k*10,by = 10)
mu_s <- mu_s/mean(mu_s)
s <- c(2,.7,.8,1)

toy_data <- sapply(1:1e2, function(x){
  out <- sapply(1:k, function(i){
    rnorm(samps[i],mean = mu_s[i], sd = s[i])
  })
  # out <- do.call(c, out)
  # out <- t(out)
  return(out)
})
toy_data <- t(toy_data)

s <- s/sum(s)
s[k] <- 1e-2


toy_data2 <- sapply(1:1e2, function(x){
  out <- sapply(1:k, function(i){
    rnorm(samps[i],mean = mu_s[i], sd = .2)
  })
  # out <- do.call(c, out)
  # out <- t(out)
  return(out)
})
toy_data2 <- t(toy_data2)
# toy_data <- rbind(toy_data[1:100,], toy_data2[1:50,], toy_data[1:10,], toy_data2[1:20,], toy_data[1:90,],  toy_data2[1:30,])
# toy_data <- rbind(toy_data, toy_data2)
# toy_data <- toy_data[-c(100:110),]


# colnames(toy_data) <- rep(paste("type", seq_len(4)), each = 10)

heatmap(toy_data, Colv = NA, Rowv = NA, col =  RColorBrewer::brewer.pal(k, "Spectral"), labRow = NA, labCol = NA)

heatmap(toy_data2, Colv = NA, Rowv = NA, col =  RColorBrewer::brewer.pal(k, "Spectral"), labRow = NA, labCol = NA)
cc <- RColorBrewer::brewer.pal(k, "Spectral")

heatmap(toy_data2, Colv = NA, Rowv = NA, col =  c("#E2CE92","#B5D293","#C7D093"), labRow = NA, labCol = NA)

ggplot(data.table("type" = rep(letters[1:4], 2),
                  vals = c(10,40,20,30, 30,20,40,10),
                  "pop" = rep(paste0("Population ", 1:2), each = 4)),
       aes(x = type, y= vals, fill = type)) +
  geom_col(width = .4) +
  facet_grid(pop~.) +
  ylab("") + xlab("") +
  transparent_legend + remove_grid +
  theme(aspect.ratio = 1/2,
        legend.position = "none",
        legend.title = element_blank(),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13))
ggsave(file.path(".", "figs","toy_pops.png"), dpi = "retina")

ggplot(data.table("type" = rep(letters[1:4], 2),
                  vals = c(10,40,20,30, 30,100,40,50),
                  "pop" = rep(paste0("Population ", 1:2), each = 4)),
       aes(x = type, y= vals, fill = type)) +
  geom_col(width = .4) +
  facet_grid(pop~.) +
  ylab("") + xlab("") +
  transparent_legend + remove_grid +
  theme(aspect.ratio = 1/2,
        legend.position = "none",
        legend.title = element_blank(),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13))
ggsave(file.path(".", "figs","toy_pops_2.png"), dpi = "retina")


pp <- c(0.1, .4, .5)

g_m <- 100

ll <- c(10, 6, 2)
ll2 <- rev(ll)

sum(pp * ll * g_m)
sum(rev(pp) * ll2 * g_m)

ll2 * g_m
ll * g_m
