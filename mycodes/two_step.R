# Preamble ----------------------------------------------------------------

# Load libraries
packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "MuMIn",
           "RColorBrewer",  "mixtools","matrixStats", "sva", "ggplot2", "sva", "nlme")
loaded_libs <- sapply(packs, library, character.only = T, warn.conflicts = F, quietly = F)

mixes_abbas <- fread(file.path(".", "cell_data", "abbas_et_al_data", "GSE11103_matrix_mixtures.txt"))
cell_type_data <- fread(file.path(".", "cell_data", "abbas_et_al_data", "GSE11103_matrix_pure.txt"))
nc <- ncol(cell_type_data)
first <- seq(1,12, by = 3)
second <- first + 1
third <- first + 2

# account for first col
first <- first + 1
second <- second + 1
third <- third + 1

first <- cell_type_data[, c(1, first), with = F]
second <- cell_type_data[, c(1, second), with = F]
third <- cell_type_data[, c(1, third), with = F]

all <- seq(1,12, by = 3) + 1
mysubset <- lapply(all, function(x) {
  names(cell_type_data)[x: (x+2)]
})

# calculate E(Var(Y|X))

var_all <- apply(cell_type_data[, -1], 1, var)
var_all <- cbind(cell_type_data[, 1],var_all)


pcs <- prcomp(t(cell_type_data[, -1]), scale. = T)
screeplot(pcs)

x <- rep(c("Jurkat", "IM-9", "Raji", "THP-1"), each = 3)
pc_beta <- lm(pcs$x ~ x)
pc_beta <- summary(pc_beta)
p_vals <- sapply(pc_beta, function(x){
  pf(q = x$fstatistic[1], df1 = x$fstatistic[2], df2 =x$fstatistic[3], lower.tail = F)
})

idx_pc <- which(p.adjust(p_vals , method = "bonferroni")<= 5e-2)
x <- pcs$x[, idx_pc]


pc_beta <- lm(t(cell_type_data[, -1]) ~ x)
pc_beta <- summary(pc_beta)
p_vals <- sapply(pc_beta, function(x){
  pf(q = x$fstatistic[1], df1 = x$fstatistic[2], df2 =x$fstatistic[3], lower.tail = F)
})
rm(x)

cell_type_data_var <-  cbind(cell_type_data[, 1],
                             apply(cbind(apply(cell_type_data[, all[1]:(all[1]+2), with = F], 1, var),
                                         apply(cell_type_data[, all[2]:(all[2]+2), with = F], 1, var),
                                         apply(cell_type_data[, all[3]:(all[3]+2), with = F], 1, var),
                                         apply(cell_type_data[, all[4]:(all[4]+2), with = F], 1, var)), 1, mean))


names(cell_type_data_var)[2] <- "E(Var(Y|X))"


cell_type_data_test_file <- file.path(".", "cell_data", "cell_type_data_test.txt")
if(file.exists(cell_type_data_test_file)){
  cell_type_data_test <- fread(cell_type_data_test_file)
}else{
  c <- 0
  cell_type_data_test <- apply(cell_type_data[, -1], 1, function(y){
    
    x <- rep(c("Jurkat", "IM-9", "Raji", "THP-1"), each = 3)
    dt <- data.table(y, x, t = rep(1:3, 4))
    # ft <- gls(y ~ x,
    #             correlation = corSymm(form = ~ t | x),
    #             data=dt,
    #             method="REML")
    c <<- c+1
    # y <<- y
    ft <- lme(y ~ x, 
              random = ~1|x,
              correlation = corSymm(form = ~ t | x),
              data=dt,
              method="REML",
              control = lmeControl(opt = "optim"))
    r.sqrd <- r.squaredGLMM(ft) %>% c
    ft <- car::Anova(ft)
    out <- ft[3]
    print(c/nrow(cell_type_data))
    out <- c(out, r.sqrd)
    # ft <- kruskal.test(y, g = as.factor(x))
    # out <- ft$p.value
    return(out)
  })
  cell_type_data_test <- do.call(rbind, cell_type_data_test)
  cell_type_data_test <- cbind(cell_type_data[, 1], cell_type_data_test)
  
  cell_type_data_test[, p.adjst:= p.adjust(`Pr(>Chisq)`, method = "BH")]
  fwrite(cell_type_data_test, file = cell_type_data_test_file)
}

summary(unlist(cell_type_data_test$`Pr(>Chisq)`))




cell_type_data <- cbind(cell_type_data[, 1],
                        apply(cell_type_data[, all[1]:(all[1]+2), with = F], 1, mean),
                        apply(cell_type_data[, all[2]:(all[2]+2), with = F], 1, mean),
                        apply(cell_type_data[, all[3]:(all[3]+2), with = F], 1, mean),
                        apply(cell_type_data[, all[4]:(all[4]+2), with = F], 1, mean))
names(cell_type_data) <- c("!Sample_title", "Jurkat", "IM-9", "Raji", "THP-1")

# calculate Var(E(Y|X)) 
cell_type_data_var <-  cbind(cell_type_data_var,
                             apply(cell_type_data[, -c(1), with = F], 1, var))
names(cell_type_data_var) <-  c("!Sample_title", "E(Var(Y|X))", "Var(E(Y|X))")
cell_type_data_var[,  total_var:= `E(Var(Y|X))` + `Var(E(Y|X))`  ]
cell_type_data_var[, prop := ifelse(`E(Var(Y|X))` > 0, `Var(E(Y|X))`/total_var, 0)]
cell_type_data_var[, `:=`(f_test=cell_type_data_test$`Pr(>Chisq)`,
                          r2 = unlist(cell_type_data_test$R2m),
                          p_vals = p_vals)]


with(cell_type_data_var, plot(prop, (f_test), main = "Relationship between p values from DE and p", ylab = "F-test unadjst p vals", xlab = "p"))
with(cell_type_data_var, plot(prop, (r2), main = "Relationship between r2 and p", ylab = "r2", xlab = "p"))
with(cell_type_data_var, plot(r2, (f_test), main = "Relationship between p values from DE and p", ylab = "F-test unadjst p vals", xlab = "p"))
with(cell_type_data_var, plot(r2, (p_vals), main = "Relationship between r2 and p_vals", ylab = "p_vals", xlab = "p"))
with(cell_type_data_var, plot(p_vals, (f_test), main = "Relationship between p values from DE and p_vals", ylab = "F-test unadjst p vals", xlab = "p_vals"))

with(cell_type_data_var[prop > .95,], plot(prop, (r2), main = "Relationship between r2 and p", ylab = "r2", xlab = "p"))

summary(cell_type_data_var$r2)
summary(cell_type_data_var$prop)
cell_type_data_var[, prop := r2]
cell_type_data_var[prop < .2,]

# Analyses ----------------------------------------------------------------

cell_type_data_lambda <-  cbind(cell_type_data[, 1],
                                1/apply(cell_type_data[, -1], 1, mean))

cell_type_data_mu <- cell_type_data
cell_type_data <- inner_join(cell_type_data, cell_type_data_lambda) %>% data.table
cell_type_data_mat <- cell_type_data[, -1]
n <- ncol(cell_type_data_mat)
mysubset <- lapply(names(cell_type_data_mat)[1:(n-1)], function(x) {c(x, names(cell_type_data_mat)[n])})
cell_type_data_mat <- cell_type_data_mat[, lapply(mysubset, function(x) (get(x[1])) * (get(x[2])) )]
names(cell_type_data_mat) <- c("Jurkat", "IM-9", "Raji", "THP-1")

# prep mixes data
analysis_dt <- inner_join(mixes_abbas, cell_type_data_lambda) %>% data.table
analysis_dt <- analysis_dt[, -1]
n <- ncol(analysis_dt)
mysubset <- lapply(names(analysis_dt)[1:(n-1)], function(x) {c(x, names(analysis_dt)[n])})
analysis_dt <- analysis_dt[, lapply(mysubset, function(x) (get(x[1])) * (get(x[2])) )]
names(analysis_dt) <- names(mixes_abbas)[-1]

# truth
MixA <- c("Jurkat" = 2.5, "IM-9" = 1.25, "Raji" =  2.5,"THP-1" = 3.75)/10
MixB <-  c("Jurkat" = 0.5, "IM-9" = 3.17, "Raji" = 4.75, "THP-1" = 1.58)/10
MixC <-  c("Jurkat" = 0.1, "IM-9" = 4.95, "Raji" = 1.65, "THP-1" = 3.3)/10
MixD  <- c("Jurkat" = 0.02, "IM-9" = 3.33, "Raji" = 3.33, "THP-1" = 3.33)/10


# subest to informative genes
p <- .95

Mix <- rbind(MixA, MixB, MixC, MixD)
Mix <- data.table("Mix" = paste0("Mix", LETTERS[1:4]), Mix)

my_sim_p_abbas2 <- function(p, use.sva = F, iter = 1, tol = 1e-6, eval = T, use.sum = NULL, scaling = T, prop_data = cell_type_data_var, dt = analysis_dt, lm2_dt = cell_type_data_mat, all_dt = mixes_abbas){
  # require(nnls)
  require(sva)
  # Subset to informative genes from reference
  
  if(!is.null(use.sum)){
    if (use.sum == "pval") {
      use.pval = T; use.sva = T
    }else if (use.sum == "r2"){
      use.r2 = T; use.sva = T
    }
    
  }else{
    use.pval = F; use.r2 = F
  }
  
  
  # Estimate
  cutoff <- quantile(prop_data$prop, probs = c(1-p,p))
  cutoff[2] <- ifelse(p == 0, 0, cutoff[2])
  info_idx <- which(prop_data$prop >= cutoff[2] )
  
  # Get initial estimates of pi
  sample_ests_LM2 <- apply(dt, 2, function(y){
    
    # LM
    y <- y[info_idx]
    my_params <- lm2_dt[info_idx, ]
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
    
    
    beta <- solve(crossprod(x)) %*% crossprod(x, y_star)
    beta <- c(beta)
    # if(scaling){
    #   beta <- c(beta, 1 - sum(beta))
    #   names(beta) <- names(my_vals)[2:(n)]
    # }
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
    # w <- all_dt[,-1]
    # w  <- var(w)
    # w  <- solve(w )
    # w  <- chol(w )
    # # account for whitening in SVA step
    # w_mixes_abbas <- tcrossprod(data.matrix(all_dt[,-1]), w) 
    w_mixes_abbas <- t(scale(t(data.matrix(all_dt[,-1]))))
  }
  
  
  
  # SSVA to eliminate noise
  x <- sample_ests_LM2[, -1]
  x <- cbind(1, x)
  if(use.sva){
    mod1 <-  x
    mod0 <- cbind(mod1[,1])
    if (scaling){
      analysis_dt2_sva <- apply(w_mixes_abbas, 2, function(x) {exp(x) })
      svseq <- svaseq(analysis_dt2_sva,mod1,mod0, constant = 0)$sv
    }else{
      svseq <- svaseq(data.matrix(all_dt[,-1]),mod1,mod0, constant = 1)$sv
    }
    
  }
  
  if (scaling){
    gene_ests_LM2 <- apply(w_mixes_abbas, 1, function(y){
      
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
      
      # account for betas < 0
      # If any beta < 0 set to zero or a very small value
      if(any(b <= 0) & !scaling){
        if(use.sva){
          b[b <= 0] <- 0
        }else{
          # ft <- nnls(x,y)
          # beta <- coef(ft)
          # beta <- beta[1:n]
          # beta_1 <- sapply(beta[-1], function(b) {b <- b + beta[1]; return(b)} )
          # b <- c(beta[1], beta_1)
          # names(b) <- colnames(sample_ests_LM2)
          b[b <= 0] <- 0
        }
        
      }
      
      
      
      if(use.r2 & use.sva){
        ft.anova <- summary(ft)$adj.r.squared
        b <- c(b, "r2" = ft.anova)
      }else if(use.pval & use.sva){
        ft.anova <-  anova(ft)
        ft.anova <- ft.anova$`Pr(>F)`[1]
        b <- c(b, "pval" = ft.anova)
      }
      
      
      return(b)
      
    })
  }else{
    gene_ests_LM2 <- apply(all_dt[,-1], 1, function(y){
      
      # LM
      n <- ncol(sample_ests_LM2)
      
      
      if(use.sva){
        ft <- lm(y ~ -1 + x + svseq)
      }else{
        
        ft <- lm(y ~ -1 + x )
        
        # wts <- 1/fitted(lm(abs(residuals(ft)) ~ -1 + x))^2
        # model.1 <- lm(y ~ -1 + x, weights=wts)
        # wts <- 1/fitted(lm(abs(residuals(model.1)) ~ -1 + x))^2
        # model.2 <- lm(y ~ -1 + x, weights=wts)
        # wts <- 1/fitted(lm(abs(residuals(model.2)) ~ -1 + x))^2
        # 
        # err <- coef(model.2) - coef(model.1)
        # err <- crossprod(err)/length(err)
        # while (err > tol) {
        #   model.1 <- lm(y ~ -1 + x, weights=wts)
        #   wts <- 1/fitted(lm(abs(residuals(model.1)) ~ -1 + x))^2
        #   
        #   model.2 <- lm(y ~ -1 + x, weights=wts)
        #   
        #   err <- coef(model.2) - coef(model.1)
        #   err <- crossprod(err)/length(err)
        #   
        #   if(err < tol){
        #     ft <<- model.2
        #   }
        #   
        # }
        
      }
      
      
      
      
      
      beta <- coef(ft)
      beta <- beta[1:n]
      beta_1 <- sapply(beta[-1], function(b) {b <- b + beta[1]; return(b)} )
      b <- c(beta[1], beta_1)
      names(b) <- colnames(sample_ests_LM2)
      
      # account for betas < 0
      # If any beta < 0 set to zero or a very small value
      if(any(b <= 0) & !scaling){
        if(use.sva){
          b[b <= 0] <- 0
        }else{
          # ft <- nnls(x,y)
          # beta <- coef(ft)
          # beta <- beta[1:n]
          # beta_1 <- sapply(beta[-1], function(b) {b <- b + beta[1]; return(b)} )
          # b <- c(beta[1], beta_1)
          # names(b) <- colnames(sample_ests_LM2)
          b[b <= 0] <- 0
        }
        
      }
      
      
      
      if(use.r2 & use.sva){
        ft.anova <- summary(ft)$adj.r.squared
        b <- c(b, "r2" = ft.anova)
      }else if(use.pval & use.sva){
        ft.anova <-  anova(ft)
        ft.anova <- ft.anova$`Pr(>F)`[1]
        b <- c(b, "pval" = ft.anova)
      }
      
      
      return(b)
      
    })
  }
  
  
  
  gene_ests_LM2 <- t(gene_ests_LM2)
  if(use.pval){
    ncol_gene <- ncol(gene_ests_LM2)
    pval_gene <- gene_ests_LM2[, ncol_gene]
    qval_gene <- p.adjust(pval_gene, method = "bon")
    qval_gene <- which(qval_gene <= 5e-2)
    info_idx1 <- intersect(info_idx, qval_gene)
    if(length(info_idx1) > 0){
      info_idx <- info_idx1
      # stop("Not enough genes")
    }
    gene_ests_LM2 <- gene_ests_LM2[, -ncol_gene]
  }else if(use.r2){
    ncol_gene <- ncol(gene_ests_LM2)
    r2_gene <- gene_ests_LM2[, ncol_gene]
    r2_gene <- which(qval_gene >= p)
    info_idx1 <- intersect(info_idx, qval_gene)
    if(length(info_idx) > 0){
      info_idx <- info_idx1
      # stop("Not enough genes")
    }
    gene_ests_LM2 <- gene_ests_LM2[, -ncol_gene]
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
      
      
      # if(scaling){
      #   y_star <- scale(y_star)
      #   x <- scale(x[, -1])
      # }
      
      
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
  
  # Iterate until convergence
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
        if (scaling){
          analysis_dt2_sva <- apply(w_mixes_abbas, 2, function(x) {exp(x) })
          svseq <- svaseq(analysis_dt2_sva,mod1,mod0, constant = 1e-5)$sv
        }else{
          svseq <- svaseq(data.matrix(all_dt[,-1]),mod1,mod0, constant = 1)$sv
        }
        
      }
      
      if(scaling){
        gene_ests_LM2 <- apply(w_mixes_abbas, 1, function(y){
          
          
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
          if(any(b <= 0) & !scaling){
            if(use.sva ){
              b[b <= 0] <- 1e-3#0.1
            }else{
              # ft <- nnls(x,y)
              # beta <- coef(ft)
              # beta <- beta[1:n]
              # beta_1 <- sapply(beta[-1], function(b) {b <- b + beta[1]; return(b)} )
              # b <- c(beta[1], beta_1)
              # names(b) <- colnames(sample_ests_LM2)
              b[b <= 0] <- 1e-3#0.1
            }
            
          }
          return(b)
          
        })
      }else{
        gene_ests_LM2 <- apply(all_dt[,-1], 1, function(y){
          
          
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
          if(any(b <= 0) & !scaling){
            if(use.sva ){
              b[b <= 0] <- 1e-3#0.1
            }else{
              # ft <- nnls(x,y)
              # beta <- coef(ft)
              # beta <- beta[1:n]
              # beta_1 <- sapply(beta[-1], function(b) {b <- b + beta[1]; return(b)} )
              # b <- c(beta[1], beta_1)
              # names(b) <- colnames(sample_ests_LM2)
              b[b <= 0] <- 1e-3#0.1
            }
            
          }
          return(b)
          
        })
      }
     
      gene_ests_LM2 <- t(gene_ests_LM2)
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
          
          
          # if(scaling){
          #   y_star <- scale(y_star)
          #   x <- scale(x[, -1])
          # }
          
          
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
  
  sample_ests_LM2 <- apply(sample_ests_LM2, 2, round, 4)
  sample_ests_LM2 <- data.table("ID" = sapply(strsplit(row.names(sample_ests_LM2), " "), "[", 1),
                                "Mix" = sapply(strsplit(row.names(sample_ests_LM2), " "), "[", 2),
                                sample_ests_LM2)
  
  
  
  if(eval){
    
    # truth
    MixA <- c("Jurkat" = 2.5, "IM-9" = 1.25, "Raji" =  2.5,"THP-1" = 3.75)/10
    MixB <-  c("Jurkat" = 0.5, "IM-9" = 3.17, "Raji" = 4.75, "THP-1" = 1.58)/10
    MixC <-  c("Jurkat" = 0.1, "IM-9" = 4.95, "Raji" = 1.65, "THP-1" = 3.3)/10
    MixD  <- c("Jurkat" = 0.02, "IM-9" = 3.33, "Raji" = 3.33, "THP-1" = 3.33)/10
    
    
    
    Mix <- rbind(MixA, MixB, MixC, MixD)
    Mix <- data.table("Mix" = paste0("Mix", LETTERS[1:4]), Mix)
    
    
    # Evaluate
    
    pear <- sapply(1:nrow(sample_ests_LM2), function(i){
      dt <- sample_ests_LM2[i,]
      m <- dt$Mix
      dt <- dt[, -c(1:2), with = F]
      dt <- unlist(dt)
      t <- Mix[Mix == m, -1, with = F]
      t <- unlist(t)
      r <- cor(dt, t)
      return(r)
    })
    
    pear_rmse <- sapply(1:nrow(sample_ests_LM2), function(i){
      dt <- sample_ests_LM2[i,]
      m <- dt$Mix
      dt <- dt[, -c(1:2), with = F]
      dt <- unlist(dt)
      t <- Mix[Mix == m, -1, with = F]
      t <- unlist(t)
      r <- crossprod(dt, t)/length(t)
      r <- c(r)
      return(r)
    })
    
    pear <- tapply(pear, rep(paste0("Mix", LETTERS[1:4]), each = 3), mean)
    pear_rmse <- tapply(pear_rmse, rep(paste0("Mix", LETTERS[1:4]), each = 3), mean)
    pear <- data.table(type = names(pear), pearson = (pear), rmse = pear_rmse,"p" = p, ngenes = length(info_idx))
    return(pear)
  }else{
    return(sample_ests_LM2)
  }
  
  
}



old.mem <- memory.limit()
memory.limit(40e3)

# With scaling
end <- .999
iter <- 1
jump <- 1e-3

my_sim_p_abbas2(p, iter = iter, scaling = F)

p_sims_ts <- lapply(seq(0, end, by = jump), my_sim_p_abbas2, use.sva = F, iter = iter, scaling = T)
p_sims_ts <- do.call(rbind, p_sims_ts)
p_sims_ts[, `:=`("scale" = T, "SVA" = F)]

p_sims_ts_with_sva <- lapply(seq(0,end, by = jump), my_sim_p_abbas2, use.sva = T, iter = iter, scaling = T)
p_sims_ts_with_sva <- do.call(rbind, p_sims_ts_with_sva)
p_sims_ts_with_sva[, `:=`("scale" = T, "SVA" = T)]

# No scaling
p_sims_ts_no_scale_no_sva <- lapply(seq(0,end, by = jump), my_sim_p_abbas2, use.sva = F, iter = iter, scaling = F)
p_sims_ts_no_scale_no_sva <- do.call(rbind, p_sims_ts_no_scale_no_sva)
p_sims_ts_no_scale_no_sva[, `:=`("scale" = F, "SVA" = F)]

p_sims_ts_no_scale <- lapply(seq(0,end, by = jump), my_sim_p_abbas2, use.sva = T, iter = iter, scaling = F)
p_sims_ts_no_scale <- do.call(rbind, p_sims_ts_no_scale)
p_sims_ts_no_scale[, `:=`("scale" = F, "SVA" = T)]

p_sims_all <- rbind(p_sims_ts, p_sims_ts_with_sva, p_sims_ts_no_scale_no_sva, p_sims_ts_no_scale)

fwrite(p_sims_all, file.path(".", "sim_results", "p_sims_TS_with_filter_abbas_all.txt"))
fwrite(p_sims_ts, file.path(".", "sim_results", "p_sims_fullTS_with_filter_abbas.txt"))


p_sims_ts <- fread(file.path(".", "sim_results", "p_sims_fullTS_with_filter_abbas.txt"))

memory.limit(old.mem)

{
  # 
  # p_sims_ts <- lapply(seq(0,.99, by = 1e-2), my_sim_p_abbas2, use.sva = F, iter = 1, scaling = F)
  # p_sims_ts <- do.call(rbind, p_sims_ts)
  # fwrite(p_sims_ts, file.path(".", "sim_results", "p_sims_TS_with_filter_old_abbas.txt"))
  # 
  # p_sims_ts <- lapply(seq(0,.99, by = 1e-2), my_sim_p_abbas2, use.sva = F, iter = 1)
  # p_sims_ts <- do.call(rbind, p_sims_ts)
  # fwrite(p_sims_ts, file.path(".", "sim_results", "p_sims_TS_with_filter_abbas.txt"))
  # 
}

# Plot --------------------------------------------------------------------
p_sims_ts_old <- fread(file.path(".", "sim_results", "p_sims_TS_abbas.txt"))
p_sims_ts_full <- fread(file.path(".", "sim_results", "p_sims_fullTS_with_filter_abbas.txt"))
p_sims_ts <- fread(file.path(".", "sim_results", "p_sims_TS_with_filter_abbas.txt"))
transparent_legend =  theme(
  legend.background = element_rect(fill ="transparent"),
  legend.key = element_rect(fill = "transparent",
                            color = "transparent")
)

remove_grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))

no_x_axis_label <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

cols <- RColorBrewer::brewer.pal(length(unique(p_sims_ts$type)), "Dark2")
names(cols) <- unique(p_sims_ts$type)

ggplot(data = p_sims_ts, aes(x = p, y = pearson, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  scale_y_continuous(breaks = seq(.75, 1, by = .05)) + 
  ylab("Pearson correlation coefficient") + xlab("CDF") +
  ggtitle("Relationship between cutoff and model fit\n 1 Iter") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.9,0.20),
        legend.title = element_blank())
ggsave(file.path(".", "figs","effect_of_cutoff_2step_1iter.png"), dpi = "retina")

ggplot(data = p_sims_ts_full, aes(x = p, y = pearson, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  scale_y_continuous(breaks = seq(.75, 1, by = .05)) + 
  ylab("Pearson correlation coefficient") + xlab("CDF") +
  ggtitle("Relationship between cutoff and model fit \n Full") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.9,0.20),
        legend.title = element_blank())
ggsave(file.path(".", "figs","effect_of_cutoff_2step_full.png"), dpi = "retina")


ggplot(data = p_sims_ts_old, aes(x = p, y = pearson, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  scale_y_continuous(breaks = seq(.65, 1, by = .05)) + 
  ylab("Pearson correlation coefficient") + xlab("CDF") +
  ggtitle("Relationship between cutoff and model fit\n 1 Iter") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.9,0.20),
        legend.title = element_blank())
ggsave(file.path(".", "figs","effect_of_cutoff_2step_1iter.png"), dpi = "retina")


sample_ests_ts_info <- my_sim_p_abbas2(0, iter = 2, scaling = T, use.sva = T, eval = F)

pear <- sapply(1:nrow(sample_ests_ts_info), function(i){
  dt <- sample_ests_ts_info[i,]
  m <- dt$Mix
  dt <- dt[, -c(1:2), with = F]
  dt <- unlist(dt)
  t <- Mix[Mix == m, -1, with = F]
  t <- unlist(t)
  r <- cor(dt, t)
  return(r)
})
pear <- tapply(pear, rep(paste0("Mix", LETTERS[1:4]), each = 3), mean)
Mix <- inner_join(sample_ests_ts_info[, .(Mix)],Mix) %>% data.table
plot_dt_ts <- rbind(cbind(sample_ests_ts_info[, .(Mix, Estimate = Jurkat)],  Mix[, .(Truth = Jurkat, type = "Jurkat")]),
                    cbind(sample_ests_ts_info[, .(Mix, Estimate = `IM-9`)],  Mix[, .(Truth = `IM-9`, type = "IM-9")]),
                    cbind(sample_ests_ts_info[, .(Mix, Estimate = Raji)],  Mix[, .(Truth = Raji, type = "Raji")]),
                    cbind(sample_ests_ts_info[, .(Mix, Estimate = `THP-1`)],  Mix[, .(Truth = `THP-1`, type = "THP-1")]))
r_plot_dt_ts <- with(plot_dt_ts, cor(x = Truth, y = Estimate))
r_plot_dt_ts <- round(r_plot_dt_ts, 2)
cols <- RColorBrewer::brewer.pal(length(unique(plot_dt_ts$type)), "Dark2")
names(cols) <- unique(plot_dt_ts$type)
ggplot(data = plot_dt_ts, aes(x = Truth, y = Estimate, color = type, shape = type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x = 0.05, y = .4, label = paste("r", "=", r_plot_dt_ts)) + 
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 8),
        legend.position = c(0.8,0.5),
        legend.title = element_blank())
ggsave(file.path(".", "figs", "ts_abbas.png"), dpi = "retina")

# Spiked data -------------------------------------------------------------

mixes_abbas_spike <- fread(file.path(".", "cell_data", "abbas_et_al_data", "FinalMixtureMatrix.Spike.Abbas.coloncancer.subset.txt"))
mixes_abbas_hct <- fread(file.path(".", "cell_data", "abbas_et_al_data", "FinalMixtureMatrix.Abbas.HCT116.Res30.subset.txt"))

all <- seq(1, 160, by = 5) + 1


mysubset <- lapply(all, function(x) {
  names(mixes_abbas_spike)[x: (x+4)]
})

keep_idx <- which(mixes_abbas$`!Sample_title` %in% mixes_abbas_spike$Probe)
cell_type_data_var <- cell_type_data_var[keep_idx, ]
cell_type_data_mat <- cell_type_data[keep_idx, c("Jurkat", "IM-9", "Raji", "THP-1"), with = F]
cell_type_data <- cell_type_data[keep_idx,]
# all_dt <- lapply(mysubset, function(cn){
#   mixes_abbas_spike[, c("Probe", cn), with = F]
# })
# 
# out <- lapply(all_dt, function(analysis_dt){
#   my_sim_p_abbas2(p = 0, use.sva = T, iter = 1, tol = 1e-6, eval = T, use.sum = NULL, scaling = T, prop_data = cell_type_data_var, 
#                    dt = analysis_dt[, -1], lm2_dt = cell_type_data_mat, all_dt = analysis_dt)
# })
# dim(mixes_abbas_spike)


analysis_dt <- mixes_abbas_spike#all_dt[[1]]

out <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 1e3, tol = 1e-10, eval = F, use.sum = NULL, scaling = T, 
                prop_data = cell_type_data_var, all_dt = analysis_dt, raw_dt = cell_type_data)
t_info <- sapply(strsplit(colnames(mixes_abbas_spike)[-1], "\\("), "[", 2)
t_info <- sapply(strsplit(t_info, "\\)"), "[", 1)
s_info <-  sapply(strsplit(t_info, ","), "[", 1)
t_info <-  sapply(strsplit(t_info, ","), "[", 2)
t_info <- trimws(t_info);s_info <- trimws(s_info)

out[, `:=`(spike_info = sapply(strsplit(s_info, " "), "[", 1),
           tumor_info = sapply(strsplit(t_info, " "), "[", 1))]

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
  m <- paste0(x[1], x[2])
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


# Mixed/ Unknown content -------------------------------------------------------------------

all_objs <- ls()
all_objs <- all_objs[!all_objs %in% c("mixes_abbas", "cell_type_data_var", "cell_type_data",
                                      "eval_mix", "my_sim_p_abbas3")]
rm(list = all_objs)
mixes_abbas_hct <- fread(file.path(".", "cell_data", "abbas_et_al_data", "FinalMixtureMatrix.Abbas.HCT116.Res30.subset.txt"))
ciber_sig <- fread(file.path(".", "cell_data", "abbas_et_al_data", "GSE11103_matrix_classes.GSE11103_matrix_pure.bm.K999.0.txt"))
ciber_sig[, V6:=NULL]
names(ciber_sig) <- names(cell_type_data)
all <- seq(1, 64, by = 4) + 1


keep_idx <- which(mixes_abbas$`!Sample_title` %in% mixes_abbas_hct$Probeset)
cell_type_data_var <- cell_type_data_var[keep_idx, ]
cell_type_data_var <- cell_type_data_var[match(ciber_sig$`!Sample_title`,`!Sample_title`),]
cell_type_data_mat <- cell_type_data[keep_idx, c("Jurkat", "IM-9", "Raji", "THP-1"), with = F]

cell_type_data <- cell_type_data[keep_idx,]
cell_type_data <- cell_type_data[match(ciber_sig$`!Sample_title`,`!Sample_title`),]
# all_dt <- lapply(mysubset, function(cn){
#   mixes_abbas_spike[, c("Probe", cn), with = F]
# })
# 
# out <- lapply(all_dt, function(analysis_dt){
#   my_sim_p_abbas2(p = 0, use.sva = T, iter = 1, tol = 1e-6, eval = T, use.sum = NULL, scaling = T, prop_data = cell_type_data_var, 
#                    dt = analysis_dt[, -1], lm2_dt = cell_type_data_mat, all_dt = analysis_dt)
# })
# dim(mixes_abbas_spike)


analysis_dt <- mixes_abbas_hct[match(ciber_sig$`!Sample_title`,`Probeset`),]

# Without SVA
{
  # c_v <- rep(c(rep(0, 4), rep(1, 12)), 4)
  out <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 1e4, tol = 1e-10, eval = F, use.sum = NULL, scaling = T,
                         prop_data = cell_type_data_var, all_dt = analysis_dt, raw_dt = cell_type_data)
  # out <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 1e4, tol = 1e-10, eval = F, use.sum = NULL, scaling = T, 
  #                        prop_data = cell_type_data_var, all_dt = analysis_dt, raw_dt = ciber_sig, method = "TMM")
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
      mr <- mean(dt/t)
      out <- c("pearson" = r, "rmse" = d, "mean_ratio" = mr) 
    }else{
      out <- c("pearson" = NA, "rmse" = NA, "mean_ratio" = NA) 
    }
    
    return(out)
  })
  cor_out <- t(cor_out)
  out <- cbind(out, cor_out)
  fwrite(out, file.path(".", "sim_results", "mixes_out_no_sva.txt"))
}
outa <- out
# With SVA
{
  out <- my_sim_p_abbas3(p = 0, use.sva = T, iter = 1e4, tol = 1e-10, eval = F, use.sum = NULL, scaling = T,
                         prop_data = cell_type_data_var, all_dt = analysis_dt, raw_dt = cell_type_data)
  # out <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 1e4, tol = 1e-10, eval = F, use.sum = NULL, scaling = T, 
  #                        prop_data = cell_type_data_var, all_dt = analysis_dt, raw_dt = ciber_sig, method = "TMM")
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
      mr <- mean(dt/t)
      out <- c("pearson" = r, "rmse" = d, "mean_ratio" = mr) 
    }else{
      out <- c("pearson" = NA, "rmse" = NA, "mean_ratio" = NA) 
    }
    
    return(out)
  })
  cor_out <- t(cor_out)
  out <- cbind(out, cor_out)
  fwrite(out, file.path(".", "sim_results", "mixes_out.txt"))
}
cbind(out$mean_ratio, outa$mean_ratio)


# LM
{
  
  out <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 0, tol = 1e-10, eval = F, use.sum = NULL, scaling = T,
                         prop_data = cell_type_data_var, all_dt = analysis_dt, raw_dt = cell_type_data)
  # out <- my_sim_p_abbas3(p = 0, use.sva = F, iter = 1e4, tol = 1e-10, eval = F, use.sum = NULL, scaling = T, 
  #                        prop_data = cell_type_data_var, all_dt = analysis_dt, raw_dt = ciber_sig, method = "TMM")
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
      mr <- mean(dt/t)
      out <- c("pearson" = r, "rmse" = d, "mean_ratio" = mr) 
    }else{
      out <- c("pearson" = NA, "rmse" = NA, "mean_ratio" = NA) 
    }
    
    return(out)
  })
  cor_out <- t(cor_out)
  out <- cbind(out, cor_out)
  fwrite(out, file.path(".", "sim_results", "mixes_out_lm.txt"))
}