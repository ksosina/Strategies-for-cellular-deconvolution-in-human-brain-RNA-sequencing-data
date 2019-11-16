est_pi_c <- function(my_data_all = rna_gene_x_cells, yes_plot = F, scaling = F, cell_type_info = cell_type_data, ...){
  
  
  # Generate and use training data
  data_split <- gene_datasets(...)
  my_data_train <- my_data_all[, c(names(my_data)[1], data_split$train), with = F]
  my_vals <- est_ref_params(my_data_train, cell_type_info)$out
  
  # Load libs
  packs <- c("data.table", "dplyr", "mixtools", "SummarizedExperiment")
  libs_loaded <- sapply(packs, require, character.only = T)
  
  # Prep validation data
  my_val_data <- my_data_all[, c(names(my_data)[1], data_split$valid), with = F]
  my_val_data <- my_val_data[genes %in% my_vals$genes]
  y <- apply(my_val_data[,-1], 1, mean)
  my_val_data <- apply(my_val_data[,-1], 2, mean)
  
  
  # EM approach
  means <- apply(my_vals[,-1], 2, mean)
  sds <- apply(my_vals[,-1], 2, sd)
  wait1 <- normalmixEM(my_val_data, lambda = .2, mu = means, sigma =  sds, arbvar = T, maxrestarts = 1e3, maxit = 5e3)
  
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
    n.max <- which(names(s) == names(s.max))
    r <- 2:n
    r <- r[!r %in% n.max]
    
    # print(c(n.max,r))
    
    # Obtain all columns except the max
    mysubset <- lapply(names(my_vals)[r], function(x) {c(x, names(my_vals)[n.max])})
    
    
    x <- my_vals[, lapply(mysubset, function(x) get(x[1])/sd(get(x[1])) - get(x[2])/sd(get(x[2]))  )]
    # x <- scale(x, center = F)
    
    # y_star <- (y - unlist(my_vals[,n, with = F])/sd(unlist(my_vals[,n, with = F])) )
    y_star <- (y - unlist(my_vals[,n.max, with = F])/sd(unlist(my_vals[,n.max, with = F])) )
    
  }else{
    
    x <- my_vals[, lapply(mysubset, function(x) get(x[1]) - get(x[2])  )]
    y_star <- (y - unlist(my_vals[,n, with = F]))
    
    
  }
  
  x <- data.matrix(cbind(1,x))
  
  # print(lm(y_star ~ -1 + x[,-1]))
  
  
  # print(apply(x, 2, sd)[-1])
  
  
  # same syntax here
  # x <- sapply(2:(n-1), function(i){
  #   {unlist(my_vals[,i]) - unlist(my_vals[,n])}
  # })
  # x <- cbind(1, t(x))
  
  
  
  
  beta <- solve(crossprod(x)) %*% crossprod(x, y_star)
  beta <- c(beta)
  
  beta <- c(beta[-1], 1 - sum(beta[-1]))
  names(beta) <- names(my_vals)[2:(n)]
  if (scaling)
    names(beta) <- names(my_vals)[c(r, n.max)]
  
  
  # if(scaling){
  #   beta <- c(beta, 1 - sum(beta))
  #   names(beta) <- names(my_vals)[2:(n)]
  # }else{
  #   beta <- c(beta[-1], 1 - sum(beta[-1]))
  #   names(beta) <- names(my_vals)[2:(n)]
  # }
  
  
  
  # Truth 
  all_pi_c <- est_ref_params(my_data_all, cell_type_info)$pi_c
  
  
  re <- list(EM = wait1, LM = beta, truth = unique(all_pi_c))
  return(re)
}
