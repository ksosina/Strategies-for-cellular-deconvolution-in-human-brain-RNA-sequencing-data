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
                          r2 = unlist(cell_type_data_test$R2m))]


with(cell_type_data_var, plot(prop, (f_test), main = "Relationship between p values from DE and p", ylab = "F-test unadjst p vals", xlab = "p"))
with(cell_type_data_var, plot(prop, (r2), main = "Relationship between r2 and p", ylab = "r2", xlab = "p"))
with(cell_type_data_var, plot(r2, (f_test), main = "Relationship between p values from DE and p", ylab = "F-test unadjst p vals", xlab = "p"))

cell_type_data_var[, prop := r2]
cell_type_data_var[prop < .2,]

# Analyses ----------------------------------------------------------------


# Prep lamda
# # use log
# cell_type_data_lambda <-  cbind(cell_type_data[, 1],
#                                 1/apply(cell_type_data[, -1], 1, mean))

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




# estimate props
sample_ests_LM <- apply(analysis_dt, 2, function(y){
  
  # LM
  my_params <- cell_type_data_mat
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
  
  return(beta)
  
})
sample_ests_LM <- t(sample_ests_LM)

# estimate props
sample_ests_LM1 <- apply(mixes_abbas[,-1], 2, function(y){
  
  # LM
  my_params <- cell_type_data_mu[, -1]
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
  beta[beta < 0] <- 0
  beta[beta >1] <- 1
  
  return(beta)
  
})
sample_ests_LM1 <- t(sample_ests_LM1)

# subest to informative genes
p <- .95
cutoff <- quantile(cell_type_data_var$prop, probs = c(1-p,p))
info_idx <- which(cell_type_data_var$prop >= cutoff[2] & cell_type_data_var$`E(Var(Y|X))` > 0)


# Using lambda 
sample_ests_LM_info <- apply(analysis_dt, 2, function(y){
  
  # LM
  y <- y[info_idx]
  my_params <- cell_type_data_mat[info_idx, ]
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
  
  ft <- lm(y_star ~ -1 + x, weights = cell_type_data_var$prop[info_idx])
  # beta <- solve(crossprod(x)) %*% crossprod(x, y_star)
  beta <- coef(ft)
  beta <- c(beta)
  # if(scaling){
  #   beta <- c(beta, 1 - sum(beta))
  #   names(beta) <- names(my_vals)[2:(n)]
  # }
  beta <- c(beta[-1], 1 - sum(beta[-1]))
  names(beta) <- names(my_params)[1:(n)]
  beta[beta<0] <- 0
  beta[beta>1] <- 1
  
  return(beta)
  
})
sample_ests_LM_info <- t(sample_ests_LM_info)
sample_ests_LM_info <- data.table("ID" = sapply(strsplit(row.names(sample_ests_LM_info), " "), "[", 1),
                                  "Mix" = sapply(strsplit(row.names(sample_ests_LM_info), " "), "[", 2),
                                  sample_ests_LM_info)
sample_ests_LM_info

# Using normal ls 
my_params <- cell_type_data_mu[info_idx, -1]
n <- ncol(my_params)

# Obtain all columns except the last
mysubset <- lapply(names(my_params)[1:(n-1)], function(x) {c(x, names(my_params)[n])})
x <- my_params[, lapply(mysubset, function(x) get(x[1])-get(x[2]))]
x <- data.matrix(cbind(1,x))

sample_ests_LM1_info_a <- lm(data.matrix(mixes_abbas[info_idx,-1]) ~ 0 + x, offset = unlist(my_params[,n, with = F]), weights = cell_type_data_var$prop[info_idx])
sample_ests_LM1_info_a <- t(sample_ests_LM1_info_a$coefficients)
colnames(sample_ests_LM1_info_a) <- names(my_params)[1:(n)]
sample_ests_LM1_info_a <- apply(sample_ests_LM1_info_a, 1, function(beta){
  beta <- c(beta[-1], 1 - sum(beta[-1]))
  names(beta) <- names(my_params)[1:(n)]
  beta[beta < 0] <- 0
  beta[beta >1] <- 1
  beta <- beta/ sum( beta )
  return(beta)
})
sample_ests_LM1_info_a <- t(sample_ests_LM1_info_a)

sample_ests_LM1_info_b <- lm(data.matrix(mixes_abbas[info_idx,-1]) ~ 0 + x, offset = unlist(my_params[,n, with = F]), weights = (cell_type_data_var$w_2[info_idx]))
sample_ests_LM1_info_b <- t(sample_ests_LM1_info_b$coefficients)
colnames(sample_ests_LM1_info_b) <- names(my_params)[1:(n)]
sample_ests_LM1_info_b <- apply(sample_ests_LM1_info_b, 1, function(beta){
  beta <- c(beta[-1], 1 - sum(beta[-1]))
  names(beta) <- names(my_params)[1:(n)]
  beta[beta < 0] <- 0
  beta[beta >1] <- 1
  beta <- beta/ sum( beta )
  return(beta)
})
sample_ests_LM1_info_b <- t(sample_ests_LM1_info_b)


sample_ests_LM1_info_c <- lm(data.matrix(mixes_abbas[,-1]) ~ 0 + x, offset = unlist(my_params[,n, with = F]))
sample_ests_LM1_info_c <- t(sample_ests_LM1_info_c$coefficients)
colnames(sample_ests_LM1_info_c) <- names(my_params)[1:(n)]
sample_ests_LM1_info_c <- apply(sample_ests_LM1_info_c, 1, function(beta){
  beta <- c(beta[-1], 1 - sum(beta[-1]))
  names(beta) <- names(my_params)[1:(n)]
  beta[beta < 0] <- 0
  beta[beta >1] <- 1
  beta <- beta/ sum( beta )
  return(beta)
})
sample_ests_LM1_info_c <- t(sample_ests_LM1_info_c)


sample_ests_LM1_info <- apply(mixes_abbas[,-1], 2, function(y){
  
  # LM
  y <- y[info_idx]
  my_params <- cell_type_data_mu[info_idx, -1]
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
  beta[beta < 0] <- 0
  beta[beta >1] <- 1
  
  return(beta)
  
})
sample_ests_LM1_info <- t(sample_ests_LM1_info)
sample_ests_LM1_info <- data.table("ID" = sapply(strsplit(row.names(sample_ests_LM1_info), " "), "[", 1),
                                  "Mix" = sapply(strsplit(row.names(sample_ests_LM1_info), " "), "[", 2),
                                  sample_ests_LM1_info)
sample_ests_LM1_info

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

eval_mix(sample_ests_LM_info)
eval_mix(sample_ests_LM1_info_a)
eval_mix(sample_ests_LM1_info_b)
eval_mix(sample_ests_LM1_info_c)

pear1 <- sapply(1:nrow(sample_ests_LM1_info), function(i){
  dt <- sample_ests_LM1_info[i,]
  m <- dt$Mix
  dt <- dt[, -c(1:2), with = F]
  dt <- unlist(dt)
  t <- Mix[Mix == m, -1, with = F]
  t <- unlist(t)
  r <- cor(dt, t)
  return(r)
})

pear1 <- tapply(pear1, rep(paste0("Mix", LETTERS[1:4]), each = 3), mean)
pear1 <- data.table(type = names(pear1), pearson = (pear1))

#sim across all values of p
p_sims <- lapply(seq(0,.999, by = 1e-3), function(p){
  # print(p)
  # Estimate
  cutoff <- quantile(cell_type_data_var$prop, probs = c(1-p,p))
  cutoff[2] <- ifelse(p == 0, 0, cutoff[2])
  info_idx <- which(cell_type_data_var$prop >= cutoff[2] )
  
  sample_ests_LM_info <- apply(analysis_dt, 2, function(y){
    
    # LM
    y <- y[info_idx]
    my_params <- cell_type_data_mat[info_idx, ]
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
    beta[beta < 0] <- 0
    beta[beta > 1] <- 1
    
    return(beta)
    
  })
  sample_ests_LM_info <- t(sample_ests_LM_info)
  sample_ests_LM_info <- data.table("ID" = sapply(strsplit(row.names(sample_ests_LM_info), " "), "[", 1),
                                    "Mix" = sapply(strsplit(row.names(sample_ests_LM_info), " "), "[", 2),
                                    sample_ests_LM_info)
  

  # Evaluate
  pear <- eval_mix(sample_ests_LM_info)
  pear <- data.table(pear, "p" = p, ngenes = length(info_idx))
  return(pear)
})
p_sims <- do.call(rbind, p_sims)


p_sims_LM1 <- lapply(seq(0,.999, by = 1e-3), function(p){
  # print(p)
  # Estimate
  cutoff <- quantile(cell_type_data_var$prop, probs = c(1-p,p))
  cutoff[2] <- ifelse(p == 0, 0, cutoff[2])
  info_idx <- which(cell_type_data_var$prop >= cutoff[2] )
  
  sample_ests_LM1_info <- apply(mixes_abbas[,-1], 2, function(y){
    
    # LM
    y <- y[info_idx]
    my_params <- cell_type_data_mu[info_idx, -1]
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
    beta[beta < 0] <- 0
    beta[beta >1] <- 1
    
    return(beta)
    
  })
  sample_ests_LM1_info <- t(sample_ests_LM1_info)
  sample_ests_LM1_info <- data.table("ID" = sapply(strsplit(row.names(sample_ests_LM1_info), " "), "[", 1),
                                    "Mix" = sapply(strsplit(row.names(sample_ests_LM1_info), " "), "[", 2),
                                    sample_ests_LM1_info)
  
  
  # Evaluate
  pear <- eval_mix(sample_ests_LM1_info)
  pear <- data.table(pear, "p" = p, ngenes = length(info_idx))
  return(pear)
})
p_sims_LM1 <- do.call(rbind, p_sims_LM1)

# two step
my_sim_p_abbas <- function(p, use.sva = T, iter = 1, tol = 1e-6, eval = T){
  # Subset to informative genes
  # Estimate
  cutoff <- quantile(cell_type_data_var$prop, probs = c(1-p,p))
  cutoff[2] <- ifelse(p == 0, 0, cutoff[2])
  info_idx <- which(cell_type_data_var$prop >= cutoff[2] )
  
  # Get initial estimates of pi
  sample_ests_LM2 <- apply(analysis_dt, 2, function(y){
    
    # LM
    y <- y[info_idx]
    my_params <- cell_type_data_mat[info_idx, ]
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
      beta[beta < 0] <- 0
      beta[beta > 0] <- beta[beta >0] /sum( beta[beta >0] )
    }
    
    
    return(beta)
    
  })
  sample_ests_LM2 <- t(sample_ests_LM2)
  sample_ests_LM2_old <- sample_ests_LM2
  
  
  # SSVA to eliminate noise
  x <- sample_ests_LM2[, -1]
  x <- cbind(1, x)
  if(use.sva){
    mod1 <-  x
    mod0 <- cbind(mod1[,1])
    svseq <- svaseq(data.matrix(mixes_abbas[,-1]),mod1,mod0, constant = 1)$sv
  }
  
  gene_ests_LM2 <- apply(mixes_abbas[,-1], 1, function(y){
    
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
  
  # Get new estimates of pi
  sample_ests_LM2 <- apply(mixes_abbas[,-1], 2, function(y){
    
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
      beta[beta < 0] <- 0
      beta[beta >0] <- beta[beta >0] /sum( beta[beta >0] )
    }
    
    
    return(beta)
    
  })
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
      gene_ests_LM2 <- apply(mixes_abbas[,-1], 1, function(y){
        
        # LM
        n <- ncol(sample_ests_LM2)
        # if(use.sva){
        #   ft <- lm(y ~ -1 + x + svseq)
        # }else{
        #   ft <- lm(y ~ -1 + x )
        # }
        ft <- lm(y ~ -1 + x )
        beta <- coef(ft)
        beta <- beta[1:n]
        beta_1 <- sapply(beta[-1], function(b) {b <- b + beta[1]; return(b)} )
        b <- c(beta[1], beta_1)
        names(b) <- colnames(sample_ests_LM2)
        return(b)
        
      })
      gene_ests_LM2 <- t(gene_ests_LM2)
      sample_ests_LM2 <- apply(mixes_abbas[,-1], 2, function(y){
        
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
          beta[beta < 0] <- 0
          beta[beta >0] <- beta[beta >0] /sum( beta[beta >0] )
        }
        
        return(beta)
        
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
  
  
  
  # Evaluate
  sample_ests_LM2 <- data.table("ID" = sapply(strsplit(row.names(sample_ests_LM2), " "), "[", 1),
                                    "Mix" = sapply(strsplit(row.names(sample_ests_LM2), " "), "[", 2),
                                sample_ests_LM2)
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
  
  pear <- tapply(pear, rep(paste0("Mix", LETTERS[1:4]), each = 3), mean)
  pear <- data.table(type = names(pear), pearson = (pear), "p" = p, ngenes = length(info_idx))
  
  if(eval){
    return(pear)
  }else{
    return(sample_ests_LM2)
  }
  
  
}
my_sim_p_abbas2 <- function(p, use.sva = T, iter = 1, tol = 1e-6, eval = T){
  # Subset to informative genes from reference
  # Estimate
  cutoff <- quantile(cell_type_data_var$prop, probs = c(1-p,p))
  cutoff[2] <- ifelse(p == 0, 0, cutoff[2])
  info_idx <- which(cell_type_data_var$prop >= cutoff[2] )
  
  # Get initial estimates of pi
  sample_ests_LM2 <- apply(analysis_dt, 2, function(y){
    
    # LM
    y <- y[info_idx]
    my_params <- cell_type_data_mat[info_idx, ]
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
      beta[beta < 0] <- 0
      beta[beta > 0] <- beta[beta >0] /sum( beta[beta >0] )
    }
    
    
    return(beta)
    
  })
  sample_ests_LM2 <- t(sample_ests_LM2)
  sample_ests_LM2_old <- sample_ests_LM2
  
  
  # SSVA to eliminate noise
  x <- sample_ests_LM2[, -1]
  x <- cbind(1, x)
  if(use.sva){
    mod1 <-  x
    mod0 <- cbind(mod1[,1])
    svseq <- svaseq(data.matrix(mixes_abbas[,-1]),mod1,mod0, constant = 1)$sv
  }
  
  gene_ests_LM2 <- apply(mixes_abbas[,-1], 1, function(y){
    
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
    
    # test
    
    if(use.sva){
      ft.anova <-  anova(ft)
      ft.anova <- ft.anova$`Pr(>F)`[1]
    }else{
      ft.anova <-  anova(ft)
      ft.anova <- ft.anova$`Pr(>F)`[1]
    }
    b <- c(b, "pval" = ft.anova)
    return(b)
    
  })
  gene_ests_LM2 <- t(gene_ests_LM2)
  ncol_gene <- ncol(gene_ests_LM2)
  pval_gene <- gene_ests_LM2[, ncol_gene]
  qval_gene <- p.adjust(pval_gene, method = "bon")
  qval_gene <- which(qval_gene <= 5e-2)
  info_idx <- intersect(info_idx, qval_gene)
  if(length(info_idx) == 0){
    stop("Not enough genes")
  }
  gene_ests_LM2 <- gene_ests_LM2[, -ncol_gene]
  
  # Get new estimates of pi
  sample_ests_LM2 <- apply(mixes_abbas[,-1], 2, function(y){
    
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
      beta[beta < 0] <- 0
      beta[beta >0] <- beta[beta >0] /sum( beta[beta >0] )
    }
    
    
    return(beta)
    
  })
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
      gene_ests_LM2 <- apply(mixes_abbas[,-1], 1, function(y){
        
        # LM
        n <- ncol(sample_ests_LM2)
        # if(use.sva){
        #   ft <- lm(y ~ -1 + x + svseq)
        # }else{
        #   ft <- lm(y ~ -1 + x )
        # }
        ft <- lm(y ~ -1 + x )
        beta <- coef(ft)
        beta <- beta[1:n]
        beta_1 <- sapply(beta[-1], function(b) {b <- b + beta[1]; return(b)} )
        b <- c(beta[1], beta_1)
        names(b) <- colnames(sample_ests_LM2)
        return(b)
        
      })
      gene_ests_LM2 <- t(gene_ests_LM2)
      sample_ests_LM2 <- apply(mixes_abbas[,-1], 2, function(y){
        
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
          beta[beta < 0] <- 0
          beta[beta >0] <- beta[beta >0] /sum( beta[beta >0] )
        }
        
        return(beta)
        
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
  
  
  
  # Evaluate
  sample_ests_LM2 <- data.table("ID" = sapply(strsplit(row.names(sample_ests_LM2), " "), "[", 1),
                                "Mix" = sapply(strsplit(row.names(sample_ests_LM2), " "), "[", 2),
                                sample_ests_LM2)
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
  
  pear <- tapply(pear, rep(paste0("Mix", LETTERS[1:4]), each = 3), mean)
  pear <- data.table(type = names(pear), pearson = (pear), "p" = p, ngenes = length(info_idx))
  
  if(eval){
    return(pear)
  }else{
    return(sample_ests_LM2)
  }
  
  
}
p_sims_ts <- lapply(seq(0,.99, by = 1e-2), my_sim_p_abbas, use.sva = F, iter = 2)
p_sims_ts <- do.call(rbind, p_sims_ts)

sample_ests_ts_info <- my_sim_p_abbas(p, iter = 2, eval = F, use.sva = F)

Mix <- inner_join(sample_ests_LM_info[, .(Mix)],Mix) %>% data.table

# Post normalization
ref_mean <- apply(cell_type_data[, -c(1,6), with = F], 1, mean)
yExprs_normfactors <- edgeR::calcNormFactors(cbind(data.matrix(mixes_abbas[, -1]), ref_mean), method = "TMM")
n_y <- length(yExprs_normfactors)
cell_type_data_lambda <-  cbind(cell_type_data[, 1],
                                1/apply(cell_type_data[, -c(1,6), with = F], 1, function(x) {
                                  x <- x * yExprs_normfactors[n_y]; return(mean(x))}))

yExprs_normfactors <- yExprs_normfactors[-n_y]
yExprs_normfactors <- replicate(nrow(mixes_abbas)  , yExprs_normfactors)
yExprs_normfactors <- t(yExprs_normfactors)
yExprs_2 <- mixes_abbas[,-1]*yExprs_normfactors

# prep mixes data
analysis_dt_tmm <- cbind(yExprs_2, cell_type_data_lambda[,-1]) %>% data.table
n <- ncol(analysis_dt_tmm)
mysubset <- lapply(names(analysis_dt_tmm)[1:(n-1)], function(x) {c(x, names(analysis_dt_tmm)[n])})
analysis_dt_tmm <- analysis_dt_tmm[, lapply(mysubset, function(x) (get(x[1])) * (get(x[2])) )]
names(analysis_dt_tmm) <- names(yExprs_2)

# estimate props
tmm_LM <- apply(analysis_dt_tmm, 2, function(y){
  
  # LM
  my_params <- cell_type_data_mat
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
  
  return(beta)
  
})
tmm_LM <- t(tmm_LM)
tmm_LM

# subest to informative genes
cutoff <- quantile(cell_type_data_var$prop, probs = c(1-p,p))
info_idx <- which(cell_type_data_var$prop >= cutoff[2] & cell_type_data_var$`E(Var(Y|X))` > 0)

tmm_LM_info <- apply(analysis_dt_tmm, 2, function(y){
  
  # LM
  y <- y[info_idx]
  my_params <- cell_type_data_mat[info_idx, ]
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
  beta[beta < 0] <- 0
  names(beta) <- names(my_params)[1:(n)]
  
  return(beta)
  
})
tmm_LM_info <- t(tmm_LM_info)
tmm_LM_info <- data.table("ID" = sapply(strsplit(row.names(tmm_LM_info), " "), "[", 1),
                                  "Mix" = sapply(strsplit(row.names(tmm_LM_info), " "), "[", 2),
                                  tmm_LM_info)
tmm_LM_info

Mix_tmm <- rbind(MixA, MixB, MixC, MixD)
Mix_tmm <- data.table("Mix" = paste0("Mix", LETTERS[1:4]), Mix_tmm)

pear_tmm <- sapply(1:nrow(tmm_LM_info), function(i){
  dt <- tmm_LM_info[i,]
  m <- dt$Mix
  dt <- dt[, -c(1:2), with = F]
  dt <- unlist(dt)
  t <- Mix_tmm[Mix == m, -1, with = F]
  t <- unlist(t)
  r <- cor(dt, t)
  return(r)
})

pear_tmm <- tapply(pear_tmm, rep(paste0("Mix", LETTERS[1:4]), each = 3), mean)
pear_tmm <- data.table(type = names(pear_tmm), pearson = (pear_tmm))
pear_tmm

#sim across all values of p
p_sims_tmm <- lapply(seq(0,.99, by = 1e-2), function(p){
  # print(p)
  # Estimate
  cutoff <- quantile(cell_type_data_var$prop, probs = c(1-p,p))
  cutoff[2] <- ifelse(p == 0, 0, cutoff[2])
  info_idx <- which(cell_type_data_var$prop >= cutoff[2] )
  
  tmm_LM_info <- apply(analysis_dt_tmm, 2, function(y){
    
    # LM
    y <- y[info_idx]
    my_params <- cell_type_data_mat[info_idx, ]
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
    
    return(beta)
    
  })
  tmm_LM_info <- t(tmm_LM_info)
  tmm_LM_info <- data.table("ID" = sapply(strsplit(row.names(tmm_LM_info), " "), "[", 1),
                                    "Mix" = sapply(strsplit(row.names(tmm_LM_info), " "), "[", 2),
                                    tmm_LM_info)
  
  
  # Evaluate
  pear_tmm <- sapply(1:nrow(tmm_LM_info), function(i){
    dt <- tmm_LM_info[i,]
    m <- dt$Mix
    dt <- dt[, -c(1:2), with = F]
    dt <- unlist(dt)
    t <- Mix_tmm[Mix == m, -1, with = F]
    t <- unlist(t)
    r <- cor(dt, t)
    return(r)
  })
  
  pear_tmm <- tapply(pear_tmm, rep(paste0("Mix", LETTERS[1:4]), each = 3), mean)
  pear_tmm <- data.table(type = names(pear_tmm), pearson = (pear_tmm), "p" = p, ngenes = length(info_idx))
  return(pear_tmm)
})
p_sims_tmm <- do.call(rbind, p_sims_tmm)

Mix_tmm <- inner_join(tmm_LM_info[, .(Mix)],Mix) %>% data.table

fwrite(p_sims, file.path(".", "sim_results", "p_sims_abbas.txt"))
fwrite(p_sims_LM1, file.path(".", "sim_results", "p_sims_LM1_abbas.txt"))
fwrite(p_sims_ts, file.path(".", "sim_results", "p_sims_TS_abbas.txt"))
# plot --------------------------------------------------------------------


transparent_legend =  theme(
  legend.background = element_rect(fill ="transparent"),
  legend.key = element_rect(fill = "transparent",
                            color = "transparent")
)

remove_grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))

no_x_axis_label <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())


plot_dt <- rbind(inner_join(sample_ests_LM_info[, .(Mix, Estimate = Jurkat)],  Mix[, .(Mix, Truth = Jurkat, type = "Jurkat")]),
                               inner_join(sample_ests_LM_info[, .(Mix, Estimate = `IM-9`)],  Mix[, .(Mix, Truth = `IM-9`, type = "IM-9")]),
                               inner_join(sample_ests_LM_info[, .(Mix, Estimate = Raji)],  Mix[, .(Mix, Truth = Raji, type = "Raji")]),
                               inner_join(sample_ests_LM_info[, .(Mix, Estimate = `THP-1`)],  Mix[, .(Mix, Truth = `THP-1`, type = "THP-1")])) %>% data.table
plot_dt_LM1 <- rbind(inner_join(sample_ests_LM1_info[, .(Mix, Estimate = Jurkat)],  Mix[, .(Mix, Truth = Jurkat, type = "Jurkat")]),
                     inner_join(sample_ests_LM1_info[, .(Mix, Estimate = `IM-9`)],  Mix[, .(Mix, Truth = `IM-9`, type = "IM-9")]),
                     inner_join(sample_ests_LM1_info[, .(Mix, Estimate = Raji)],  Mix[, .(Mix, Truth = Raji, type = "Raji")]),
                     inner_join(sample_ests_LM1_info[, .(Mix, Estimate = `THP-1`)],  Mix[, .(Mix, Truth = `THP-1`, type = "THP-1")])) %>% data.table
plot_dt_ts <- rbind(inner_join(sample_ests_ts_info[, .(Mix, Estimate = Jurkat)],  Mix[, .(Mix, Truth = Jurkat, type = "Jurkat")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = `IM-9`)],  Mix[, .(Mix, Truth = `IM-9`, type = "IM-9")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = Raji)],  Mix[, .(Mix, Truth = Raji, type = "Raji")]),
                    inner_join(sample_ests_ts_info[, .(Mix, Estimate = `THP-1`)],  Mix[, .(Mix, Truth = `THP-1`, type = "THP-1")])) %>% data.table
fwrite(plot_dt, file.path(".", "sim_results", "plot_dt_abbas.txt"))
fwrite(plot_dt_LM1, file.path(".", "sim_results", "plot_dt_LM1_abbas.txt"))
fwrite(plot_dt_ts, file.path(".", "sim_results", "plot_dt_TS_abbas.txt"))


cibersort <- fread(file.path(".", "cell_data", "cibersort_data","CIBERSORT.Output_Job2.txt"))

cibersort <- data.table("ID" = sapply(strsplit(cibersort$`Input Sample`, " "), "[", 1),
                          "Mix" = sapply(strsplit(cibersort$`Input Sample`, " "), "[", 2),
                          cibersort[,-1])

plot_dt_cibersort <- rbind(inner_join(cibersort[, .(Mix, Estimate = Jurkat)],  Mix[, .(Mix, Truth = Jurkat, type = "Jurkat")]),
                           inner_join(cibersort[, .(Mix, Estimate = `IM-9`)],  Mix[, .(Mix, Truth = `IM-9`, type = "IM-9")]),
                           inner_join(cibersort[, .(Mix, Estimate = Raji)],  Mix[, .(Mix, Truth = Raji, type = "Raji")]),
                           inner_join(cibersort[, .(Mix, Estimate = `THP-1`)],  Mix[, .(Mix, Truth = `THP-1`, type = "THP-1")]))
fwrite(plot_dt_cibersort, file.path(".", "sim_results", "cibersort_abbas.txt"))

plot_dt_tmm <- rbind(inner_join(tmm_LM_info[, .(Mix, Estimate = Jurkat)],  Mix_tmm[, .(Mix, Truth = Jurkat, type = "Jurkat")]),
                     inner_join(tmm_LM_info[, .(Mix, Estimate = `IM-9`)],  Mix_tmm[, .(Mix, Truth = `IM-9`, type = "IM-9")]),
                     inner_join(tmm_LM_info[, .(Mix, Estimate = Raji)],  Mix_tmm[, .(Mix, Truth = Raji, type = "Raji")]),
                     inner_join(tmm_LM_info[, .(Mix, Estimate = `THP-1`)],  Mix_tmm[, .(Mix, Truth = `THP-1`, type = "THP-1")]))

cols <- RColorBrewer::brewer.pal(length(unique(plot_dt$type)), "Dark2")
names(cols) <- unique(plot_dt$type)

r_plot_dt <- with(plot_dt, cor(x = Truth, y = Estimate))
r_plot_dt_tmm <- with(plot_dt_tmm, cor(x = Truth, y = Estimate))
r_plot_dt_ts <- with(plot_dt_ts, cor(x = Truth, y = Estimate))
r_plot_dt_cibersort <- with(plot_dt_cibersort, cor(x = Truth, y = Estimate))
r_plot_dt <- round(r_plot_dt, 2)
r_plot_dt_tmm <- round(r_plot_dt_tmm, 2)
r_plot_dt_ts <- round(r_plot_dt_ts, 2)
r_plot_dt_cibersort <- round(r_plot_dt_cibersort, 2)

ggplot(data = plot_dt, aes(x = Truth, y = Estimate, color = type, shape = type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x = 0.05, y = .4, label = paste("r", "=", r_plot_dt)) + 
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 8),
        legend.position = c(0.8,0.5),
        legend.title = element_blank())
ggsave(file.path(".", "figs", "abbas_95.png"), dpi = "retina")



ggplot(data = plot_dt_cibersort, aes(x = Truth, y = Estimate, color = type, shape = type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x = 0.05, y = .4, label = paste("r", "=", r_plot_dt_cibersort)) + 
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 8),
        legend.position = c(0.8,0.5),
        legend.title = element_blank())

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
ggsave(file.path(".", "figs", "ts_abbas_95.png"), dpi = "retina")

ggplot(data = plot_dt_tmm, aes(x = Truth, y = Estimate, color = type, shape = type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 8),
        legend.position = c(0.8,0.5),
        legend.title = element_blank())
ggsave(file.path(".", "figs","abbas_95_tmm.png"), dpi = "retina")


cols <- RColorBrewer::brewer.pal(length(unique(p_sims$type)), "Dark2")
names(cols) <- unique(p_sims$type)

ggplot(data = p_sims, aes(x = p, y = pearson, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  ylab("Pearson correlation coefficient") + xlab("CDF") +
  ggtitle("Relationship between cutoff and model fit") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.35,0.45),
        legend.title = element_blank())
ggsave(file.path(".", "figs","effect_of_cutoff.png"), dpi = "retina")


ggplot(data = p_sims_tmm, aes(x = p, y = pearson, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  ylab("Pearson correlation coefficient") + xlab("CDF") +
  ggtitle("Relationship between cutoff and model fit") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.35,0.45),
        legend.title = element_blank())
ggsave(file.path(".", "figs","effect_of_cutoff_tmm.png"), dpi = "retina")
