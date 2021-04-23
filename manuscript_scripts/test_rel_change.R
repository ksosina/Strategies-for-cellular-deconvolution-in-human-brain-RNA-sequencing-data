# Preamble ----------------------------------------------------------------

packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer", 
           "mixtools","matrixStats", "MuSiC", "convert", "xbioc", "ggplot2", "sva", "plotly",
           "doParallel", "parallel")
libs_loaded <- sapply(packs, library, character.only = T)

load(file.path(".", "all_results", "all_analysis.RData"))


# Functions ---------------------------------------------------------------

boot_rel_change <- function(dt1, dt2, B = 1e5, seed = NULL, show.progress = T){
  require(data.table)
  require(dplyr)
  # check same samples 
  s1 <- dt1$samples;s2 <- dt2$samples
  s1 <- sort(s1);s2 <- sort(s2)
  if(!identical(s1, s2)){
    stop("Different samples in each pop")
  }
  s <- s1
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  
  if(show.progress){
    pb <- txtProgressBar(min = 0, max = B , width = NA, style = 3)
    pb_i <- 0
  }
  
  bt_dt <- lapply(1:B, function(r){
    dt_samp <- data.table("samples" = sample(s, size = 223, replace = T))
    # dt_samp1 <- inner_join(dt_samp, dt1, by = "samples") %>% data.table
    # dt_samp2 <- inner_join(dt_samp, dt2, by = "samples") %>% data.table
    
    idx1 <- which(dt1$samples %in% dt_samp$samples)
    times1 <- tapply(dt_samp$samples, dt_samp$samples, length)
    times1 <- times1[match(names(times1), dt1$samples[idx1])]
    # identical(dt1$samples[idx1], names(times1))
    
    idx1 <- rep(idx1, times = times1) # account for sampling w replacement
    
    dt_samp1 <- dt1[idx1]
    
    idx2 <- which(dt2$samples %in% dt_samp$samples)
    # times2 <- tapply(dt2$samples, dt_samp$samples, length)
    times2 <- times1[match(names(times1), dt2$samples[idx2])]
    # identical(dt2$samples[idx2], names(times2))
    idx2 <- rep(idx2, times = times2) # account for sampling w replacement
    
    dt_samp2 <- dt2[idx2]
    
    mytable1 <- dt_samp1[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                          RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = type]
    
    
    mytable2 <- dt_samp2[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                          RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = type]
    
    if(show.progress){
      pb_i <<- pb_i+1
      setTxtProgressBar(pb, pb_i, title = paste(round(pb_i/B)*100,"% done\n"))
    }
    
   out <- data.table("rel_change_rsqrd" = (mytable2$R_sqrd - mytable1$R_sqrd)/mytable1$R_sqrd,
                     "rel_change_rmse" = (mytable2$RMSE - mytable1$RMSE)/mytable1$RMSE)
   
   
   
   
   
   return(out)
    
  
  })
   
  bt_dt <- do.call(rbind, bt_dt)
  
  m1 <- dt1[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                           RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = type]
  
  m2 <- dt2[, .(R_sqrd = cor(est, `Houseman DNAm-based`)^2,
                           RMSE = sqrt(mean((est-`Houseman DNAm-based`)^2))), by = type]
  
  obsrd_dt <- data.table("rel_change_rsqrd" = (m2$R_sqrd - m1$R_sqrd)/m1$R_sqrd,
                    "rel_change_rmse" = (m2$RMSE - m1$RMSE)/m1$RMSE)
  
  
  # Larger deviations of each bt sample from its mean (i.e obsrvd R sqrd or rmse)
  # See Two Guidelines for Bootstrap Hypothesis Testing (Hall 1991)
  p1 <- sum(abs(bt_dt$rel_change_rsqrd - obsrd_dt$rel_change_rsqrd)  >= abs(obsrd_dt$rel_change_rsqrd) ) + 1
  p1 <- p1/(B+1)
  p2 <- sum(abs(bt_dt$rel_change_rmse - obsrd_dt$rel_change_rmse) >= abs(obsrd_dt$rel_change_rmse) ) + 1
  p2 <- p2/(B+1)
  
  res <- list("R_sqrd" = c("observed relative_change" = abs(obsrd_dt$rel_change_rsqrd),
                           "pval" = p1),
              "RMSE" = c("observed relative_change" = abs(obsrd_dt$rel_change_rmse),
                         "pval" = p2))
  
  return(res)
}


# Testing relative change -------------------------------------------------

# tt <- inner_join(music_est_nnls, dar_nnls) %>% 
#   inner_join(data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)) %>% data.table
# 
# 
# cor(tt[, .(`TRUE`,neu_p, `Houseman DNAm-based`)])^2
# 
# t_res <- with(tt, t.test(x = neu_p, y = `Houseman DNAm-based`, paired = T))
# t_res$p.value


darmansis_default_test_dt <- inner_join(music_est[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default]")],
                                        data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)) %>%  data.table

dar_w_nac_cs <- inner_join(my_music_est_nw[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default:Nac all genes]")],
                           data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)) %>%  data.table

dar_w_nac_cs_50 <- inner_join(my_music_est_nw_50[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default:Nac top 50 genes]")],
                              data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)) %>%  data.table

dar_w_nac_cs_25 <- inner_join(my_music_est_nw_25[, .(samples, est = `TRUE`, type = "Ref:[Darmanis];Size:[Default:Nac top 25 genes]")],
                              data.table(samples = rse_gene$SampleID, "Houseman DNAm-based" = rse_gene$NeuN_pos_DNAm)) %>% data.table


# Darmanis default vs Darmanis ref + NAc all genes cell size
dVnc_all <- boot_rel_change(darmansis_default_test_dt, dar_w_nac_cs, B = 1e4, seed = 123)

# Darmanis default vs Darmanis ref + NAc all genes cell size
dVnc_50 <- boot_rel_change(darmansis_default_test_dt, dar_w_nac_cs_50, B = 1e4, seed = 123)

# Darmanis default vs Darmanis ref + NAc all genes cell size
dVnc_25 <- boot_rel_change(darmansis_default_test_dt, dar_w_nac_cs_25, B = 1e4, seed = 123)

dVnc_all;dVnc_50;dVnc_25

# Bland Altman plots ------------------------------------------------------

library(ggplot2)
plt_dt <- rbind(darmansis_default_test_dt, dar_w_nac_cs,
                dar_w_nac_cs_50, dar_w_nac_cs_25)

ggplot(data = plt_dt, aes(x = (est + `Houseman DNAm-based`)/2, y = est - `Houseman DNAm-based`,  color = type) ) +
  geom_point(size = 3)
