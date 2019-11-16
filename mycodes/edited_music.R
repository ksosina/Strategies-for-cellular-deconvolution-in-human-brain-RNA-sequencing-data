
my_music_basis <- function (x, non.zero = TRUE, markers = NULL, clusters, samples, 
          select.ct = NULL, ct.cov = FALSE, verbose = TRUE, my_abundance = NULL) 
{
  if (!is.null(select.ct)) {
    s.ct = sampleNames(x)[as.character(pVar(x, clusters)) %in% 
                            select.ct]
    x <- x[, s.ct, drop = FALSE]
  }
  if (non.zero) {
    nz.gene = rownames(x)[(rowSums(exprs(x)) != 0)]
    x <- x[nz.gene, , drop = FALSE]
  }
  clusters <- as.character(pVar(x, clusters))
  samples <- as.character(pVar(x, samples))
  M.theta <- sapply(unique(clusters), function(ct) {
    my.rowMeans(sapply(unique(samples), function(sid) {
      y = exprs(x)[, clusters %in% ct & samples %in% sid, 
                   drop = FALSE]
      rowSums(y)/sum(y)
    }), na.rm = TRUE)
  })
  if (verbose) {
    message("Creating Relative Abudance Matrix...")
  }
  if (ct.cov) {
    nGenes = nrow(x)
    n.ct = length(unique(clusters))
    nSubs = length(unique(samples))
    Theta <- sapply(unique(clusters), function(ct) {
      sapply(unique(samples), function(sid) {
        y = exprs(x)[, clusters %in% ct & samples %in% 
                       sid, drop = FALSE]
        return(rowSums(y)/sum(y))
      })
    })
    if (!is.null(select.ct)) {
      m.ct = match(select.ct, colnames(Theta))
      Theta = Theta[, m.ct]
    }
    Sigma.ct = sapply(1:nGenes, function(g) {
      sigma.temp = Theta[nGenes * (0:(nSubs - 1)) + g, 
                         ]
      Cov.temp = cov(sigma.temp)
      Cov.temp1 = cov(sigma.temp[rowSums(is.na(Theta[nGenes * 
                                                       (0:(nSubs - 1)) + 1, ])) == 0, ])
      Cov.temp[which(colSums(is.na(sigma.temp)) > 0), 
               ] = Cov.temp1[which(colSums(is.na(sigma.temp)) > 
                                     0), ]
      Cov.temp[, which(colSums(is.na(sigma.temp)) > 0)] = Cov.temp1[, 
                                                                    which(colSums(is.na(sigma.temp)) > 0)]
      return(Cov.temp)
    })
    colnames(Sigma.ct) = rownames(x)
    if (!is.null(markers)) {
      ids <- intersect(unlist(markers), rownames(x))
      m.ids = match(ids, rownames(x))
      Sigma.ct <- Sigma.ct[, m.ids]
    }
    if (verbose) {
      message("Creating Covariance Matrix...")
    }
  }
  else {
    Sigma <- sapply(unique(clusters), function(ct) {
      apply(sapply(unique(samples), function(sid) {
        y = exprs(x)[, clusters %in% ct & samples %in% 
                       sid, drop = FALSE]
        rowSums(y)/sum(y)
      }), 1, var, na.rm = TRUE)
    })
    if (!is.null(select.ct)) {
      m.ct = match(select.ct, colnames(Sigma))
      Sigma = Sigma[, m.ct]
    }
    if (!is.null(markers)) {
      ids <- intersect(unlist(markers), rownames(x))
      m.ids = match(ids, rownames(x))
      Sigma <- Sigma[m.ids, ]
    }
    if (verbose) {
      message("Creating Variance Matrix...")
    }
  }
  S <- sapply(unique(clusters), function(ct) {
    my.rowMeans(sapply(unique(samples), function(sid) {
      y = exprs(x)[, clusters %in% ct & samples %in% sid, 
                   drop = FALSE]
      sum(y)/ncol(y)
    }), na.rm = TRUE)
  })
  if (verbose) {
    message("Creating Library Size Matrix...")
  }
  S[S == 0] = NA
  M.S = colMeans(S, na.rm = TRUE)
  if(!is.null(my_abundance)){
    my_ms_names <- names(M.S)
    my_abundance <- my_abundance[my_ms_names %in% Celltype, ]
    M.S <- my_abundance[match(my_ms_names, my_abundance$Celltype),]
    M.S <- M.S$w_mean
    names(M.S) <- my_ms_names
  }
  D <- t(t(M.theta) * M.S)
  if (!is.null(select.ct)) {
    m.ct = match(select.ct, colnames(D))
    D = D[, m.ct]
    S = S[, m.ct]
    M.S = M.S[m.ct]
    M.theta = M.theta[, m.ct]
  }
  if (!is.null(markers)) {
    ids <- intersect(unlist(markers), rownames(x))
    m.ids = match(ids, rownames(x))
    D <- D[m.ids, ]
    M.theta <- M.theta[m.ids, ]
  }
  if (ct.cov) {
    return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta, 
                Sigma.ct = Sigma.ct))
  }
  else {
    return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta, 
                Sigma = Sigma))
  }
}





# prop --------------------------------------------------------------------

my_music_prop <- function (bulk.eset, sc.eset, markers = NULL, clusters, samples, 
                           select.ct = NULL, ct.cov = FALSE, verbose = TRUE, iter.max = 1000, 
                           nu = 1e-04, eps = 0.01, inter.as.bseq = F, normalize = F, my_abundance = NULL, 
                           ...) 
{
  bulk.gene = rownames(bulk.eset)[rowMeans(exprs(bulk.eset)) != 
                                    0]
  bulk.eset = bulk.eset[bulk.gene, , drop = FALSE]
  if (is.null(markers)) {
    sc.markers = bulk.gene
  }
  else {
    sc.markers = intersect(bulk.gene, unlist(markers))
  }
  sc.basis = my_music_basis(sc.eset, non.zero = TRUE, markers = sc.markers, 
                         clusters = clusters, samples = samples, select.ct = select.ct, 
                         ct.cov = ct.cov, verbose = verbose, my_abundance = my_abundance)
  cm.gene = intersect(rownames(sc.basis$Disgn.mtx), bulk.gene)
  if (is.null(markers)) {
    if (length(cm.gene) < 0.2 * min(length(bulk.gene), nrow(sc.eset))) 
      stop("Too few common genes!")
  }
  else {
    if (length(cm.gene) < 0.2 * length(unlist(markers))) 
      stop("Too few common genes!")
  }
  if (verbose) {
    message(paste("Used", length(cm.gene), "common genes..."))
  }
  m.sc = match(cm.gene, rownames(sc.basis$Disgn.mtx))
  m.bulk = match(cm.gene, bulk.gene)
  D1 = sc.basis$Disgn.mtx[m.sc, ]
  M.S = colMeans(sc.basis$S, na.rm = T)
  if(!is.null(my_abundance)){
    my_ms_names <- names(M.S)
    my_abundance <- my_abundance[my_ms_names %in% Celltype, ]
    M.S <- my_abundance[match(my_ms_names, my_abundance$Celltype),]
    M.S <- M.S$w_mean
    names(M.S) <- my_ms_names
  }
  Yjg = relative.ab(exprs(bulk.eset)[m.bulk, ])
  N.bulk = ncol(bulk.eset)
  if (ct.cov) {
    Sigma.ct = sc.basis$Sigma.ct[, m.sc]
    if (sum(Yjg[, i] == 0) > 0) {
      D1.temp = D1[Yjg[, i] != 0, ]
      Yjg.temp = Yjg[Yjg[, i] != 0, i]
      Sigma.ct.temp = Sigma.ct[, Yjg[, i] != 0]
      if (verbose) 
        message(paste(colnames(Yjg)[i], "has common genes", 
                      sum(Yjg[, i] != 0), "..."))
    }
    else {
      D1.temp = D1
      Yjg.temp = Yjg[, i]
      Sigma.ct.temp = Sigma.ct
      if (verbose) 
        message(paste(colnames(Yjg)[i], "has common genes", 
                      sum(Yjg[, i] != 0), "..."))
    }
    lm.D1.weighted = music.iter.ct(Yjg.temp, D1.temp, M.S, 
                                   Sigma.ct.temp, iter.max = iter.max, nu = nu, eps = eps, 
                                   inter.as.bseq = inter.as.bseq, normalize = normalize)
    Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
    Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
    weight.gene.temp = rep(NA, nrow(Yjg))
    weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
    Weight.gene = cbind(Weight.gene, weight.gene.temp)
    r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
    Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
  }
  else {
    Sigma = sc.basis$Sigma[m.sc, ]
    valid.ct = (colSums(is.na(Sigma)) == 0) & (colSums(is.na(D1)) == 
                                                 0) & (!is.na(M.S))
    if (sum(valid.ct) <= 1) {
      stop("Not enough valid cell type!")
    }
    if (verbose) {
      message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
    }
    D1 = D1[, valid.ct]
    M.S = M.S[valid.ct]
    Sigma = Sigma[, valid.ct]
    Est.prop.allgene = NULL
    Est.prop.weighted = NULL
    Weight.gene = NULL
    r.squared.full = NULL
    Var.prop = NULL
    for (i in 1:N.bulk) {
      if (sum(Yjg[, i] == 0) > 0) {
        D1.temp = D1[Yjg[, i] != 0, ]
        Yjg.temp = Yjg[Yjg[, i] != 0, i]
        Sigma.temp = Sigma[Yjg[, i] != 0, ]
        if (verbose) 
          message(paste(colnames(Yjg)[i], "has common genes", 
                        sum(Yjg[, i] != 0), "..."))
      }
      else {
        D1.temp = D1
        Yjg.temp = Yjg[, i]
        Sigma.temp = Sigma
        if (verbose) 
          message(paste(colnames(Yjg)[i], "has common genes", 
                        sum(Yjg[, i] != 0), "..."))
      }
      lm.D1.weighted = music.iter(Yjg.temp, D1.temp, M.S, 
                                  Sigma.temp, iter.max = iter.max, nu = nu, eps = eps, 
                                  inter.as.bseq = inter.as.bseq, normalize = normalize)
      Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
      Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
      weight.gene.temp = rep(NA, nrow(Yjg))
      weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
      Weight.gene = cbind(Weight.gene, weight.gene.temp)
      r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
      Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
    }
  }
  colnames(Est.prop.weighted) = colnames(D1)
  rownames(Est.prop.weighted) = colnames(Yjg)
  colnames(Est.prop.allgene) = colnames(D1)
  rownames(Est.prop.allgene) = colnames(Yjg)
  names(r.squared.full) = colnames(Yjg)
  colnames(Weight.gene) = colnames(Yjg)
  rownames(Weight.gene) = cm.gene
  colnames(Var.prop) = colnames(D1)
  rownames(Var.prop) = colnames(Yjg)
  return(list(Est.prop.weighted = Est.prop.weighted, Est.prop.allgene = Est.prop.allgene, 
              Weight.gene = Weight.gene, r.squared.full = r.squared.full, 
              Var.prop = Var.prop))
}
