setwd("~/Library/CloudStorage/GoogleDrive-gianluca.sottile@community.unipa.it/Drive condivisi/Fici_Rita/04 - R/Simulazioni Rita/functions")
source(file = "rStar.R")
source(file = "rStar.l.R")
source(file = "rPar.R")
source(file = "pfpca_tmp.R")
source(file = "prcurve2.R")
source(file = "rcfggm.R")
source(file = "KL function.R")
source(file = "KLCV.R")
source(file = "AIC eBIC.R")
source(file = "lambda and rho iteration functions.R")

source("~/Dropbox/group-cglasso/08 - github/Hematopoiesis-network-inference-from-RT-qPCR-data-main/01 - RCode/datajcggm.R")
source("~/Dropbox/group-cglasso/08 - github/Hematopoiesis-network-inference-from-RT-qPCR-data-main/01 - RCode/jcglasso.R")
source("~/Dropbox/group-cglasso/08 - github/Hematopoiesis-network-inference-from-RT-qPCR-data-main/01 - RCode/jcggm.R")
source("~/Dropbox/group-cglasso/08 - github/Hematopoiesis-network-inference-from-RT-qPCR-data-main/01 - RCode/gof.R")
source("~/Dropbox/group-cglasso/08 - github/Hematopoiesis-network-inference-from-RT-qPCR-data-main/01 - RCode/jcglasso.R")
source("~/Dropbox/group-cglasso/08 - github/Hematopoiesis-network-inference-from-RT-qPCR-data-main/01 - RCode/to_graph.R")

#####################

lmbMax_fun <- function(obj, lambda.id = 1L, rho.id = 1L) {
  K <- length(obj$Z)
  seq_k <- seq_len(K)
  n <- obj$nobs
  q <- obj$npred
  p <- obj$nresp
  fk <- proportions(n)
  
  id_Y <- obj$InfoStructure$id_Y
  id_X <- obj$InfoStructure$id_X
  
  lambda_max <- max(sqrt(rowSums(sapply(seq_k, function(k) {
    X <- obj$Zipt[seq_len(n[k]), id_X, k, lambda.id, rho.id]
    Y <- obj$Zipt[seq_len(n[k]), id_Y, k, lambda.id, rho.id]
    B <- obj$B[, , k, lambda.id, rho.id]
    Ycenter <- Y - cbind(1, X) %*% B
    Sxy <- crossprod(X, Ycenter) / n[k]
    Thtyy <- obj$Tht[id_Y, id_Y, k, lambda.id, rho.id]
    # max(fk[k] * abs(Sxy %*% Thtyy / matrix(diag(Thtyy), nrow = q, ncol = p)))
    pmax(fk[k] * abs(Sxy %*% Thtyy), 0)})^2)))
  lambda_max
}
rhoMax_fun <- function(obj, lambda.id = 1L, rho.id = 1L) {
  K <- length(obj$Z)
  seq_k <- seq_len(K)
  n <- obj$nobs
  fk <- proportions(n)
  p <- obj$nresp
  U <- outer(1:p, 1:p, "<")
  
  id_Y <- obj$InfoStructure$id_Y
  id_X <- obj$InfoStructure$id_X
  
  rho_max <- max(sqrt(rowSums(sapply(seq_k, function(k) {
    X <- obj$Zipt[seq_len(n[k]), id_X, k, lambda.id, rho.id]
    Y <- obj$Zipt[seq_len(n[k]), id_Y, k, lambda.id, rho.id]
    B <- obj$B[, , k, lambda.id, rho.id]
    Ycenter <- Y - cbind(1, X) %*% B
    Syy <- crossprod(Ycenter) / n[k]
    pmax(fk[k] * abs(Syy)[U], 0)})^2)))
  
  rho_max
}
ROCTheta <- function(thetah, E, nrho, plot.it = TRUE){
  p <- dim(thetah)[1]
  Tht.sup <- matrix(0, p, p)
  Tht.sup[E] <- 1
  U <- lower.tri(Tht.sup, diag = FALSE)
  thetat_v <- Tht.sup[U]
  A <- which(abs(thetat_v) > 0)
  thetah_m <- apply(thetah, 3, function(M) M[U])
  
  temp <- sapply(1:nrho, \(r){
    tbl <- table(factor(thetat_v, levels = c(1, 0)), factor(1*(thetah_m[, r] != 0), levels = c(1, 0)))  
    c(tbl[1,1]/sum(tbl[1,]), tbl[2,1]/sum(tbl[2,]), sum(diag(tbl)) / sum(tbl))
  })
  TPR <- temp[1,]
  FPR <- temp[2,]
  Accuracy <- temp[3,]
  
  tblMin <- table(factor(thetat_v, levels = c(1, 0)), factor(1*(thetah_m[, 1]*0 != 0), levels = c(1, 0)))  
  tblMax <- table(factor(thetat_v, levels = c(1, 0)), factor(1*((thetah_m[, 1]+1) != 0), levels = c(1, 0)))  
  
  TPR.roc <- c(tblMin[1,1]/sum(tblMin[1,]), TPR, tblMax[1,1]/sum(tblMax[1,]))
  FPR.roc <- c(tblMin[2,1]/sum(tblMin[2,]), FPR, tblMax[2,1]/sum(tblMax[2,]))
  
  splFun <- suppressWarnings(splinefun(FPR.roc, TPR.roc, method = "monoH.FC"))
  if(plot.it) {
    plot(FPR.roc, TPR.roc, type = "b", xlim = c(0, 1), ylim = c(0, 1))
    curve(splFun, from = 0, to = 1, add = TRUE, col = 2, lty = 2)
    abline(0, 1, col = 3)
  }
  auc <- integrate(splFun, lower = 0, upper = 1)$value
  
  out <- list(TPR = TPR, FPR = FPR, ROC = auc, accuracy = Accuracy)
  out
}
ROCB <- function(Bh, Bt, nlambda, plot.it = TRUE){
  Bt_v <- c(Bt)
  A <- which(abs(Bt_v) > 0)
  Bh_m <- apply(Bh, 3, function(M) c(M))
  
  temp <- sapply(1:nlambda, \(l){
    tbl <- table(factor(Bt_v, levels = c(1, 0)), factor(1*(Bh_m[, l] != 0), levels = c(1, 0)))  
    c(tbl[1,1]/sum(tbl[1,]), tbl[2,1]/sum(tbl[2,]), sum(diag(tbl)) / sum(tbl))
  })
  TPR <- temp[1,]
  FPR <- temp[2,]
  Accuracy <- temp[3,]
  
  tblMin <- table(factor(Bt_v, levels = c(1, 0)), factor(1*(Bh_m[, 1]*0 != 0), levels = c(1, 0)))  
  tblMax <- table(factor(Bt_v, levels = c(1, 0)), factor(1*((Bh_m[, 1]+1) != 0), levels = c(1, 0)))  
  
  TPR.roc <- c(tblMin[1,1]/sum(tblMin[1,]), TPR, tblMax[1,1]/sum(tblMax[1,]))
  FPR.roc <- c(tblMin[2,1]/sum(tblMin[2,]), FPR, tblMax[2,1]/sum(tblMax[2,]))
  
  splFun <- suppressWarnings(splinefun(FPR.roc, TPR.roc, method = "monoH.FC"))
  if(plot.it) {
    plot(FPR.roc, TPR.roc, type = "b", xlim = c(0, 1), ylim = c(0, 1))
    curve(splFun, from = 0, to = 1, add = TRUE, col = 2, lty = 2)
    abline(0, 1, col = 3)
  }
  auc <- integrate(splFun, lower = 0, upper = 1)$value
  
  out <- list(TPR = TPR, FPR = FPR, ROC = auc, accuracy = Accuracy)
  out
}

#####################

# library
library(cglasso)
library(fda)
library(grplasso)
library(JGL)
library(ggplot2)

### setting ####
n <- 50 # 100
L <- 3
nstars.size <- 5
nstars <- 5
p <- nstars.size * nstars * (L + 1)
q <- 10
n.ThtY <- p * ( p - 1 ) / 2
n.b <- p * q

ntp <- tmp <- 30

a <- .4
B.min <- a * 2.5 # 1.00
B.max <- a * 3.5 # 1.40

b <- .2
tht.min <- b * 2   # 0.40
tht.max <- b * 2.5 # 0.50

g.ebic <- 0.5

Omg <- diag(q) * 1
phi <- create.fourier.basis(nbasis = L)
scal.terms <- c(1, 0.2)
scal.terms[1] * c(1:L)^{-scal.terms[2]}

nsim <- 100

rho.min.ratio <- 0.1
rho.max.ratio <- 1
nrho <- 10
perc.rho.seq <- seq(from = rho.max.ratio , to = rho.min.ratio , length.out = nrho)
lmbd.min.ratio <- 0.1
lmbd.max.ratio <- 1
nlambda <- 10
perc.lmbd.seq <- seq(from = lmbd.max.ratio , to = lmbd.min.ratio , length.out = nlambda)

B.supp <- matrix(0, nrow = q, ncol = p,
                 dimnames = list(paste0("X", seq_len(q)),
                                 paste0("Y", seq_len(p))))
set.seed(1)
for(k in seq_len(p)) B.supp[, k] <- sample(c(rep(0, q - 2), rep(1, 2)))
Par <- rPar(p = p, q = q, L = L, B.supp = B.supp, B.min = B.min, 
            B.max = B.max, nstars.size = nstars.size, nstars = nstars,  
            tht.min = tht.min, tht.max = tht.max, Omg = Omg, scal.terms = scal.terms, union = FALSE)
Sigma_trueY <- array(NA, dim = c(p,p,L))
B.true <- vector( mode = "list", length = L)
for( l in seq_len(L)){
  Sigma_trueY[, , l] <- Par$Sgm[-(1:q), -(1:q), l]
  sig.x <- Par$Sgm[1:q, 1:q, l]
  B.true[[l]] <-  solve(sig.x) %*% Par$Sgm[1:q, -c(1:q), l]
  print(diag(Par$Sgm[-c(1:q), -c(1:q), l] - t(B.true[[l]]) %*% sig.x %*% B.true[[l]]))
}
Tht.supp <- matrix(0, p, p)
for (l in seq_len(L)) Tht.supp <- Tht.supp + (Par$Tht[-seq_len(q), -seq_len(q), l] != 0)

par(mfrow = c(1, 2))
image(Tht.supp != 0)
image(B.supp)

pdf("~/Downloads/structure_simul.pdf", width = 16*.75, height = 8*.75)
par(mfrow=c(1, 2))
image(1 * (Tht.supp != 0), col = c("white", "black"), 
      ylab = "", cex.lab = 2, cex.main = 3, 
      main = expression(Theta[l]), axes = FALSE)
axis(1, at = c(0, .2, .4, .6, .8, 1), labels = c(0, 10, 20, 30, 40, 50)*2, cex.axis = 2)
axis(2, at = c(0, .2, .4, .6, .8, 1), labels = c(0, 10, 20, 30, 40, 50)*2, cex.axis = 2)
image(B.supp, col = c("white", "black"), 
      ylab = "", cex.lab = 2, cex.main = 3, 
      main = expression(B[l]), axes = FALSE)
axis(1, at = c(0, .2, .4, .6, .8, 1), labels = c(0, .2, .4, .6, .8, 1)*10, cex.axis = 2)
axis(2, at = c(0, .2, .4, .6, .8, 1), labels = c(0, 10, 20, 30, 40, 50)*2, cex.axis = 2)
dev.off()

Tht.supp[upper.tri(Tht.supp, TRUE)] <- 0
E.tht <- which(Tht.supp != 0, arr.ind = TRUE)

range(list2array(B.true))
range(Par$Tht[-c(1:q), -c(1:q), ])

###### Tht and Sgm [Y X]
Sgm.YX <- Tht.YX <- array(NA , dim = dim(Par$Tht))
colnames(Sgm.YX) <- colnames(Tht.YX) <- rownames(Tht.YX) <- rownames(Sgm.YX) <- c(paste0("Y", seq_len(p)), paste0("X", seq_len(q)))

for (l in seq_len(L)){
  Sgm.YX[ 1:p, 1:p, l] <- Par$Sgm[-c(1:q), -c(1:q), l]
  Sgm.YX[ -c(1:p), -c(1:p), l] <- Par$Sgm[1:q, 1:q, l]
  Sgm.YX[ 1:p, -c(1:p), l] <- Par$Sgm[-c(1:q), 1:q, l]
  Sgm.YX[ -c(1:p), 1:p, l] <- Par$Sgm[1:q, -c(1:q), l]
  Tht.YX[ , , l] <- round(solve(Sgm.YX[ , , l]), 5)
}

#####################

resAUC <- array(NA, dim = c(nsim, nlambda))
resMSE <- array(NA, dim = c(nsim, nlambda, nrho))
resKL <- array(NA, dim = c(nsim, nlambda, nrho))
resMIN <- array(NA, dim = c(nsim, nlambda, 2, 3))

resAUC2 <- array(NA, dim = c(nsim))
resMSE2 <- array(NA, dim = c(nsim, nrho))

resAUC3 <- array(NA, dim = c(nsim, nlambda))
resMSE3 <- array(NA, dim = c(nsim, nlambda, nrho))

#####################

pb <- txtProgressBar(min = 0L, max = nsim, style = 3L)
jj <- 0L
set.seed(1234)
for(h in seq_len(nsim)) {
  jj <- jj + 1L
  setTxtProgressBar(pb, jj)
  df <- rcfggm(n = n, p = p, q = q, ntp = ntp, Sgm = Par$Sgm, phi = phi)$Xi.z
  
  x <- y <-  df1 <- Sgm <- data.list <- eps.Y <- B.stim <- out <- vector(mode = "list", length = L)
  for( l in seq_len(L)) {
    x[[l]] <- df[, 1:q, l]
    y[[l]] <- df[, -c(1:q), l]
  }
  data.listY <- datajcggm(Y = y)
  data.listYX <- datajcggm(Y = y, X = x)
  
  #### JCGLASSO
  
  outMAX <- jcglasso(data = data.listYX, nrho = 1L, nlambda = 1L,
                     alpha1 = 0.0, alpha2 = 0.0, alpha3 = 0.0, nu = 0.0)
  lambda.seq <- lmbMax_fun(outMAX) * perc.lmbd.seq
  rho.seq <- rhoMax_fun(outMAX) * perc.rho.seq
  
  for(j in 1:nlambda) {
    out.lambda <- jcglasso(data = data.listYX, rho = 1e12, lambda = lambda.seq[j],
                           alpha1 = 0.0, alpha2 = 0.0, alpha3 = 0.0, nu = 0.0)
    rho.seq <- rhoMax_fun(out.lambda) * perc.rho.seq
    out.tht <- jcglasso(data = data.listYX, rho = rho.seq, lambda = lambda.seq[j],
                        alpha1 = 0.0, alpha2 = 0.0, alpha3 = 0.0, nu = 0.0, trace = 0L)
    
    thetah <- 0.0
    b.stim <- tht.stimYX <- tht.stim <- sig.stim <- B.stim <- vector(mode = "list", length = L)
    Bmatrix <- array(0.0, dim = c(q + 1, p, L, nrho))
    THTmatrix <- array(0.0, dim = c(p, p, L, nrho))
    for(l in seq_len(L)) {
      Bmatrix[,,l,] <- coef(out.tht, type = "B", class.id = l)
      b.stim[[l]] <- coef(out.tht, type = "B", class.id = l)[-1, , ]
      B.stim[[l]] <- aperm(b.stim[[l]], c(2, 1, 3))
      tht.stimYX[[l]] <- out.tht$Tht[, , l, , ]
      THTmatrix[,,l,] <- coef(out.tht, type = "Theta", class.id = l)
      tht.stim[[l]] <- coef(out.tht, type = "Theta", class.id = l)
      thetah <- thetah + tht.stim[[l]]
      sig.stim[[l]] <- out.tht$Sgm[out.tht$InfoStructure$id_Y, out.tht$InfoStructure$id_Y, l, 1, ]
    }
    tmp <- ROCTheta(thetah = thetah, E = E.tht, nrho = nrho, plot.it = FALSE)
    # print(tmp)
    tmp2 <- sapply(1:nrho, \(r) mean(sapply(1:L, \(l) norm(Par$Tht[-c(1:q), -c(1:q), l] - THTmatrix[, , l, r], "F"))))
    
    KL_loss <- KL.tht1(Theta_stim = tht.stim, Theta_true = Par$Tht[-c(1:q), -c(1:q), ], L = L,
                       npar = nrho, X = x, B_stim = B.stim, B_true = Par$B)
    # par(mfrow = c(1, 1))
    # plot(rho.seq, KL_loss, type = "l")
    # KL_loss2 <- KL.tht_B(ThetaXY_stim = tht.stimYX, ThetaXY_true = Par$Tht, L = L, npar = nrho, p = p, q = q)
    # plot(rho.seq, KL_loss2, type = "l")
    
    outKLCV <- KLCV.gauss2(tht.stimYX = tht.stimYX, b.stim = b.stim,
                           data = data.listYX, L = L, whole.tht = FALSE)
    # outKLCV2 <- KLCV.fixed.X(tht.stim = tht.stim, eps.Y = eps.Y, L)
    
    outAIC <- AIC(out.tht)$value_gof
    outEBIC <- BIC(out.tht, g = g.ebic, type = "FD")$value_gof
    # par(mfrow = c(1, 3))
    # plot(rho.seq, outKLCV, type = "l")
    # plot(rho.seq, outAIC, type = "l")
    # plot(rho.seq, outEBIC, type = "l")
    
    id.min.klcv <- which.min(outKLCV)
    id.min.aic <- which.min(outAIC)
    id.min.ebic <- which.min(outEBIC)
    
    resAUC[jj, j] <- tmp$ROC
    resMSE[jj, j, ] <- tmp2
    resKL[jj, j, ] <- KL_loss
    resMIN[jj, j, 1, ] <- KL_loss[c(id.min.klcv, id.min.aic, id.min.ebic)]
    resMIN[jj, j, 2, ] <- tmp$accuracy[c(id.min.klcv, id.min.aic, id.min.ebic)]
  }
  
  cat("\n\n")
  print(cbind(AUC = resAUC[jj, ], resMSE[jj, , ]))
  cat("\n")
  cbind(resMIN[jj, , 1, ], NA, resMIN[jj, , 2, ])
  
  #### JGL
  
  out.JGL.max <- jcglasso(data = data.listY, rho = 1E12, alpha1 = 0.0)
  rho.seq <- rhoMax_fun(out.JGL.max) * perc.rho.seq
  outJGL <- jcglasso(data = data.listY, rho = rho.seq, alpha1 = 0.0, trace = 0L)
  
  thetah <- 0.0
  tht.stimYX <- tht.stim <- sig.stim <- vector(mode = "list", length = L)
  Bmatrix <- array(0.0, dim = c(q + 1, p, L, nrho))
  THTmatrix <- array(0.0, dim = c(p, p, L, nrho))
  for(l in seq_len(L)) {
    Bmatrix[1L, , l, ] <- coef(outJGL, type = "B", class.id = l)
    THTmatrix[, , l, ] <- coef(outJGL, type = "Theta", class.id = l)
    tht.stimYX[[l]] <- outJGL$Tht[, , l, , ]
    tht.stim[[l]] <- coef(outJGL, type = "Theta", class.id = l)
    thetah <- thetah + tht.stim[[l]]
    sig.stim[[l]] <- outJGL$Sgm[outJGL$InfoStructure$id_Y, outJGL$InfoStructure$id_Y, l, 1, ]
  }
  tmp <- ROCTheta(thetah = thetah, E = E.tht, nrho = nrho, plot.it = FALSE)
  tmp2 <- sapply(1:nrho, \(r) mean(sapply(1:L, \(l) norm(Par$Tht[-c(1:q), -c(1:q), l] - THTmatrix[, , l, r], "F"))))
  
  resAUC2[jj] <- tmp$ROC
  resMSE2[jj, ] <- tmp2
  
  cat("\n\n")
  print(c(AUC = resAUC2[jj], resMSE2[jj, ]))
  
  #### CGLASSO
  
  b.stim <- array(0.0, dim = c(q, p, nrho, L, nlambda))
  B.stim <- array(0.0, dim = c(p, q, nrho, L, nlambda))
  thetah <- array(0.0, dim = c(p, p, nrho, L, nlambda))
  tht.stimYX <- array(0.0, dim = c(p + q, p + q, nrho, L, nlambda))
  sig.stim <- tht.stim <- array(0.0, dim = c(p, p, nrho, L, nlambda))
  Bmatrix <- array(0.0, dim = c(q + 1, p, L, nrho, nlambda))
  THTmatrix <- array(0.0, dim = c(p, p, L, nrho, nlambda))
  for(k in 1:L) {
    data.listYXk <- datajcggm(Y = y[k], X = x[k])
    outCGL <- jcglasso(data = data.listYXk, nrho = 1L, nlambda = 1L)
    lambda.seq <- lmbMax_fun(outCGL) * perc.lmbd.seq
    for(j in 1:nlambda){
      outCGL2 <- jcglasso(data = data.listYXk, rho = 1E12, lambda = lambda.seq[j], nu = 0.0, 
                          alpha1 = 0.0, alpha2 = 0.0, alpha3 = 0.0)  
      rho.seq <- rhoMax_fun(outCGL2) * perc.rho.seq
      outCGL3 <- jcglasso(data = data.listYXk, rho = rho.seq, lambda = lambda.seq[j], nu = 0.0, 
                          alpha1 = 0.0, alpha2 = 0.0, alpha3 = 0.0, trace = 0L)  
      
      Bmatrix[, , k, , j] <- coef(outCGL3, type = "B")
      b.stim[, , , k, j] <- coef(outCGL3, type = "B")[-1, , ]
      B.stim[, , , k, j] <- aperm(coef(outCGL3, type = "B")[-1, , ], c(2, 1, 3))
      THTmatrix[, , k, , j] <- coef(outCGL3, type = "Theta")
      thetah[, , , k, j] <- coef(outCGL3, type = "Theta")
      tht.stimYX[, , , k, j] <- outCGL3$Tht[, , 1L, 1L, ]
      tht.stim[, , , k, j] <- coef(outCGL3, type = "Theta")
      sig.stim[, , , k, j] <- outCGL3$Sgm[outCGL3$InfoStructure$id_Y, outCGL3$InfoStructure$id_Y, 1L, 1L, ]
    }
  }
  
  tmp <- sapply(1:nlambda, \(j) {
    thth <- 1 * ((thetah[, , , 1, j] != 0) | (thetah[, , , 2, j] != 0) | (thetah[, , , 3, j] != 0))
    ROCTheta(thetah = thth, E = E.tht, nrho = length(rho.seq), plot.it = FALSE)$ROC
  })
  
  tmp2 <- sapply(1:nlambda, \(j) {
    sapply(1:nrho, \(r) mean(sapply(1:L, \(l) norm(Par$Tht[-c(1:q), -c(1:q), l] - THTmatrix[, , l, r, j], "F"))))
  })
  
  # sapply(1:nlambda, \(j) {
  #   thth <- 1 * ((thetah3[, , , 1, j] != 0) & (thetah3[, , , 2, j] != 0) & (thetah3[, , , 3, j] != 0))
  #   ROCTheta(thetah = thth, E = E.tht, nrho = length(rho.seq), plot.it = FALSE)$ROC
  # })
  
  resAUC3[jj, ] <- tmp
  resMSE3[jj, , ] <- t(tmp2)
  
  cat("\n\n")
  print(cbind(AUC = resAUC3[jj, ], resMSE3[jj, , ]))
}

rbind(colMeans(resAUC), mean(resAUC2), colMeans(resAUC3))
rbind(apply(apply(resMSE, 3, colMeans), 2, min),
      min(colMeans(resMSE2)), apply(apply(resMSE3, 3, colMeans), 2, min))

apply(resMIN[, , 1, ], 2, colMeans)
apply(resMIN[, , 2, ], 2, colMeans)


dfAUC <- data.frame(value = c(colMeans(resAUC), rep(mean(resAUC2), nlambda), colMeans(resAUC3)), 
                    nu = rep(round(perc.lmbd.seq, 2), 3L),
                    method = factor(rep(c("fGGRM", "Zapata", "Naive"), each = nlambda), 
                                    levels = c("fGGRM", "Zapata", "Naive")))

p1 <- ggplot(dfAUC) + 
  geom_point(aes(x = nu, y = value, shape = method), size = 3) +
  geom_line(aes(x = nu, y = value, linetype = method)) +
  theme_classic() + ylab(expression(Area~under~the~curve~of~hat(Theta))) + xlab(expression(nu / nu[max])) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = .4), 
        legend.text = element_text(size = 14), strip.text = element_text(size = 12),
        legend.position = "inside", legend.title = element_blank(), legend.position.inside = c(.8,.8))
p1

# dfAUC2 <- data.frame(value = c(c(pmin(resAUC, 1.000)), rep(c(resAUC2), nlambda), c(resAUC3)), 
#                      lambda = factor(rep(rep(round(perc.lmbd.seq, 2), each = nsim), 3L), 
#                                      levels = round(rev(perc.lmbd.seq), 2)),
#                      method = factor(rep(c("fGGRM", "Zapata", "Naive"), each = nlambda*nsim), 
#                                      levels = c("fGGRM", "Zapata", "Naive")))
# 
# p1bis <- ggplot(dfAUC2) + 
#   geom_boxplot(aes(x = lambda, y = value, colour = method), outliers = FALSE, notch = TRUE) + 
#   theme_classic() + ylab(expression(Area~under~the~curve~of~hat(Theta))) + xlab(expression(lambda / lambda[max])) + 
#   theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12),
#         axis.text.x = element_text(angle = 45, vjust = .4), legend.text = element_text(size = 14),
#         legend.position = "inside", legend.title = element_blank(), legend.position.inside = c(.8,.8))
# p1bis

dfMSE2 <- data.frame(value = c(c(apply(resMSE, 3, colMeans)), rep(colMeans(resMSE2), nrho), c(apply(resMSE3, 3, colMeans))), 
                     nu = factor(rep(rep(round(perc.lmbd.seq, 2), each = nrho), 3L), 
                                 levels = round(rev(perc.lmbd.seq), 2),
                                 labels = c(expression(nu / nu[max] == 0.1),
                                            expression(nu / nu[max] == 0.2),
                                            expression(nu / nu[max] == 0.3),
                                            expression(nu / nu[max] == 0.4),
                                            expression(nu / nu[max] == 0.5),
                                            expression(nu / nu[max] == 0.6),
                                            expression(nu / nu[max] == 0.7),
                                            expression(nu / nu[max] == 0.8),
                                            expression(nu / nu[max] == 0.9),
                                            expression(nu / nu[max] == 1.0))),
                     rho = rep(rep(round(perc.rho.seq, 2), nrho), 3L),
                     method = factor(rep(c("fGGRM", "Zapata", "Naive"), each = nlambda*nrho), 
                                     levels = c("fGGRM", "Zapata", "Naive")))

p2bis <- ggplot(subset(dfMSE2, nu %in% c(expression(nu / nu[max] == 0.2),
                                         expression(nu / nu[max] == 0.4),
                                         expression(nu / nu[max] == 0.6),
                                         expression(nu / nu[max] == 0.8)))) + 
  geom_point(aes(x = rho, y = value, shape = method), size = 3) +
  geom_line(aes(x = rho, y = value, linetype = method)) +
  facet_wrap(nu ~ ., labeller = labeller(nu = label_parsed)) + 
  theme_classic() + ylab(expression(Mean~Squared~Error~of~hat(Theta))) + xlab(expression(rho / rho[max])) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = .4), 
        legend.text = element_text(size = 14), strip.text = element_text(size = 12),
        legend.position = "right", legend.title = element_blank())
p2bis

ggsave("~/Downloads/AUC_theta_n100.pdf", plot = p1, device = "pdf", scale = .7, width = 12, height = 10, units = "in")
ggsave("~/Downloads/MSE_theta_n100.pdf", plot = p2bis, device = "pdf", scale = .7, width = 12, height = 10, units = "in")

gridExtra::grid.arrange(p1, p2bis, nrow = 1L)
# gridExtra::grid.arrange(p1bis, p2bis, nrow = 1L)

