vcov.lmerMod <- function(object, full = TRUE, ...) {
  if (!is(object, "lmerMod")) stop("estfun.lmerMod() only works for lmer() models.")
  if (object@devcomp$dims[10] != 0) stop("estfun.lmerMod() only works for ML estimation.")
  
  
  ## preparation for Block 4:
  ## get all elements by getME and exclude multiple random effect models.
  parts <- getME(object, "ALL")
  
  ## shortcut, y-X\beta
  yXbe <- parts$y - tcrossprod(parts$X, t(parts$beta))
  ## total variance
  Zlam <- tcrossprod(parts$Z, parts$Lambdat)
  V <- (tcrossprod(Zlam, Zlam) + diag(1, parts$n)) * (parts$sigma)^2
  M <- solve(chol(V))
  invV <- tcrossprod(M, M)
    
  ## block 1 fixhes 
  if (full == FALSE) {
    fixvar <- solve(tcrossprod(crossprod(parts$X, invV), t(parts$X)))
    fixvar  
  } else {
    fixhes <- tcrossprod(crossprod(parts$X, invV), t(parts$X))
  
    ## length = (variance covariance parameter in G + residual variance.)
    uluti <- length(parts$theta)
    devV <- vector("list", (uluti + 1))
  
    ## get d(G)/d(sigma), faciliated by d(Lambda). 
    devLambda <- vector("list", uluti)
    score_varcov <- matrix(NA, nrow = length(parts$y), ncol = uluti)
    for (i in 1:uluti) {
      devLambda[[i]] <- parts$Lambda
      devLambda[[i]][which(devLambda[[i]] != parts$theta[i])] <- 0
      devLambda[[i]][which(parts$Lambda == parts$theta[i])] <- 1
      devLambda[[i]] <- forceSymmetric(devLambda[[i]], uplo = "L")
      devV[[i]] <- tcrossprod(tcrossprod(parts$Z, t(devLambda[[i]])), parts$Z)
    }
    devV[[(uluti+1)]] <- diag(1, nrow(parts$X))
  
    ## Block 4: sigma's second derivative.
    ranhes <- matrix(NA, nrow = (uluti + 1), ncol = (uluti + 1))
    ## combinations (allow repeat) to indicate entries
    entries <- rbind(matrix(rep(1: (uluti + 1), each = 2),
      (uluti + 1), 2, byrow = TRUE), t(combn((uluti + 1), 2)))
      entries <- entries[order(entries[,1], entries[,2]), ]
  
    for (i in 1: nrow(entries)) {
      ranhes[lower.tri(ranhes, diag = TRUE)][i] <- as.numeric((1/2) *
        sum(diag(tcrossprod(tcrossprod(crossprod(invV, devV[[entries[i,1]]]), 
        invV), t(devV[[entries[i,2]]])))))
    }
    ranhes <- forceSymmetric(ranhes, uplo = "L")
  
    ## block 2 and block 3: second derivative of sigma and beta.
    varcov_beta <- matrix(NA, length(devV), length(parts$beta))
    for (j in 1 : (length(devV))) {
      varcov_beta[j,] <- as.vector(tcrossprod(crossprod(parts$X, 
        (tcrossprod(crossprod(invV, devV[[j]]), invV))), t(yXbe)))
    }
  
    ## Organize full_varcov
    full_varcov <- solve(rbind(cbind(fixhes, t(varcov_beta)),
      cbind(varcov_beta, ranhes)))
  
    colnames(full_varcov) <- c(names(parts$fixef), paste("cov",
      names(parts$theta), sep="_"), "residual")

    return(full_varcov)
  }
}
