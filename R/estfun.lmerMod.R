estfun.lmerMod <- function(x, ...) {
  if (!is(x, "lmerMod")) stop("estfun.lmerMod() only works for lmer() models.")
  if (!is.null(x@call$weights)) stop ("Models with weights specification is currently not supported.")
  if (length(grep("cbind", x@call$formula))!=0) stop ("Models with cbind specification is currently not supported.")
  
  ## warning for high correlations. 
  cor <- attr(lme4::VarCorr(x)[[1]], "correlation")
  if(any(abs(cor[lower.tri(cor)]) > .9)) warning("Correlations > |.9| detected. Scores of random (co)variances may be unstable.")
  
  ## error for variance equal to 0 or non-positive definite
  ndim = nrow(cor)
  if (ndim == 1 & unclass(VarCorr(x))[[1]]  == 0) stop ("Random effect's variance is close to 0. Cannot compute derivatives.")
  if (ndim > 1 & matrixcalc::is.positive.definite(unclass(VarCorr(x))[[1]])) stop ("Random effect's variance covariance matrix is not positive definite. Cannot compute derivatives.")
  
  ## get all elements by getME and exclude multiple random effect models.
  parts <- getME(x, "ALL")
  
  ## obtain dot arguments
  dotdotdot <- list(...)
  if("level" %in% names(dotdotdot)){
    level <- dotdotdot$level
  } else {
    if(length(parts$l_i) > 1L){
      ## multiple cluster variables, so need level=1
      level <- 1L
    } else {
      level <- 2L
    }
  }
  
  if("ranpar" %in% names(dotdotdot)){
    ranpar <- dotdotdot$ranpar
  } else {
    ranpar <- "var"
  }  

  ## checks
  if(!(level %in% c(1L, 2L))) stop("invalid 'level' argument supplied")
  if (length(parts$l_i) > 1L & level == 2L) stop("Multiple cluster variables detected. Supply 'level=1' argument to estfun.lmerMod().")
  
  ## prepare shortcuts
  uluti <- length(parts$theta)
  yXbe <- parts$y - tcrossprod(parts$X, t(parts$beta))
  Zlam <- tcrossprod(parts$Z, parts$Lambdat)
  V <- (tcrossprod(Zlam, Zlam) + diag(1, parts$n)) * (parts$sigma)^2
  M <- solve(chol(V))
  invV <- tcrossprod(M, M)
  yXbesoV <- crossprod(yXbe, invV)
  LambdaInd <- parts$Lambda
  LambdaInd@x[] <- seq(1:uluti)
  invVX <- crossprod(parts$X, invV)
  Pmid <- solve(crossprod(parts$X, t(invVX)))
  P <- invV - tcrossprod(crossprod(invVX, Pmid), t(invVX))
  
  ## adapt from Stroup book page 131, last eq,
  ## score for fixed effects parameter: score_beta=XR^{-1}(y-Zb-X\beta)
  score_beta <- t(crossprod(parts$X, invV)) * (yXbe)
  
  ## prepare for score of variance covariance parameters. 
  ## get d(G)/d(sigma), faciliated by d(Lambda).
  devLambda <- vector("list", uluti)
  score_varcov <- matrix(NA, nrow = length(parts$y), ncol = uluti)
  devV <- vector ("list", (uluti + 1))
  
  for (i in 1:uluti) {
    devLambda[[i]] <- forceSymmetric(LambdaInd==i, uplo = "L")
    devV[[i]] <- tcrossprod(tcrossprod(parts$Z, t(devLambda[[i]])), parts$Z)
  }
  devV[[(uluti+1)]] <- diag(1,nrow(parts$X))
  
  ## score for variance covariance parameter
  score_varcov <- matrix(NA, nrow = nrow(parts$X), ncol = (uluti + 1))
  
  ## ML estimates. 
  if (parts$devcomp$dims[['REML']] == 0) {
    for (j in 1:length(devV)) {
      score_varcov[,j] <- as.vector(-(1/2) * diag(crossprod(invV, devV[[j]])) +
      t((1/2) * tcrossprod(tcrossprod(yXbesoV, t(devV[[j]])), invV)) *
      (yXbe))
    }
  }

  ## REML estimates
  if (parts$devcomp$dims[['REML']] > 0) {
    for (j in 1:length(devV)) {
      score_varcov[,j] <- as.vector(-(1/2) * diag(crossprod(P, devV[[j]])) +
      t((1/2) * tcrossprod(tcrossprod(yXbesoV, t(devV[[j]])), invV)) *
      (yXbe))
    }
  }

  ## re-parameterization for sd/cor
  if (ranpar == "var"){
    score_varcov <- score_varcov 
  } else if (ranpar == "sd"){
      sdcormat <- as.data.frame(VarCorr(x,comp = "Std.Dev"),
                                order = "lower.tri")
      sdcormat$sdcor2[which(is.na(sdcormat$var2))] <-
        sdcormat$sdcor[which(is.na(sdcormat$var2))]*2
      sdcormat$sdcor2[which(!is.na(sdcormat$var2))] <-
        sdcormat$vcov[which(!is.na(sdcormat$var2))]/
          sdcormat$sdcor[which(!is.na(sdcormat$var2))]
      score_varcov <- sweep(score_varcov, MARGIN = 2, sdcormat$sdcor2, `*`)
    } else {
        stop("ranpar needs to be var or sd for lmerMod object.")
  }
  
  
  ## Organize score matrix
  score <- cbind(as.matrix(score_beta), as.matrix(score_varcov))
  colnames(score) <- c(names(parts$fixef),
                       paste("cov", names(parts$theta), sep="_"), "residual")

  ## Clusterwise scores if level==2
  if (level == 2) {
    #index <- rep(1:parts$l_i, (parts$n/parts$l_i))
    index <- parts$flist[[1]]
    #index <- index[order(index)]
    scoretemp <- aggregate(x = score, by = list(index), FUN = sum)
    score <- as.matrix(scoretemp[,-1])
    rownames(score) <- scoretemp[,1]
  }
  
  return(score)
}
