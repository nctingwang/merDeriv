vcov.lmerMod <- function(object, ranpar = "var", ...) {
  
  if (!is(object, "lmerMod")) stop("estfun.lmerMod() only works for lmer() models.")
  
  dotdotdot <- list(...)
  if("full" %in% names(dotdotdot)){
    full <- dotdotdot$full
  } else {
    full <- FALSE
  }
  if ("information" %in% names(dotdotdot)) {
    information <- dotdotdot$information
  } else {
    information <- "expected"
  }
  if(!(full %in% c("TRUE", "FALSE"))) stop("invalid 'full' argument supplied")
  if(!(information %in% c("expected", "observed"))) stop("invalid 'information' argument supplied")
  
  ## preparation for short cuts:
  ## get all elements by getME and exclude multiple random effect models.
  parts <- getME(object, "ALL")
  
  ## shortcut, y-X\beta
  yXbe <- parts$y - tcrossprod(parts$X, t(parts$beta))
  ## total variance
  uluti <- length(parts$theta)
  Zlam <- tcrossprod(parts$Z, parts$Lambdat)
  V <- (tcrossprod(Zlam, Zlam) + Diagonal(parts$n, 1)) * (parts$sigma)^2
  M <- solve(chol(V))
  invV <- tcrossprod(M, M)
  LambdaInd <- parts$Lambda
  LambdaInd@x[] <- seq(1:uluti) 
  invVX <- crossprod(parts$X, invV)
  Pmid <- solve(crossprod(parts$X, t(invVX)))
  P <- invV - tcrossprod(crossprod(invVX, Pmid), t(invVX))
  
  ## block 1 fixhes
  fixvar <- solve(tcrossprod(crossprod(parts$X, invV), t(parts$X)))
  if (full == FALSE) {
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
      devLambda[[i]] <- forceSymmetric(LambdaInd==i, uplo = "L")
      devV[[i]] <- tcrossprod(tcrossprod(parts$Z, t(devLambda[[i]])), parts$Z)
    }
    devV[[(uluti+1)]] <- Diagonal(nrow(parts$X), 1)
    
    ## Block 4: sigma's second derivative.
    ranhes <- matrix(NA, nrow = (uluti + 1), ncol = (uluti + 1))
    ## combinations (allow repeat) to indicate entries
    entries <- rbind(matrix(rep(1: (uluti + 1), each = 2),
                            (uluti + 1), 2, byrow = TRUE), t(combn((uluti + 1), 2)))
    entries <- entries[order(entries[,1], entries[,2]), ]
    ## ML estimates
    if (object@devcomp$dims[10] == 0) {
      if(information == "expected") {
        ranhes[lower.tri(ranhes, diag = TRUE)] <- 
          apply(entries, 1, function(x) as.numeric((1/2) *
                lav_matrix_trace(tcrossprod(tcrossprod(crossprod(invV,
                devV[[x[1]]]), invV), t(devV[[x[2]]])))))
      }
      if(information == "observed") {
        ranhes[lower.tri(ranhes, diag = TRUE)] <- apply(entries, 1, 
          function(x) -as.numeric((1/2) *
          lav_matrix_trace(tcrossprod(tcrossprod(crossprod(invV,
          devV[[x[1]]]), invV), t(devV[[x[2]]])))) +
          tcrossprod((tcrossprod((crossprod(yXbe,
          tcrossprod(tcrossprod(crossprod(invV,
          devV[[x[1]]]), invV), t(devV[[x[2]]])))), invV)), t(yXbe)))
      }      
    }
    ## REML estimates
    if (object@devcomp$dims[10] == 2|object@devcomp$dims[10] == 1) {
      if(information == "expected") {
        ranhes[lower.tri(ranhes, diag = TRUE)] <- apply(entries, 1, 
          function(x) as.numeric((1/2) * lav_matrix_trace(tcrossprod(
          tcrossprod(crossprod(P, devV[[x[1]]]), P), t(devV[[x[2]]])))))
      }
      if(information == "observed") {
        ranhes[lower.tri(ranhes, diag = TRUE)] <- apply(entries, 1, 
          function(x) -as.numeric((1/2) * lav_matrix_trace(
          tcrossprod(tcrossprod(crossprod(P, devV[[x[1]]]), P), 
          t(devV[[x[2]]])))) + tcrossprod((tcrossprod((crossprod(
          yXbe, tcrossprod(tcrossprod(crossprod(invV,
          devV[[x[1]]]), P), t(devV[[x[2]]])))), invV)), t(yXbe)))             
      }      
    }
    
    ranhes <- forceSymmetric(ranhes, uplo = "L")
    
    ## block 2 and block 3: second derivative of sigma and beta.
    if (information == "expected"){
      varcov_beta <- matrix(0, length(devV), length(parts$beta))
    }
    if (information == "observed"){
      varcov_beta <- matrix(NA, length(devV), length(parts$beta))
      for (j in 1 : (length(devV))) {
        varcov_beta[j,] <- as.vector(tcrossprod(crossprod(parts$X, 
        (tcrossprod(crossprod(invV, devV[[j]]), invV))), t(yXbe)))
      }
    }
    ## Organize full_varcov
    ## reprameterize sd and var
    if (ranpar == "var") {
        ranhes <- ranhes
    } else if (ranpar == "sd") {
        sdcormat <- as.data.frame(VarCorr(x,comp = "Std.Dev"), order = "lower.tri")
        sdcormat$sdcor2[which(is.na(sdcormat$var2))] <- sdcormat$sdcor[which(is.na(sdcormat$var2))]*2
        sdcormat$sdcor2[which(!is.na(sdcormat$var2))] <- sdcormat$vcov[which(!is.na(sdcormat$var2))]/
        sdcormat$sdcor[which(!is.na(sdcormat$var2))]
        score_varcov <- sweep(score_varcov, MARGIN = 2, sdcormat$sdcor, `*`)
      }
    }
    full_varcov <- solve(rbind(cbind(fixhes, t(varcov_beta)),
                               cbind(varcov_beta, ranhes)))
    
    colnames(full_varcov) <- c(names(parts$fixef), paste("cov",
                             names(parts$theta), sep = "_"), "residual")
    
    callingFun <- try(deparse(sys.call(-2)), silent = TRUE)
    if(length(callingFun) > 1) callingFun <- paste(callingFun, collapse="")
    if(!inherits(callingFun, "try-error") & grepl("summary.merMod", callingFun)){
      return(fixvar)
    } else {
      return(full_varcov)
    }
  }
}
