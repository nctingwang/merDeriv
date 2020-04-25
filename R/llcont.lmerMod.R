llcont.lmerMod <- function(x, ...) {
  if (!is(x, "lmerMod")) stop("llcont.lmerMod() only works for lmer() models.")
  if (x@devcomp$dims[10] != 0) stop("llcont.lmerMod() only works for ML estimation.")
  if (!is.null(x@call$weights)) stop ("Models with weights specification is currently not supported.")
  if (length(grep("cbind", x@call$formula))!=0) stop ("Models with cbind specification is currently not supported.")  

  dotdotdot <- list(...)
  if("level" %in% names(dotdotdot)){
    level <- dotdotdot$level
  } else {
    level <- 2
  }
  if(!(level %in% c(1L, 2L))) stop("invalid 'level' argument supplied")
  
  ## error for variance equal to 0 or non-positive definite
  cor <- attr(lme4::VarCorr(x)[[1]], "correlation")
  ndim = nrow(cor)
  if (ndim == 1 & unclass(VarCorr(x))[[1]]  == 0) stop ("Random effect's variance is close to 0. Please check the model estimation")
  if (ndim > 1 & !matrixcalc::is.positive.definite(as.matrix(vcov.merMod(x), method = c("chol")))) stop ("Random effect's variance covariance matrix is not positive definite. Please check the model estimation")
  
  ## get all elements by getME and exclude multiple random effect models.
  parts <- getME(x, "ALL")
  
  ## prepare shortcuts
  yXbe <- parts$y - tcrossprod(parts$X, t(parts$beta))
  Zlam <- tcrossprod(parts$Z, parts$Lambdat)
  V <- (tcrossprod(Zlam, Zlam) + diag(1, parts$n)) * (parts$sigma)^2
  M <- solve(chol(V))
  invV <- tcrossprod(M, M)
  yXbesoV <- crossprod(yXbe, invV)
  
  ## get each observation's log likelihood
  ll <- -(1/2) * log(2*pi) - (1/2) * log(diag(lu(V)@U)) - 
    (1/2) * t(yXbesoV) * yXbe
  ll <- as.numeric(ll)

  ## Clusterwise if level==2
  if (level == 2) {
    index <- parts$flist[[1]]
    lltmp <- aggregate(x = ll, by = list(index), FUN = sum)
    ll <- lltmp[,-1]  
    names(ll) <- lltmp[,1]
  }
  
  return(ll)
}
