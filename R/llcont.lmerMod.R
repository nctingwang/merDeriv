llcont.lmerMod <- function(object, ...) {
  if (!is(object, "lmerMod")) stop("estfun.lmerMod() only works for lmer() models.")
  if (object@devcomp$dims[10] != 0) stop("estfun.lmerMod() only works for ML estimation.")
  
  ## get all elements by getME and exclude multiple random effect models.
  parts <- getME(object, "ALL")
  
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
  
  return(as.numeric(ll))
}
