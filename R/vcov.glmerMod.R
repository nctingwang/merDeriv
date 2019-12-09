vcov.glmerMod <- function(object, ranpar = "var", ...) {
  dotdotdot <- list(...)
  if("full" %in% names(dotdotdot)){
    full <- dotdotdot$full
  } else {
    full <- FALSE
  }
  
  if (full == FALSE) {
    full_vcov <- vcov.merMod(object)
  } else {
    ## Hessian was based on deviance function, which is the 
    ## -2*LogLik. That's why divided by -2
    full_vcov_noorder <- -solve(object@optinfo$derivs$Hessian/(-2))
  
    ## Block order in Hessian was theta, beta. Reorganize to 
    ## put fixed parameter block first to match with score 
    ## matrix order.
  
    ## count parameter numbers
    p <- nrow(full_vcov_noorder)
    pfix <- length(object@beta)
    pran <- length(object@theta)
    ## reorder four blocks
    full_vcov <- matrix(NA, nrow(full_vcov_noorder), ncol(full_vcov_noorder))
    full_vcov[1:pfix, 1:pfix] <- full_vcov_noorder[(pran + 1):p, (pran + 1):p]
    full_vcov[(pfix + 1):p, (pfix + 1): p] <- full_vcov_noorder[1:pran, 1:pran]
    full_vcov[(pfix + 1):p, 1:pfix] <- full_vcov_noorder[1:pran, (pran + 1): p]
    full_vcov[1:pfix, (pfix + 1): p] <- full_vcov_noorder[(pran + 1): p, 1:pran]
    

    ## reparameterize for sd and var for random variance/covariance parameters.
    if (pran ==1 & (ranpar == "sd"|ranpar == "theta")) {
       full_vcov <- full_vcov
    } else if (ranpar == "var") {
      library("numDeriv")
      useSc <- attr(x,"useSc")
      dd <- lme4:::devfun2(object,useSc=TRUE,signames=FALSE, transfuns = 
                             list(from.chol = Cv_to_Sv,
                              to.sd = identity))
      vdd <- as.data.frame(VarCorr(object,comp = "Std.Dev"), order = "lower.tri")
      ranpars <- vdd[,"sdcor"]
      pars <- c(ranpars,object@beta)
      hh1 <- hessian(dd,pars)
      vv2 <- -2*solve(hh1)
        } else {
        stop("ranpar needs to be theta, sd or var")
    }
    } 

    ## name the matrix
    parts <- getME(object, c("fixef", "theta"))
    colnames(full_vcov) <- c(names(parts$fixef), paste("cov",
                                                         names(parts$theta), sep="_"))
  }
  return(full_vcov)
}

