vcov.glmerMod <- function(object, ...) {
  dotdotdot <- list(...)
  if("full" %in% names(dotdotdot)){
    full <- dotdotdot$full
  } else {
    full <- FALSE
  }

  if("ranpar" %in% names(dotdotdot)){
    ranpar <- dotdotdot$ranpar
  } else {
    ranpar <- "var"
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
    if (ranpar == "theta") {
       full_vcov <- full_vcov
      } else if (ranpar == "sd"| ranpar == "var") {
        dd <- devfun2(object,useSc=FALSE,signames=TRUE)
        nvp <- length(attr(dd,"thopt"))
        pars <- attr(dd,"optimum")
        pars <- pars[!is.na(names(pars))] 
        hh <- hessian(dd, pars)/(-2)
        
      if (ranpar == "var"){
        sdcormat <- as.data.frame(VarCorr(object,comp = "Std.Dev"), order = "lower.tri")
        sdcormat$sdcor2[which(is.na(sdcormat$var2))] <- (1/2)*(sdcormat$sdcor[which(is.na(sdcormat$var2))])^(-1/2)
        sdcormat$sdcor2[which(!is.na(sdcormat$var2))] <- (-1)*(sdcormat$vcov[which(!is.na(sdcormat$var2))]/
          sdcormat$sdcor[which(!is.na(sdcormat$var2))])^(-1)
        hh[((pran + 1):p), (1:pran)] <- sweep(as.matrix(hh[((pran + 1):p), (1:pran)]), 
                                                 MARGIN = 1, sdcormat$sdcor2, `*`)
        hh[(1:pran), ((pran + 1):p)] <- sweep(as.matrix(hh[(1:pran), ((pran + 1):p)]), 
                                                 MARGIN = 1, sdcormat$sdcor2, `*`)
        ## ranhes reparameterization
        if (pran == 1){
            entries = matrix(1, 1, 1)
            weight <- (sdcormat$sdcor2)^2
        } else {
          entries <- rbind(matrix(rep(1: pran, each = 2),
                                pran, 2, byrow = TRUE), 
                         t(combn(pran, 2)))
          entries <- entries[order(entries[,1], entries[,2]), ]
          weight <- apply(entries, 1, function(x) sdcormat$sdcor2[x[1]] * sdcormat$sdcor2[x[2]])
        }
        hh[1:pran, 1:pran][lower.tri(hh[1:pran, 1:pran], 
                                     diag = TRUE)] <- weight * hh[1:pran, 1:pran][lower.tri(hh[1:pran, 1:pran], diag = TRUE)]

      }
        full_vcov_noorder <- -solve(hh)
        full_vcov <- matrix(NA, nrow(full_vcov_noorder), ncol(full_vcov_noorder))
        full_vcov[1:pfix, 1:pfix] <- full_vcov_noorder[(pran + 1):p, (pran + 1):p]
        full_vcov[(pfix + 1):p, (pfix + 1): p] <- full_vcov_noorder[1:pran, 1:pran]
        full_vcov[(pfix + 1):p, 1:pfix] <- full_vcov_noorder[1:pran, (pran + 1): p]
        full_vcov[1:pfix, (pfix + 1): p] <- full_vcov_noorder[(pran + 1): p, 1:pran]
    } else {
           warning("ranpar needs to be theta, sd or var")
    }
      

    ## name the matrix
    parts <- getME(object, c("fixef", "theta"))
    colnames(full_vcov) <- c(names(parts$fixef), paste("cov",
                                                         names(parts$theta), sep="_"))
  }
  return(full_vcov)
}

