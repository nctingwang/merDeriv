vcov.glmerMod <- function(object, ...) {
  if (object@call$family!="binomial" & object@call$family!="poisson") stop("family has to be binomial or poisson) 

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
    if (length(getME(object, "l_i")) > 1L) stop("Multiple cluster variables detected. This type of model is currently not supported for full vcov.")
    
    ## Hessian was based on deviance function, which is the 
    ## -2*LogLik. That's why divided by -2
    if (object@devcomp$dims[['nAGQ']] == 0L) stop("For full vcov, nAGQ of at least 1 is required.")
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
    }
      
    if (ranpar == "sd") {
        dd <- devfun2(object,useSc=FALSE,signames=TRUE)
        nvp <- length(attr(dd,"thopt"))
        pars <- attr(dd,"optimum")
        pars <- pars[!is.na(names(pars))] 
        hh <- hessian(dd, pars)/(-2)

        full_vcov_noorder <- -solve(hh)
        full_vcov <- matrix(NA, nrow(full_vcov_noorder),
          ncol(full_vcov_noorder))
        full_vcov[1:pfix, 1:pfix] <-
          full_vcov_noorder[(pran + 1):p, (pran + 1):p]
        full_vcov[(pfix + 1):p, (pfix + 1): p] <-
          full_vcov_noorder[1:pran, 1:pran]
        full_vcov[(pfix + 1):p, 1:pfix] <-
          full_vcov_noorder[1:pran, (pran + 1): p]
        full_vcov[1:pfix, (pfix + 1): p] <-
            full_vcov_noorder[(pran + 1): p, 1:pran]
     }
        
    if (ranpar == "var"){
        dd <- devfun2(object,useSc=FALSE,signames=TRUE)
        nvp <- length(attr(dd,"thopt"))
        pars <- attr(dd,"optimum")
        pars <- pars[!is.na(names(pars))] 
        hh <- hessian(dd, pars)/(-2)
           
        sdcormat <- as.data.frame(VarCorr(object,comp = "Std.Dev"),
          order = "lower.tri")
        sdcormat$sdcor2[which(is.na(sdcormat$var2))] <-
          (1/2)*(sdcormat$sdcor[which(is.na(sdcormat$var2))])^(-1/2)
        sdcormat$sdcor2[which(!is.na(sdcormat$var2))] <- (-1)*
          (sdcormat$vcov[which(!is.na(sdcormat$var2))]/
             sdcormat$sdcor[which(!is.na(sdcormat$var2))])^(-1)
        hh[((pran + 1):p), (1:pran)] <- sweep(as.matrix(hh[((pran + 1):p),
          (1:pran)]), MARGIN = 2, sdcormat$sdcor2, `*`)
        hh[(1:pran), ((pran + 1):p)] <- t(hh[((pran + 1):p), (1:pran)])
        ## ranhes reparameterization
        if (pran == 1){
            entries = matrix(1, 1, 1)
            weight <- (sdcormat$sdcor2)^2
        } else {
          entries <- rbind(matrix(rep(1: pran, each = 2),
            pran, 2, byrow = TRUE), t(combn(pran, 2)))
          entries <- entries[order(entries[,1], entries[,2]), ]
          weight <- apply(entries, 1, function(x)
            sdcormat$sdcor2[x[1]] * sdcormat$sdcor2[x[2]])
        }
 
        hh[1:pran, 1:pran][lower.tri(hh[1:pran, 1:pran], diag = TRUE)] <-
          weight * hh[1:pran, 1:pran][lower.tri(hh[1:pran, 1:pran],
                                                diag = TRUE)]
        if (pran > 1){
          hh[1:pran, 1:pran] <- as.matrix(forceSymmetric(hh[1:pran, 1:pran],
            uplo = "L"))
        }
        
        full_vcov_noorder <- -solve(hh)
        full_vcov <- matrix(NA, nrow(full_vcov_noorder),
          ncol(full_vcov_noorder))
        full_vcov[1:pfix, 1:pfix] <-
          full_vcov_noorder[(pran + 1):p, (pran + 1):p]
        full_vcov[(pfix + 1):p, (pfix + 1): p] <-
          full_vcov_noorder[1:pran, 1:pran]
        full_vcov[(pfix + 1):p, 1:pfix] <-
          full_vcov_noorder[1:pran, (pran + 1): p]
        full_vcov[1:pfix, (pfix + 1): p] <-
          full_vcov_noorder[(pran + 1): p, 1:pran]
    }
    if (!(ranpar %in% c("sd", "theta", "var"))){
       stop("ranpar needs to be sd, theta or var for glmerMod object.")
    }
      
    ## name the matrix
    parts <- getME(object, c("fixef", "theta"))
    colnames(full_vcov) <- c(names(parts$fixef), paste("cov",
      names(parts$theta), sep="_"))
  }
  return(full_vcov)
}

