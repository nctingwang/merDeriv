estfun.glmerMod <- function(x, ranpar = "var",...){
  ## log-likelihood contributions of a glmer() model with
  ## one grouping variable (no crossed or nested random effects
  ## allowed for now). Much code comes from Yves.
  if (!is(x, "glmerMod")) stop("estfun.glmerMod() only works for glmer() models.")
  ## check multiple groups.
  if (length(getME(x, "l_i")) > 1L) stop("Multiple cluster variables detected. This type of model is currently not supported.")
  if (length(x@theta) > 1) warning ("scores may be not accurate due to the fact that nAGQ = 1 is implemented in lme4 model estimation with multiple random effects")
  if (!is.null(x@call$weights)) stop ("Models with weights specification is currently not supported.")
  if (length(grep("cbind", x@call$formula))!=0) stop ("Models with cbind specification is currently not supported.")  
    
  ## extract nAGQ used in model fit, unless overridden by ...
  ddd <- list(...)
  if ("nAGQ" %in% names(ddd)){
    ngq <- ddd$nAGQ
  } else {
    ngq <- x@devcomp$dims[7]
  }
  
  ## 1a. obtain random effect predictions + sds from predict()
  ##    these become etamns and etasds below, removing
  ##    need for "adaptive" quadrature.
  ## (etamns are random effect means, etasds are random
  ## effect sds)
  fe.pred <- predict(x, re.form = NA)
  re.modes <- ranef(x, condVar = TRUE)
  re.vars <- vector("list", length(re.modes))
  for(i in 1:length(re.modes)){
    re.vars[[i]] <- attr(re.modes[[i]], "postVar")
  }
  re.b <- getME(x, "b")

  ## TODO: fix this if multiple grouping variables
  grps <- getME(x, "flist")[[1]]
  grpnm <- levels(grps)
  VarCov <- unclass(VarCorr(x))

  ## 1b. obtain family info
  fam <- x@call$family
  if(class(fam) == "name" | class(fam) == "character"){
    ## no link specified
    fam <- do.call(as.character(fam), list())
  } else {
    ## link specified
    fam <- eval(fam)
  }
  
  ## 2. obtain prediction for each observation,
  ##    including only fixed effects (random effects
  ##    added during quadrature). also observed y
  preds <- predict(x)
  Data <- getME(x, "y")
  ## Z matrix and X
  Z <- getME(x, "Z")
  X <- getME(x, "X")
  ## Prepare for random score. 
  parts <- getME(x, c("theta", "Lambda", "Lind"))
  uluti <- length(parts$theta)
  iLambda <- solve(parts$Lambda)
  devLambda <- vector("list", uluti)
  LambdaInd <- parts$Lambda
  LambdaInd@x[] <- seq(1:uluti)
  
  for (k in 1:uluti) {
    devLambda[[k]] <- LambdaInd==k
  }

  ## 3. Quadrature
  N <- nobs(x)
  ndim <- sapply(VarCov, nrow)
  ## FIXME this has length > 1 for crossed
  J <- getME(x, "l_i")
  #lik <- numeric(J)
  score <- matrix(NA, J, (ncol(X) + uluti))
  ## to contain pdf
  out <- matrix(NA, J, (ncol(X) + uluti + 1))  
        
  ## quadrature points:
  lav_integration_gauss_hermite <- getFromNamespace("lav_integration_gauss_hermite", "lavaan")
  if(ndim == 1){
    XW <- lav_integration_gauss_hermite(n    = ngq,
                                                 ndim = ndim)
  } else {
    XW <- lav_integration_gauss_hermite(n    = ngq,
                                                 ndim = ndim,
                                                 dnorm = TRUE)
  }

  for(j in 1:J){
    if(ndim == 1){
      etasds <- sqrt(as.numeric(re.vars[[1]][,,j]))
      etamns <- re.modes[[1]][j,]
    
      w.star <- sqrt(2) * etasds * dnorm(etasds * (sqrt(2)*XW$x) + etamns, 
                rep(0, ndim), sqrt(VarCov[[1]])) * exp(XW$x^2) * XW$w

      x.star <- etasds * (sqrt(2)*XW$x) + etamns

      out[j,] <- (t(score.prod(S = x.star,
                  Xi = as.matrix(X[grps==grpnm[j],]),
                  Y = Data[grps==grpnm[j]],
                  fe.pred = fe.pred[grps==grpnm[j]],
                  Zi = as.matrix(Z[grps==grpnm[j],]),
                  re.modes = re.modes[[1]],
                  grp = as.numeric(grpnm[j]), fam = fam,
                  devLambda = devLambda, Lambda = parts$Lambda,
                  iLambda = iLambda,
                  formula = x@call$formula, frame = x@frame)) %*% w.star)
      if ((ranpar == "sd")| (ranpar == "theta")){
          score[j,] <- out[j,-1]/out[j,1]
      } else {
          if (ranpar == "var"){
            ## get variance and sd 
            sdcormat <- as.data.frame(VarCorr(x,comp = "Std.Dev"), order = "lower.tri")
            ## chain rule
            sdcormat$sdcor2[which(is.na(sdcormat$var2))] <- (1/2) * 
              sdcormat$sdcor[which(is.na(sdcormat$var2))]^(-1/2)
            out[j,(ncol(X) + 1)] <- out[j,(ncol(X) + 1)] * sdcormat$sdcor2
            
            score[j,] <- out[j,-1]/out[j,1]
          } else{
          stop("ranpar needs to be sd, theta or var for glmerMod object.")
        }
      }

    } else {
      ## from integration3_cfa.R (multivariate version)
      ## FIXME: if >1 grouping var, length(re.vars) > 1
      C <- t(chol(re.vars[[1]][,,j]))
      etamns <- re.modes[[1]][j,]

      x.star <- t(as.matrix(C %*% t(XW$x) + as.numeric(etamns)))
      lav_mvnorm_dmvnorm <- getFromNamespace("lav_mvnorm_dmvnorm", "lavaan")
      w.star <- XW$w * (2*pi)^(ndim/2) * det(C) * exp(0.5 * apply(XW$x, 1, crossprod)) *     
        lav_mvnorm_dmvnorm(x.star, Mu = rep(0, ndim), Sigma = VarCov[[1]], log = FALSE)

      out[j,] <- (t(score.prod(S = x.star,
                               Xi = as.matrix(X[grps==grpnm[j],]),
                               Y = Data[grps==grpnm[j]],
                               fe.pred = fe.pred[grps==grpnm[j]],
                               Zi = as.matrix(Z[grps==grpnm[j],]),
                               re.modes = re.modes[[1]],
                               grp = as.numeric(grpnm[j]), fam = fam, 
                               devLambda = devLambda, Lambda = parts$Lambda,
                               iLambda = iLambda,
                               formula = x@call$formula, frame = x@frame)) %*% w.star)

          score[j,] <- out[j,-1]/out[j,1]
    }
  }
  if (ranpar=="theta"){
    score <- score
  } else {
    if (ranpar == "var") {
    ## create weight matrix
    require(matrixcalc)
    d0 <- (diag(1,nrow=ndim^2) + commutation.matrix(r=ndim)) %*% 
      (parts$Lambda[(1:ndim), (1:ndim)] %x% diag(1,nrow=ndim))
    L <- elimination.matrix(ndim)
    d1 <- L %*% d0 %*% t(L)
    dfin <- solve(d1)
    
    score[, ((ncol(X)+1):ncol(score))] <- as.matrix(score[, ((ncol(X)+1):ncol(score))] %*% dfin)
    } else {
      if (ranpar == "sd"){
        require(matrixcalc)
        d0 <- (diag(1,nrow=ndim^2) + commutation.matrix(r=ndim)) %*% 
          (parts$Lambda[(1:ndim), (1:ndim)] %x% diag(1,nrow=ndim))
        L <- elimination.matrix(ndim)
        d1 <- L %*% d0 %*% t(L)
        dfin <- solve(d1)
        
        score[, ((ncol(X)+1):ncol(score))] <- as.matrix(score[, ((ncol(X)+1):ncol(score))] %*% dfin)
        ## parameterize to sd and corr
        sdcormat <- as.data.frame(VarCorr(x,comp = "Std.Dev"), order = "lower.tri")
        sdcormat$sdcor2[which(is.na(sdcormat$var2))] <- sdcormat$sdcor[which(is.na(sdcormat$var2))]*2
        sdcormat$sdcor2[which(!is.na(sdcormat$var2))] <- sdcormat$vcov[which(!is.na(sdcormat$var2))]/
          sdcormat$sdcor[which(!is.na(sdcormat$var2))]
        score[, ((ncol(X)+1):ncol(score))] <- sweep(score[, ((ncol(X)+1):ncol(score))], MARGIN = 2, sdcormat$sdcor2, `*`)
      }
    }
  }
  score
}



score.prod <- function(S, Xi, Y = NULL, fe.pred, Zi, re.modes, grp, fam, devLambda, Lambda, iLambda, formula, frame, j) {
  Y <- as.numeric(Y)
  
  # number of quadrature points
  nQ <- if(is.matrix(S)) NROW(S) else length(S)
  
  # FIXME binomial n; pull this out of fam or fit?
  n <- 1
  ## contain glm score and likelihood
  out <- matrix(NA, nQ, (ncol(Xi) + length(devLambda) + 1))

  ## for non-canonical link, aphi is obtained from glm$dispersion
  if(fam$link != do.call(fam$family, list())$link){
    bracketrm <- gsub("\\(.*","", formula[[3]])[2]
    formglm <- paste0(paste0(as.character(formula)[2],
      as.character(formula)[1]), bracketrm)
    aphi <- summary(glm(formula(formglm), frame, family= fam[[1]]))$dispersion
  }
   tmpre <- as.matrix(re.modes)
 

  for(i in 1:nQ) {
    # add eta to fixed effects
    # TODO handle multiple grouping variables
    eta <- if(is.matrix(S)) t(S[i,,drop = FALSE]) else S[i]
    tmpre[grp,] <- eta
    linkyhat <- as.numeric(fe.pred + Zi %*% as.numeric(t(tmpre)))
    
    # use inverse link function on yhat
    yhat <- fam$linkinv(linkyhat)
    Zi_resid <- crossprod(Zi, as.matrix(Y - yhat))
    u  <- iLambda %*% as.vector(t(tmpre))
    iLam_dL <- lapply(devLambda, function(x) x %*% u)
    
    # score matrix with canonical link.
    ranscore <- rep(NA, length(devLambda))
    if(fam$link == do.call(fam$family, list())$link){
      glmscore <- crossprod(Xi, as.matrix(Y - yhat))
      ## random 
      for(k in 1:length(devLambda)){
        ranscore[k] <- crossprod(iLam_dL[[k]], Zi_resid)
      }
    } else {
      ## non-canonical link. 
      ## invD and invV. 
      invD <- diag(fam$mu.eta(linkyhat))
      invV <- aphi * diag(1/fam$variance(yhat))
      glmscore <- tcrossprod(tcrossprod(crossprod(Xi, invD),
        invV), t(as.matrix(Y-yhat)))
      Zi_resid <- tcrossprod(tcrossprod(crossprod(Zi, invD),
        invV), t(as.matrix(Y-yhat)))

      for(k in 1: length(devLambda)){
        ranscore[k] <- crossprod(iLam_dL[[k]], Zi_resid)
      }
    }
    #ranscoresum <- colSums(ranscore)
    ## add one to contain pdf.   
    likscore <- c(1, glmscore, ranscore)  

    # FIXME does fam$aic() work this way for other families
    # besides binomial?
    #logly <- -fam$aic(Y, n, yhat, wt = 1)/2      
    logly <- t(likscore) * exp(-fam$aic(Y, n, yhat, wt = 1)/2)
      
    #if(fam[1] =="poisson"){
    #  logly <- t(likscore) * exp(sum(dpois(Y, yhat, log = TRUE)))
    #}
    
    out[i,] <- colSums(logly)
  }
  
  out
}


