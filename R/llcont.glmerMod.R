llcont.glmerMod <- function(x, ...){
  ## log-likelihood contributions of a glmer() model with
  ## one grouping variable (no crossed or nested random effects
  ## allowed for now). Much code is taken from Yves.

  if (!is(x, "glmerMod")) stop("llcont.glmerMod() only works for glmer() models.")
  ## check multiple groups.
  if (length(getME(x, "l_i")) > 1L) stop("Multiple cluster variables detected. This type of model is currently not supported.")
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

  ## 3. Quadrature
  N <- nobs(x)
  ndim <- sapply(VarCov, nrow)
  ## FIXME this has length > 1 for crossed
  J <- getME(x, "l_i")
  lik <- numeric(J)

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

  ## Z matrix
  Z <- getME(x, "Z")
  ## nonzero entries:
  ##Zrow <- Z@i + 1
  ##Zcol <- rep(seq_along(diff(Z@p)), diff(Z@p))

  for(j in 1:J){
    if(ndim == 1){
      etasds <- sqrt(as.numeric(re.vars[[1]][,,j]))
      etamns <- re.modes[[1]][j,]
    
      w.star <- sqrt(2) * etasds * dnorm(etasds * (sqrt(2)*XW$x) + etamns,
        rep(0, ndim), sqrt(VarCov[[1]])) * exp(XW$x^2) * XW$w

      x.star <- etasds * (sqrt(2)*XW$x) + etamns

      lik[j] <- sum(ly.prod(x.star, Data[grps==grpnm[j]],
                  fe.pred[grps==grpnm[j]],
                  as.matrix(Z[grps==grpnm[j],]),
                  re.modes[[1]],
                  as.numeric(grpnm[j]), fam) * w.star )
    } else {
      ## from integration3_cfa.R (multivariate version)
      ## FIXME: if >1 grouping var, length(re.vars) > 1
      C <- t(chol(re.vars[[1]][,,j]))
      etamns <- re.modes[[1]][j,]

      x.star <- t(as.matrix(C %*% t(XW$x) + as.numeric(etamns)))
      lav_mvnorm_dmvnorm <- getFromNamespace("lav_mvnorm_dmvnorm", "lavaan")
      w.star <- XW$w * (2*pi)^(ndim/2) * det(C) *
        exp(0.5 * apply(XW$x, 1, crossprod)) *
        lav_mvnorm_dmvnorm(x.star, Mu = rep(0, ndim),
        Sigma = VarCov[[1]], log = FALSE)

      lik[j] <- sum(ly.prod(S = x.star,
                  Y = Data[grps==grpnm[j]],
                  fe.pred = fe.pred[grps==grpnm[j]],
                  Zi = as.matrix(Z[grps==grpnm[j],]),
                  re.modes = re.modes[[1]],
                  grp = as.numeric(grpnm[j]),
                  fam = fam) * w.star )
    }
  }

  log(lik)
}

  
## Code here taken from Yves's function from integration1_cfa.R
## single observation, multiple values (quadrature points) for eta
## x are quadrature points
## Y is the data
ly.prod <- function(S, Y = NULL, fe.pred, Zi, re.modes, grp, fam) {
    Y <- as.numeric(Y)

    # number of quadrature points
    nQ <- if(is.matrix(S)) NROW(S) else length(S)

    # FIXME binomial n; pull this out of fam or fit?
    n <- 1

    out <- numeric(nQ)

    for(i in 1:nQ) {
        # add eta to fixed effects
        # TODO handle multiple grouping variables
        eta <- if(is.matrix(S)) t(S[i,,drop = FALSE]) else S[i]
        tmpre <- as.matrix(re.modes)
        tmpre[grp,] <- eta
        linkyhat <- as.numeric(fe.pred + Zi %*% as.numeric(t(tmpre)))
        
        # use inverse link function on yhat
        yhat <- fam$linkinv(linkyhat)
        
        # log-likelihood
        # FIXME does fam$aic() work this way for other families
        # besides binomial?
        logly <- -fam$aic(Y, n, yhat, wt = 1)/2
        #logly <- dbinom(Y, n, yhat, log = TRUE)

        out[i] <- exp(sum(logly))
    }

    out
}


