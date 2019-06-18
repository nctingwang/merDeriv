estfun.glmerMod <- function(x){
    if (!is(x, "glmerMod")) stop("estfun.glmerMod() only works for glmer() models.")
    if (length(x@theta) > 1 warning ("scores may be not accurate due to the
fact that nAGQ = 1 is implemented in lme4 model estimation with multiple random effects")
  ## log-likelihood contributions of a glmer() model with
  ## one grouping variable (no crossed or nested random effects
  ## allowed for now). Much code is taken from Yves.
    
  ## extract nAGQ used in model fit. 
  ngq <- x@devcomp$dims[7]  
  
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
  if(ndim == 1){
    XW <- lavaan:::lav_integration_gauss_hermite(n    = ngq,
                                                 ndim = ndim)
  } else {
    XW <- lavaan:::lav_integration_gauss_hermite(n    = ngq,
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
                  formula = fit@call$formula, frame = fit@frame)) %*% w.star)

      score[j,] <- out[j,-1]/out[j,1]
    } else {
      ## from integration3_cfa.R (multivariate version)
      ## FIXME: if >1 grouping var, length(re.vars) > 1
      C <- t(chol(re.vars[[1]][,,j]))
      etamns <- re.modes[[1]][j,]

      x.star <- t(as.matrix(C %*% t(XW$x) + as.numeric(etamns)))
      w.star <- XW$w * (2*pi)^(ndim/2) * det(C) * exp(0.5 * apply(XW$x, 1, crossprod)) * lavaan:::lav_mvnorm_dmvnorm(x.star, Mu = rep(0, ndim), Sigma = VarCov[[1]], log = FALSE)

      out[j,] <- (t(score.prod(S = x.star,
                               Xi = as.matrix(X[grps==grpnm[j],]),
                               Y = Data[grps==grpnm[j]],
                               fe.pred = fe.pred[grps==grpnm[j]],
                               Zi = as.matrix(Z[grps==grpnm[j],]),
                               re.modes = re.modes[[1]],
                               grp = as.numeric(grpnm[j]), fam = fam, 
                               devLambda = devLambda, Lambda = parts$Lambda,
                               iLambda = iLambda,
                               formula = fit@call$formula, frame = fit@frame)) %*% w.star)

      score[j,] <- out[j,-1]/out[j,1]
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
    eta <- if(is.matrix(S)) t(S[i,,drop = FALSE]) else x[i]
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


if(FALSE){
  library("lme4")
  source("llcont.glmerMod.R")
  source("estfun.glmerMod.R")

  
  ## fit Rasch model in mirt
  library("mirt")
  data <- expand.table(LSAT7)[1:200,]
  mod <- mirt(data, 1, itemtype="Rasch")

  ## fit model in lme4
  library("reshape2")
  data$person <- as.numeric(rownames(data))
  datalong <- melt(data, id = c("person"))
  fit <- glmer(value ~ -1 + variable + (1|person), family = binomial, 
               data = datalong, nAGQ = 20)
  ## score from lme4
  score <- estfun.glmerMod(fit)
  colSums(score)
  
  ## compare score obtained from integrate()
  ## defin jointdist() and fun2()
  jointdist <- function(theta, x, alph, bet){
    res <- matrix(NA, length(theta), length(alph))
    for(j in 1:length(theta)){
        res[j,] <- dbinom(x, size=1, prob=plogis(alph*theta[j] + bet))
    }
    apply(res, 1, prod) * dnorm(theta)
  }

   # deriv wrt slope
  fun1 <- function(theta, x, alph, bet, xvec, avec, bvec){
    (theta * (x - plogis(alph*theta + bet))) *
        jointdist(theta, xvec, avec, bvec)
  }
 
   ## deriv wrt intercept
   fun2 <- function(theta, x, alph, bet, xvec, avec, bvec){
     (x - plogis(alph*theta + bet)) * jointdist(theta, xvec, avec, bvec)
   }

  ## construct the score matrix by integrate()
  data <- expand.table(LSAT7)[1:200,]
  score2 <- matrix(NA, nrow(data), (2*ncol(data)))
  ## ML estimates of item parameters:
  itempars <- as.numeric(sapply(coef(mod), function(x) x[1:2]))
  # equal slops: sqrt(1.023)
  itempars[((1:5)*2 - 1)] <- itempars[((1:5)*2 - 1)]*sqrt(itempars[12]) 

  avec <- itempars[((1:5)*2 - 1)] # equal slopes
  bvec <- itempars[(1:5)*2] # intercepts
  for(i in 1:nrow(data)){
      normconst <- integrate(jointdist, -6, 6, x= as.numeric(data[i,]),
                             alph=avec, bet=bvec)$value

      for(j in 1:5){
        par1loc <- (j-1)*2 + 1

        score2[i,j] <- (1/normconst) * integrate(fun2, -6, 6,
                                             x = as.numeric(data[i,j]),
                                             alph = itempars[par1loc],
                                             bet = itempars[(par1loc+1)],
                                             xvec = as.numeric(data[i,]),
                                             avec = avec, bvec = bvec)$value
        score2[i,(5+j)] <- (1/normconst) * integrate(fun1, -6, 6, x=as.numeric(data[i,j]),
                                                   alph=itempars[par1loc],
                                                   bet=itempars[(par1loc+1)],
                                                   xvec = as.numeric(data[i,]), avec=avec,
                                                   bvec=bvec)$value
     }
  }

  ## compare these two methods
  ## fixed mathed
  plot(score[,c(1:5)], score2[,c(1:5)])
  ## random matched
  plot(score[,6], apply(score2[,c(6:10)], 1, sum))

  
  ## compare with estfun.AllModelClass in mirt. 
  library("mirt")
  data <- expand.table(LSAT7)
  mod <- mirt(data, 1, itemtype="Rasch")

  ## fit model in lme4
  library("reshape2")
  ## extract from mirt
  data$person <- as.numeric(rownames(data))
  datalong <- melt(data, id = c("person"))
  fit <- glmer(value ~ -1 + variable + (1|person), family = binomial, 
               data = datalong, nAGQ = 20)
  ## score from lme4
  score <- estfun.glmerMod(fit)
  mod1 <- mirt(expand.table(LSAT7), 1, SE = TRUE, SE.type = "crossprod")
  sc1 <- estfun.AllModelClass(mod1)
  ## reorganize sc1 to corresponds to the parameters
  fix <- sc1[,c(2, 4, 6, 8, 10)]
  random <- apply(sc1[,c(1,3, 5, 7, 9)], 1, sum)
  plot(score[,1:5], fix)
  plot(score[,6], random)
  plot(score, cbind(fix,random))

  ###########
  ## non-canonical link check
  fit2 <- glmer(value ~ -1 + variable + (1|person), family =
                binomial(link = "cloglog"), 
                data = datalong, nAGQ = 20)
  fixscore2 <- estfun.glmerMod(fit2)
  colSums(fixscore2)
  

#####################################################
  ## fit poisson

  # other families (poisson, etc). Make sure
  # aic() functions work similarly. 
  simfun <- function(ng = 20, nr = 100, fsd = 1, indsd = 0.2, b = c(1,
     2)) {
     ntot <- nr * ng
     b.reff <- rnorm(ng, sd = fsd)
     b.rind <- rnorm(ntot, sd = indsd)
     x <- runif(ntot)
     dd <- data.frame(x, f = factor(rep(c(1:ng), each = nr)),
         obs = 1:ntot)
     dd$eta0 <- model.matrix(~x, data = dd) %*% b
     dd$eta <- with(dd, eta0 + b.reff[f] + b.rind[obs])
     dd$mu <- exp(dd$eta)
     dd$y <- with(dd, rpois(ntot, lambda = mu))
     dd
  }
  dd <- simfun()
  m0 <- glmer(y ~ x + (1|f), family = "poisson", data = dd, nAGQ = 20)
  score <- estfun.glmerMod(m0)
  colSums(score)

 

 ## Models with more than 1 random effects can not have nAGQ > 1. Thus scores
 ## and log likelihood may not be precise. 
##############################################
  ## 3-dimensional, correlated random effects
  data(finance, package="smdata")
  fit2 <- glmer(corr ~ jmeth + (jmeth | item), data=finance,
                family=binomial)
  ll2 <- llcont.glmerMod(fit2)
  c(sum(ll2), logLik(fit2))
  

  estfun.glmerMod(fit3)
  colSums(sc3)


}
