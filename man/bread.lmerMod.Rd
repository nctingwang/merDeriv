\name{bread.lmerMod}
\alias{bread.lmerMod}
\title{Extract Bread Component for Huber-White Sandwich Estimator of Linear Mixed Effects Models}
\description{
  This function calculates the bread component of the Huber-White 
  sandwich estimator (variance covariance matrix multiplied by 
  the number of clusters) for a linear mixed effects model 
  of class \code{\linkS4class{lmerMod}}.  
}
\usage{
\method{bread}{lmerMod}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{\linkS4class{lmerMod}}.}
  \item{\dots}{additional arguments, including \code{full} 
    \code{information}, and \code{ranpar} (\code{full = FALSE},
    \code{information = "expected"} and \code{ranpar = "var"}
    are default; see details).}
}
\value{
  A \emph{p} by \emph{p} "bread" matrix for the Huber-White sandwich estimator 
  (variance-covariance matrix multiplied by the number of clusters), where
  \emph{p} represents the number of parameters.  
  If \code{full = FALSE}, returns the 
  variance-covariance matrix of only fixed effect 
  parameters. If \code{full = TRUE} , returns the variance-covariance matrix 
  for all fitted parameters (including fixed effect parameters, 
  random effect (co)variances, and residual variance. 
  If \code{information = "expected"}, the variance-covariance matrix 
  is based on the inversion of Fisher information matrix.   
  If \code{information = "observed"}, the variance-covariance matrix
  is based on the observed Fisher information, which is the negative
  of Hessian matrix.  If \code{ranpar = "var"}, the random effects are 
  parameterized as variance/covariance; If \code{ranpar = "sd"}, 
  the random effects are parameterized as standard deviation/correlation. 
}
\references{
  Wang, T. & Merkle, E. C. (2018). Derivative Computations and Robust Standard
  Errors for Linear Mixed Effects Models in lme4. \emph{Journal of Statistical 
  Software}, \bold{87}(1), 1-16. \doi{10.18637/jss.v087.c01}

  Zeileis, A. (2006). Object-Oriented Computation of Sandwich Estimators. 
  Journal of Statistical Software, \bold{16}(9), 1-16. 
  \url{https://www.jstatsoft.org/v16/i09/}
}
\examples{
\dontrun{
# The sleepstudy example
lme4fit <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE)

# bread component for all parameters
bread(lme4fit, full = TRUE, information = "expected", ranpar = "var")
}
}
