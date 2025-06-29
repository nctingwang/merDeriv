\name{vcov.lmerMod}
\alias{vcov.lmerMod}
\title{Extract Variance-Covariance Matrix of all Parameters for Linear 
  Mixed Effects Models}
\description{This function calculates the variance-covariance 
  matrix for all parameters (fixed, random effect, and residual) in a linear 
  mixed effects model of class \code{\link[lme4]{lmerMod-class}}.  
}
\usage{
\method{vcov}{lmerMod}(object, \dots)
}
\arguments{
  \item{object}{An object of class \code{\link[lme4]{lmerMod-class}}.}
  \item{\dots}{additional arguments, including \code{full}, 
     \code{information} and \code{ranpar} (\code{full = FALSE},  
     \code{information = "expected"} and \code{ranpar = "var"} 
     are default; see details).}
}
\value{
  A \emph{p} by \emph{p} variance-covariance matrix, where \emph{p}
  represents the number of parameters. 
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
}
\examples{
\dontrun{
# The sleepstudy example
lme4fit <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE)

# variance covariance matrix for all parameters
vcov(lme4fit, full = TRUE, ranpar = "var")
}
}
