\name{vcov.lmerMod}
\alias{vcov.lmerMod}
\title{Extract Variance-Covariance Matrix of all Parameters for Linear 
  Mixed Effects Models}
\description{This function calculates the variance-covariance 
  matrix for all parameters (fixed, random effect, and residual) in a linear 
  mixed effects model of class \code{\linkS4class{lmerMod}}.  
}
\usage{
\method{vcov}{lmerMod}(object, \dots)
}
\arguments{
  \item{object}{An object of class \code{\linkS4class{lmerMod}}.}
  \item{\dots}{additional arguments, including \code{full} 
  and \code{information} (\code{full = FALSE} 
  and \code{information = "expected"} 
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
  of Hessian matrix.
}
\references{
  Wang, T. & Merkle, E. C. (2016). Derivative Computations and Robust Standard
  Errors for Linear Mixed Effects Models in lme4. 
  \url{https://arxiv.org/abs/1612.04911}
}
\examples{
\dontrun{
    # The sleepstudy example
    lme4fit <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE)

    # variance covariance matrix for all parameters
    vcov(lme4fit, full = TRUE)
}
}
