\name{vcov.lmerMod}
\alias{vcov.lmerMod}
\title{Extract Variance-Covariance Matrix of all Parameters for Linear 
  Mixed Effects Models}
\description{This function calculates the variance-covariance 
  matrix for all parameters (fixed, random effect, and residual) in a linear 
  mixed effects model of class \code{\linkS4class{lmerMod}}.  
}
\usage{
\method{vcov}{lmerMod}(object, full = FALSE, information = "expected", \dots)
}
\arguments{
  \item{object}{An object of class \code{\linkS4class{lmerMod}}.}
  \item{full}{ If \code{full = FALSE} (default), returns the 
    variance-covariance matrix of only fixed effect 
    parameters. If \code{full = TRUE} , returns the variance-covariance matrix 
    for all fitted parameters (including fixed effect parameters, 
    random effect (co)variances, and residual variance).  
   }
  \item{information}{If \code{information = "expected"} (default),
    the variance-covariance matrix 
    is based on the inversion of Fisher information matrix.   
    If \code{information = "observed"}, the variance-covariance matrix
    is based on the observed Fisher information, which is the negative
    of Hessian matrix.}
  \item{\dots}{additional arguments.}
}
\value{
  A \emph{p} by \emph{p} variance-covariance matrix, where \emph{p}
  represents the number of parameters.
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
