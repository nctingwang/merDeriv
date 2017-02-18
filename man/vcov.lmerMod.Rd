\name{vcov.lmerMod}
\alias{vcov.lmerMod}
\title{Extract Variance-Covariance Matrix of all Parameters for Linear 
  Mixed Effects Models}
\description{This function calculates the variance-covariance 
  matrix for all parameters (fixed, random effect, and residual) in a linear 
  mixed effects model of class \code{\linkS4class{lmerMod}}.  
}
\usage{
  \method{vcov}{lmerMod}(object, full = TRUE, \dots)
}
\arguments{
  \item{object}{An object of class \code{\linkS4class{lmerMod}}.}
  \item{full}{If \code{full = TRUE} (default), the variance-covariance matrix 
    for all fitted parameters (including fixed effect parameters, 
    random effect (co)variances, and residual variance.  
    If \code{full = FALSE}, the variance-covariance matrix of only fixed effect 
    parameters is returned.} 
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
    vcov(lme4fit)
  }
}
