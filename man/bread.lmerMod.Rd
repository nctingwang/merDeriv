\name{bread.lmerMod}
\alias{bread.lmerMod}
\title{Extract Bread Component for Huber-White Sandwich Estimator of Linear Mixed Effects Models}
\description{
  This function calculates the bread component of Huber-White 
  sandwich estimator (variance covariance matrix multiply 
  the number of clusters) for a linear mixed effects model, 
  of class \code{\linkS4class{lmerMod}}.  
}
\usage{
  \method{bread}{lmerMod}(object, full = TRUE, \dots)
}
\arguments{
  \item{object}{An object of class \code{\linkS4class{lmerMod}}.}
  \item{full}{If \code{full = TRUE} (the default), the bread component of 
    Huber-White sandwich estimator for all 
    fitted parameters, including fixed effect parameters, random 
    effect (co)variances, and residual variance are returned.  
    If \code{full = FALSE}, the bread component for only
    fixed effect parameters are returned.} 
  \item{\dots}{additional arguments.}
}
\value{
  A \emph{p} by \emph{p} "bread" matrix for the Huber-White sandwich estimator 
  (variance-covariance matrix multiplied by the number of clusters), where
  \emph{p} represents the number of parameters. 
}
\references{
  Wang, T. & Merkle, E. C. (2016). Derivative Computations and Robust Standard
  Errors for Linear Mixed Effects Models in lme4. 
  \url{https://arxiv.org/abs/1612.04911}

  Zeileis, A. (2006). Object-Oriented Computation of Sandwich Estimators. 
  Journal of Statistical Software, \bold{16}(9), 1-16. 
  \url{http://www.jstatsoft.org/v16/i09/}
}
\examples{
  \dontrun{
    # The sleepstudy example
    lme4fit <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE)

    # bread component for all parameters
    bread(lme4fit)
  }
}
