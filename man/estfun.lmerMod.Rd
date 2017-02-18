\name{estfun.lmerMod}
\alias{estfun.lmerMod}
\title{Extract Case-wise and Cluster-wise Derivatives for Linear Mixed Effects Models}
\description{
  A function for extracting the case-wise and cluster-wise 
  derivatives of a linear mixed effects models fitted via \pkg{lme4}. 
  This function returns the case-wise and cluster-wise 
  scores, evaluated at the ML estimates.
}
\usage{
  \method{estfun}{lmerMod}(object, level = 2, \dots)
}
\arguments{
  \item{object}{An object of class \code{\linkS4class{lmerMod}}.}
  \item{level}{If \code{level = 1}, case-wise scores are returned. 
  If \code{level = 2} (the default), cluster-wise scores are returned.}
  \item{\dots}{additional arguments.}
}
\value{An \emph{n} by \emph{p} score matrix, corresponding to \emph{n} 
  observations and \emph{p} parameters. 
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

    # casewise scores
    estfun(lme4fit, level = 1)

    # clusterwise scores
    estfun(lme4fit, level = 2)
  }
}

