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
\method{estfun}{lmerMod}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{\linkS4class{lmerMod}}.}
  \item{\dots}{additional arguments, including \code{level}.  If \code{level = 1}, 
  case-wise scores are returned.  If \code{level = 2} (the default), 
  cluster-wise scores are returned.}
}
\value{An \emph{g} (or \emph{n}) by \emph{p} score matrix, corresponding to 
  \emph{g} clusters (\emph{n} observations) and \emph{p} parameters. 
  For models with multiple clustering variables (three-level models, 
  crossed random effects), an error is thrown if \code{level==2}.
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

