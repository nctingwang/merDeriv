\name{estfun.glmerMod}
\alias{estfun.glmerMod}
\title{Extract Cluster-wise Derivatives for Generalized Linear Mixed Effects Models}
\description{
  A function for extracting the cluster-wise 
  derivatives of a generalized linear mixed effects models fitted
  via \pkg{lme4}.  This function returns the cluster-wise 
  scores, evaluated at the ML estimates.
}
\usage{
\method{estfun}{glmerMod}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{\linkS4class{glmerMod}}.}
  \item{\dots}{Additional arguments.}
}
\value{A \emph{g} by \emph{p} score matrix, corresponding to 
  \emph{g} clusters and \emph{p} parameters. 
  For models with multiple clustering variables (three-level models, 
  crossed random effects), an error is thrown.
}
\references{
  Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). 
  Fitting Linear Mixed-Effects Models Using lme4. \emph{Journal of
  Statistical Software}, \bold{67}(1), 1-48. \doi{10.18637/jss.v067.i01}. 
}
\examples{
\dontrun{
# The cbpp example
lme4fit <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 data = cbpp, family = binomial)

# clusterwise scores
estfun(lme4fit)
}
}

