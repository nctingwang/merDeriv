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
  \item{x}{An object of class \code{\link[lme4]{glmerMod-class}}.}
  \item{\dots}{Additional arguments, including \code{ranpar}
    (\code{ranpar = "var"} is default; see details)}.
}
\value{A \emph{g} by \emph{p} score matrix, corresponding to 
  \emph{g} clusters and \emph{p} parameters. 
  For models with multiple clustering variables (three-level models, 
  crossed random effects), an error is thrown.
  If \code{ranpar = "var"}, the random effects are parameterized as 
  variance/covariance; If \code{ranpar = "sd"}, the random effects are 
  parameterized as standard deviation/correlation; If \code{ranpar = "theta"}, 
  the random effects are parameterized as components of Cholesky decomposition. 
}
\references{
  Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). 
  Fitting Linear Mixed-Effects Models Using lme4. \emph{Journal of
  Statistical Software}, \bold{67}(1), 1-48. \doi{10.18637/jss.v067.i01}. 
}
\examples{
\dontrun{
data(finance, package = "smdata")

lme4fit <- glmer(corr ~ jmeth + (1 | item), data = finance,
                 family = binomial, nAGQ = 20)

# clusterwise scores
estfun(lme4fit, ranpar = "var")
}
}

