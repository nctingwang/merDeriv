\name{bread.glmerMod}
\alias{bread.glmerMod}
\title{Extract Bread Component for Huber-White Sandwich Estimator of
  Generalized Linear Mixed Effects Models}
\description{
  This function calculates the bread component of the Huber-White 
  sandwich estimator (variance covariance matrix multiplied by 
  the number of clusters) for a generalized linear mixed effects model 
  of class \code{\linkS4class{glmerMod}}.  
}
\usage{
\method{bread}{glmerMod}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{\linkS4class{glmerMod}}.}
  \item{\dots}{additional arguments,
    including \code{full} and \code{ranpar} (\code{full = FALSE},
    \code{ranpar = "var"}; see details).}
}
\value{
  A \emph{p} by \emph{p} "bread" matrix for the Huber-White sandwich estimator 
  (variance-covariance matrix based on observed Fisher information
  multiplied by the number of clusters), where
  \emph{p} represents the number of parameters.  
  If \code{full = FALSE}, returns the 
  variance-covariance matrix of only fixed effect 
  parameters. If \code{full = TRUE} , returns the variance-covariance matrix 
  for all fitted parameters (including fixed effect parameters, 
  random effect (co)variances, and residual variance.
  If \code{ranpar = "var"}, the random effects are parameterized as 
  variance/covariance; If \code{ranpar = "sd"}, the random effects are 
  parameterized as standard deviation/correlation; If \code{ranpar = "theta"}, 
  the random effects are parameterized as components of Cholesky
  decomposition.  
}
\references{
  Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). 
  Fitting Linear Mixed-Effects Models Using lme4. \emph{Journal of
  Statistical Software}, \bold{67}(1), 1-48. \doi{10.18637/jss.v067.i01}. 

  Zeileis, A. (2006). Object-Oriented Computation of Sandwich Estimators. 
  Journal of Statistical Software, \bold{16}(9), 1-16. 
  \url{http://www.jstatsoft.org/v16/i09/}
}
\examples{
\dontrun{
# The cbpp example
if (require("smdata")) {
   data(finance, package = "smdata")

   lme4fit <- glmer(corr ~ jmeth + (1 | item), data = finance,
                 family = binomial, nAGQ = 20)

  # bread component for all parameters
  bread(lme4fit, full = TRUE, ranpar = "var")
 }
}
}
