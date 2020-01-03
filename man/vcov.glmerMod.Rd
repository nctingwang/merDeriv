\name{vcov.glmerMod}
\alias{vcov.glmerMod}
\title{Extract Variance-Covariance Matrix of all Parameters for Generalized 
  Linear Mixed Effects Models}
\description{This function calculates the variance-covariance 
  matrix for all parameters (fixed and random effect) in a generalized linear 
  mixed effects model of class \code{\linkS4class{glmerMod}}.  
}
\usage{
\method{vcov}{glmerMod}(object, \dots)
}
\arguments{
  \item{object}{An object of class \code{\linkS4class{glmerMod}}.}
  \item{\dots}{additional arguments, including \code{full} and \code{ranpar}
  (\code{full = FALSE} and \code{ranpar = "var"} are default; see details).}
}
\value{
  A \emph{p} by \emph{p} variance-covariance matrix, where \emph{p}
  represents the number of parameters. 
  If \code{full = FALSE}, returns the 
  variance-covariance matrix of only fixed effect 
  parameters. If \code{full = TRUE} , returns the variance-covariance matrix 
  for all fitted parameters (including fixed effect parameters and random 
  effect (co)variances). The variance-covariance matrix
  is based on the negative of Hessian matrix, which is extracted from 
  \code{lme4}.  If \code{ranpar = "var"}, the random effects are
  parameterized as variance/covariance; If \code{ranpar = "sd"},
  the random effects are parameterized as standard
  deviation/correlation; If \code{ranpar = "theta"}, 
  the random effects are parameterized as components of Cholesky decomposition. 
}
\references{
  Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). 
  Fitting Linear Mixed-Effects Models Using lme4. \emph{Journal of
  Statistical Software}, \bold{67}(1), 1-48. \doi{10.18637/jss.v067.i01}. 
}
\examples{
\dontrun{
# The cbpp example
data(finance, package="smdata")

lme4fit <- glmer(corr ~ jmeth + (1 | item), data = finance,
                 family = binomial, nAGQ = 20)

# variance covariance matrix for all parameters
vcov(lme4fit, full = TRUE, ranpar = "var")
}
}
