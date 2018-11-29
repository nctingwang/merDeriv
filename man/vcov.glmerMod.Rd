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
  \item{\dots}{additional arguments, including \code{full} 
  (\code{full = FALSE} is default; see details).}
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
  \code{lme4}. 
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

    # variance covariance matrix for all parameters
    vcov(lme4fit, full = TRUE)
}
}
