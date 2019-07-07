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
    including \code{full} (\code{full = FALSE}; see details).}
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
lme4fit <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 data = cbpp, family = binomial)

# bread component for all parameters
bread(lme4fit, full = TRUE)
}
}
