\name{llcont.lmerMod}
\alias{llcont.lmerMod}
\title{Extract Case-wise Log Likelihoods for Linear Mixed Effects Models}
\description{A function for extracting the case-wise log likelihoods
  of a linear mixed effects model fitted via \pkg{lme4}. 
  This function returns the case-wise log likelihoods, evaluated at the
  ML estimates.
}
\usage{
\method{llcont}{lmerMod}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{\link[lme4]{lmerMod-class}}.}
  \item{\dots}{additional arguments, including \code{level}
    (\code{level = 2} is default; see details).}
}
\value{If \code{level = 2}, a vector of log-likelihoods whose length is
  the number of clusters. If \code{level = 1}, a vector of length \emph{n}, 
  containing log-likelihoods for all \emph{n} observations.}
\references{
  Wang, T. & Merkle, E. C. (2018). Derivative Computations and Robust Standard
  Errors for Linear Mixed Effects Models in lme4. \emph{Journal of Statistical 
  Software}, \bold{87}(1), 1-16. \doi{10.18637/jss.v087.c01}
  
  Merkle, E. C., You, D. and Preacher, K. J., 2016. Testing Nonnested Structural 
  Equation Models. \emph{Psychological Methods}, \bold{21}(2), 151. 
  \url{https://arxiv.org/pdf/1402.6720v3}
}
\examples{
\dontrun{
# The sleepstudy example
lme4fit <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE)

# clusterwise log likelihood
llcont(lme4fit)
}  
}