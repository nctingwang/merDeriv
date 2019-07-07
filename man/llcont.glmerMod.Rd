\name{llcont.glmerMod}
\alias{llcont.glmerMod}
\title{Extract Cluster-wise Log Likelihoods for Generalized Linear Mixed Effects Models}
\description{A function for extracting the cluster-wise log likelihoods
  of a generalized linear mixed effects model fitted via \pkg{lme4}. 
  This function returns the cluster-wise log likelihoods, evaluated at the
  ML estimates.
}
\usage{
\method{llcont}{glmerMod}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{\linkS4class{glmerMod}}.}
  \item{\dots}{Additional arguments.}
}
\value{A vector of log-likelihoods whose length is
  the number of clusters. }
\references{
  Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). 
  Fitting Linear Mixed-Effects Models Using lme4. \emph{Journal of
  Statistical Software}, \bold{67}(1), 1-48. \doi{10.18637/jss.v067.i01}.
  
  Merkle, E. C., You, D. and Preacher, K. J., 2016. Testing Nonnested Structural 
  Equation Models. \emph{Psychological Methods}, \bold{21}(2), 151. 
  \url{https://arxiv.org/pdf/1402.6720v3}
}
\examples{
\dontrun{
# The cbpp example
lme4fit <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 data = cbpp, family = binomial)

# clusterwise log likelihood
llcont(lme4fit)
}  
}
