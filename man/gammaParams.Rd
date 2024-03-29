% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesCompanion.R
\name{gammaParams}
\alias{gammaParams}
\title{Compute params of gamma distribution}
\usage{
gammaParams(shape = NA, rate = NA, m = NA, sd = NA, mode = NA, verbose = FALSE)
}
\arguments{
\item{shape}{parameter}

\item{rate}{parameter}

\item{m}{mean of the distribution}

\item{sd}{standard deviation of the distribution}

\item{mode}{of the distribution}

\item{verbose}{= FALSE set to TRUE to get more explanation}
}
\value{
A vector of all possible values of parameters
}
\description{
Compute params of gamma distribution
from other meaningfull params.
Only the couples (shape, rate) or
(mean, sd) or (mean, mode) are implemented .
}
\examples{
gammaParams(m=0.5, mode=0.5) # warning
gammaParams(shape=0.6, rate=1) # OK 
gammaParams(m=0.5, sd=0.2)     # OK 
gammaParams(mode=0.5, sd=5)    # OK
}
\author{
Julien Bousquet (2022)
}
