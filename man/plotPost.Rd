% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotPost.R
\name{plotPost}
\alias{plotPost}
\title{Plot posterior MCMC
Plot the histogram of posterior MCMC, writen by
John K. Kruschke in Doing Bayesian Data Analysis, 2nd edition}
\usage{
plotPost(
  paramSampleVec,
  cenTend = "mode",
  compVal = NULL,
  ROPE = NULL,
  credMass = 0.95,
  HDItextPlace = 0.7,
  xlab = NULL,
  xlim = NULL,
  yaxt = NULL,
  ylab = NULL,
  main = NULL,
  cex = NULL,
  cex.lab = NULL,
  col = NULL,
  border = NULL,
  showCurve = FALSE,
  breaks = NULL,
  ...
)
}
\arguments{
\item{paramSampleVec}{the vector of the MCMC parameter's sample}

\item{cenTend}{"mode" by default for Central Tendency. Possible other values are "median", or "mean".}

\item{compVal}{numeric : value to compare.}

\item{ROPE}{vector of 2 limits defining the Range Of Practical Equivalence}

\item{credMass}{the mass of credibility in the High Density Interval}

\item{HDItextPlace}{by default 0.7 : position of the text "HDI".}

\item{showCurve}{a boolean. FALSE by default. If TRUE, plot the curve of density instead off histogram.}

\item{...}{other parameters given to hist() function, such xlab, xlim, yaxt, etc.}
}
\value{
Plot the histogram of the mcmcm values. Return the vector of mcmc summary values.
}
\description{
Plot posterior MCMC
Plot the histogram of posterior MCMC, writen by
John K. Kruschke in Doing Bayesian Data Analysis, 2nd edition
}
\examples{
plotPost(rnorm(2000), ROPE=c(-1,1), compVal=0)
}
\author{
John K. Kruschke (2015)
}
