% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clpDataSim.R
\name{clpDataSim}
\alias{clpDataSim}
\title{Simulate ERR data from a compound lognormal Poisson distribution}
\usage{
clpDataSim(
  R,
  N,
  reduction,
  community_mean,
  cv_between,
  cv_within,
  cv_slide,
  cv_reduction,
  overall_mean = exp(mean(log(community_mean))),
  edt = 24
)
}
\arguments{
\item{R}{The number of replicate datasets}

\item{N}{The maximum number of individuals in total}

\item{reduction}{The geometric mean reduction (may be vectorised)}

\item{community_mean}{A vector of geometric mean pre-treatment EPG in each community}

\item{cv_between}{A vector of CV reflecting variation in EPG between individuals in each community}

\item{cv_within}{Day-to-day variation in EPG within an individual}

\item{cv_slide}{Variation between faecal samples from the same individual and day}

\item{cv_reduction}{Variation in efficacy between individuals}

\item{overall_mean}{The overall mean (i.e. the mean of the distribution reflecting community means)}

\item{edt}{The egg detection threshold (24 EPG is standard for Kato-Katz)}
}
\value{
A data frame containing the simulated data
}
\description{
Simulate ERR data from a compound lognormal Poisson distribution
}
\examples{
clpDataSim(10, 10, c(0.05, 0.1, 0.15), c(800,1000), c(1.5, 1), 1, 1, 1)

}
