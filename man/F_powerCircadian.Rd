% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/F_powerCircadian.R
\name{F_powerCircadian}
\alias{F_powerCircadian}
\title{F-statistics-based Theoretical Circadian Power Calculation}
\usage{
F_powerCircadian(n, A, sigma, phi = 0, period = 24, cts = NULL, alpha = 0.05)
}
\arguments{
\item{n}{Sample size.}

\item{A}{Amplitude of the sine curve: \eqn{A * sin(2\pi/period * (phi + cts)) + C}.}

\item{sigma}{\sigma in the independent Normal error term \eqn{N(0,\sigma)}.}

\item{phi}{Phase of the sine curve. Default is 0.}

\item{period}{Period of the sine curve. Default is 24.}

\item{cts}{Circadian time design vector.}

\item{alpha}{Type I error control. Default is 0.05.}
}
\value{
Theoretucal power based on F-statistics.
}
\description{
Theoretical circadian power calculation based on F-statistics
}
\details{
Theoretically calculate power of circaidan data under indepednent Normal error assumption using F-statistics.
}
\examples{
n = 100
cts = seq(0,24,length.out = n+1)[-1]
A = 1
sigma = 1
phi = 0
C = 10
F_powerCircadian(n, A, sigma, phi=0, period = 24, cts=cts, alpha = 0.05)
}
\author{
Wei Zong, Zhiguang Huo
}
