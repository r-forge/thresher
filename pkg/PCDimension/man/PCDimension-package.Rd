\name{PCDimension-package}
\alias{PCDimension-package}
\alias{PCDimension}
\docType{package}
\title{
Finding the Number of Significant Principal Components
}
\description{
  The \code{PCDimension} package defines classes and methods to
  determine the number of statistically signifcant components in a data
  set. The main innovation includes methods to automate the graphical
  Bayesian method introduced by Auer and Gervini. However, the package
  also includes an implementaion of the broken stick method and of
  randomization-based methods introduced by ter Braak.
}
\details{
\tabular{ll}{
Package: \tab PCDimension\cr
Type: \tab Package\cr
Version: \tab 1.1.8\cr
Date: \tab 2018-01-09\cr
License: \tab Apache 2.0\cr
Depends: \tab ClassDiscovery\cr
Imports: \tab methods, stats, graphics, oompaBase, kernlab, changepoint, cpm\cr
}
Standard usage is to apply the \code{AuerGervini} function to a data set in
order to estimate the principal component dimension.
}
\author{
Kevin R. Coombes and Min Wang

Maintainer:  <krc@silicovore.com>
}
\references{
  [1] P Auer, D Gervini.
  Choosing principal components: a new graphical method based on Bayesian model selection.
  Communications in Statistics-Simulation and Computation 37 (5),
  962-977.

  [2] ter Braak CFJ. CANOCO -- a Fortran program for canonical community
  ordination by [partial] [detrended] [canonical] correspondence
  analysis, principal component analysis and redundancy analysis
  (version 2.1). Agricultural Mathematics Group, Report LWA-88- 02,
  Wageningen, 1988.

  [3] ter Braak CFJ. Update notes: CANOCO (version 3.1). Agricultural
  Mathematics Group, Wageningen, 1990.

  [4] Wang M, Kornbla SM, Coombes KR.
  Decomposing the Apoptosis Pathway Into Biologically Interpretable
  Principal Components.
  bioRxiv, 2017. <doi://10.1101/237883>.
}
\keyword{ package }
