\name{rndLambdaF}
\alias{rndLambdaF}
\title{
  Principal Component Statistics Based on Randomization
}
\description{
  Implements randomization-based procedures to estimate the number of
  principal components.
}
\usage{
rndLambdaF(data, B = 1000, alpha = 0.05)
}
\arguments{
  \item{data}{A numeric data matrix.}
  \item{B}{An integer; the number of times to scramble the data columns.}
  \item{alpha}{A real number between 0 and 1; the significance level.}
}
\details{
  The randomization procedures implemented here were first developed by
  ter Brack [1,2].  In a simulation study, Peres-Neto and colleagues
  concluded that these methods were among the best [3]. Our own
  simulations on larger data matrices find that rnd-Lambda performs
  well (comparably to Auer-Gervini, though slower), but that rnd-F works
  poorly.

  The test procedure is: (1) randomize the values with all the attribute
  columns of the data matrix; (2) perform PCA on the scrambled data
  matrix; and (3) compute the test statistics. All three steps are
  repeated a total of (B - 1) times, where B is large enough to
  guarantee accuracy when estimating p-values; in practice, B is usually
  set to 1000. In each randomization, two test statistics are computed:
  (1) the eigenvalue \eqn{\lambda_k} for the k-th principal component; and
  (2) a pseudo F-ratio computed as \eqn{\lambda_k / \sum_{i=k+1}^n \lambda_i}.
  Finally, the p-value for each k and each statistic of
  interest is estimated to be the proportion of the test statistics in
  all data sets that are greater than or equal to the one in the
  observed data matrix.
}

\value{A named vector of length two, containing the predicted number of
  principal components based on the rnd-Lambda and rnd-F statistics. }

\references{
  [1] ter Braak CFJ. CANOCO -- a Fortran program for canonical community
  ordination by [partial] [detrended] [canonical] correspondence
  analysis, principal component analysis and redundancy analysis
  (version 2.1). Agricultural Mathematics Group, Report LWA-88- 02,
  Wageningen, 1988.

  [2] ter Braak CFJ. Update notes: CANOCO (version 3.1). Agricultural
  Mathematics Group, Wageningen, 1990.

  [3] Peres-Neto PR, Jackson DA and Somers KM. How many principal components? Stopping
  rules for determining the number of non-trivial axes revisited. Computational Statistics and
  Data Analysis 2005; 49: 974--997.
}
\author{
  Kevin R. Coombes <krc@silicovore.com>, Min Wang <wang.1807@osu.edu>.
}
\seealso{
  \link{AuerGervini-class}
}
\examples{
dataset <- matrix(rnorm(200*15, 6), ncol=15)
rndLambdaF(dataset)
}
\keyword{models}
