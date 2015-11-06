\name{agDimFunction}
\alias{agDimFunction}
\alias{agDimTwiceMean}
\alias{agDimKmeans}
\alias{agDimKmeans3}
\alias{agDimLeap}
\alias{agDimSpectral}
\alias{agDimTtest}
\alias{agDimTtest2}
\alias{agDimCPT}
\alias{makeAgCpmFun}
\title{
  Divide Steps into "Long" and "Short" to Compute Auer-Gervini Dimension
}
\description{
  Auer and Gervini developed a Bayesian graphical method to determine
  the number \eqn{d} of significant principal components; a brief
  overview is included in the help for the \code{\link{AuerGervini}}
  class.  The output of their method is a step function that displays
  the maximum a posteriori (MAP) choice of \eqn{d} as a step function of
  a one-parameter family of prior distributions, and thbey recommend
  choosing the higest "long" step.  The functions described here help
  automate the process of dividing the step lengths into "long" and
  "short" classes.
}
\usage{
  agDimTwiceMean(stepLength)
  agDimKmeans(stepLength)
  agDimKmeans3(stepLength)
  agDimLeap(stepLength)
  agDimSpectral(stepLength)
  agDimTtest(stepLength, extra=0)
  agDimTtest2(stepLength)
  agDimCPT(stepLength)
  makeAgCpmFun(method)
}
\arguments{
  \item{stepLength}{A numeric vector}
  \item{method}{A character string describing a method supported by the
    \code{\link{detectChangePointBatch}} function in the \code{cpm}
    package.}
  \item{extra}{Just ignore this. Don't use it. It's a hack to avoid
    having to maintain two differfent versions of the same code.}
}
\details{
  The \code{agDimTwiceMean} function implements a simple and naive rule:
  a step is considered long if it as least twice the mean length.

  The \code{agDimKmeans} uses the \code{\link{kmeans}} algorithm with
  \eqn{k=2} to divide the step lengths into two classes.  Starting
  centers for the groups are taken to be the minimum and maximum
  values.

  The \code{agDimKmeans3} function uses \code{kmeans} with \eqn{k=3},
  using the median as the third center.  Only one of the three groups is
  considered "short".

  The \code{agDimSpectral} applies spectral clustering (as implemented
  by the \code{\link{specc}} function from the \code{kernlab} package)
  to divide the steps lengths into two groups.

  The \code{agDimLeap}, \code{agDimTtest}, and \code{agDimTtest2}
  functions implement two variants of a novel algorithm specialized for
  this particular task.  The idea is to start by sorting the step
  lengths so that \deqn{L_1 \le L_2 \le \dots \le L_n.}  Then, for each
  \eqn{i \in 3,\dots, N-1}, we compute the mean and standard deviation
  of the first \eqn{i} step lengths.  Finally, one computes the
  likelhood that \eqn{L_{i+1}} comes from the normal distribution
  defined by the first \eqn{i} lengths. If the probability that
  \eqn{L_{i+1}} is larger is less than \eqn{0.01}, then it is chosen as
  the "smallest long step".

  The novel method just described can also be viewed as a way to detect
  a particular kind of change point.  So, we also provide the
  \code{agDimCPT} function that uses the changepoint detection
  algorithm implement by the \code{\link{cpt.mean}} function in the
  \code{changepoint} package.  More generally, the \code{makeAgCpmFun}
  allows you to use any of the changepoint models implemented as part
  of the \code{\link{detectChangePointBatch}} function in the \code{cpm}
  package.
}
\value{
  Each of the functions \code{agDimTwiceMean}, \code{agDimKmeans},
  \code{agDimKmeans3}, \code{agDimSpectral}, \code{agDimTtest},
  \code{agDimTtest2}, and \code{agDimCPT} returns a logical vector whose
  length is equal to the input \code{stepLength}.  \code{TRUE} values
  identify "long" steps and \code{FALSE} values identify "short"
  steps.

  The \code{makeAgCpmFun} returns a function that takes one argument (a
  numeric \code{stepLength} vector) and returns a logical vector of the
  same length.
}
\references{
  P Auer, D Gervini.
  Choosing principal components: a new graphical method based on Bayesian model selection.
  Communications in Statistics-Simulation and Computation 37 (5), 962-977
}
\author{
  Kevin R. Coombes <krc@silicovore.com>, Min Wang <wang.1807@osu.edu>.
}
\seealso{
  The functions described here implerment different algorithms that can
  be used by the \code{\link{agDimension}} function to automatically
  compute the number of significant principal components based on the
  \code{\link{AuerGervini}} approach.  Several of these functions are
  wrappers around functions defined in other packages, including
  \code{\link{specc}} in the \code{kernlab} package,
  \code{\link{cpt.mean}} in the \code{changepoint} package, and
  \code{\link{detectChangePointBatch}} in the \code{cpm} package.
}
\examples{
# simulate variances
lambda <- rev(sort(diff(sort(c(0, 1, runif(9))))))
# apply the Auer-Gervini method
ag <- AuerGervini(lambda, dd=c(3,10))
# Review the results
summary(ag)
agDimension(ag)
agDimension(ag, agDimKmeans)
agDimension(ag, agDimSpectral)
f <- makeAgCpmFun("Exponential")
agDimension(ag, f)
}
\keyword{models}
