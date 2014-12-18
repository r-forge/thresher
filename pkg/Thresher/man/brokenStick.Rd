\name{brokenStick}
\alias{brokenStick}
\alias{bsDimension}
\title{
  The Broken Stick Method
}
\description{
  The Broken Stick model is one proposed method for estimating the
  number of statistically significant principal components. 
}
\usage{
brokenStick(k, n)
bsDimension(lambda, FUZZ = 0.005)
}
\arguments{
  \item{k}{}
  \item{n}{}
  \item{lambda}{}
  \item{FUZZ}{}
}
\details{
  The Broken Stick model is one proposed method for estimating the
  number of statistically significant principal components.  The idea is
  to model $N$ variances by taking a stick of unit length and breaking it
  into $N$ pieces by randomly slecting break points from a uniform
  distribution. 

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Jackson, D. A. (1993). Stopping rules in principal components analysis: a comparison of heuristical and statistical approaches. Ecology 74, 2204--2214.

Legendre, P. and Legendre, L. (1998) Numerical Ecology. 2nd English ed. Elsevier.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
\keyword{ models }
\keyword{ cluster }
