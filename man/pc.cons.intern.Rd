\name{pc.cons.intern}
\alias{pc.cons.intern}
\title{Internal function for conservative PC algorithm}
\description{
  This function is not intended for direct useage.  Use \code{\link{pc}}
  or \code{\link{fci}} using options for the conservative version of the
  algorithms instead.

For any unshielded triple A-B-C, consider all subsets D of the neighbors
of A and of the neighbors of C, and record the sets D for which A and C
are conditionally independent given D. If B is in none of these sets, do
nothing (it is a v-structure). If B is in all sets, do nothing (it is
not a v-structure). If B is in some but not all sets, mark the triple as
\dQuote{unfaithful}.
}

\usage{
pc.cons.intern(sk, suffStat, indepTest, alpha, verbose = FALSE,
               version.unf = c(NA, NA))
}

\arguments{
  \item{sk}{an object as returned from \code{\link{skeleton}()}.}
  \item{suffStat}{sufficient statistic: List containing all necessary
    elements for the conditional independence decisions in the
    function \code{indepTest}.}
  \item{indepTest}{predefined function for testing conditional independence. The
    function is internally called as \code{indepTest(x,y,S,suffStat)}, and
    tests conditional independence of \code{x} and \code{y} given
    \code{S}. Here, \code{x} and \code{y} are variables, and \code{S} is
    a (possibly empty) vector of variables (all variables are denoted
    by their column numbers
    in the adjacency matrix). \code{suffStat} is a list containing
    all relevant elements for the conditional independence
    decisions. The return value of \code{indepTest} is the p-value of
    the test for conditional independence.}
  \item{alpha}{significance level for the individual conditional
    independence tests.}
  \item{verbose}{If \code{TRUE}, detailed output is provided.}
  \item{version.unf}{vector of length two. First argument: Consider the
    case, where a indep c given S in the skeleton; furthermore, suppose
    that a and c are dependent given every subset of the neigbors in the
    conservative step; then, some error must have occured in the pc, but
    not necessarily in the fci (since here sepsets can contain nodes
    outside the neighborhood). If this option is 1 the corresponding triple
    is marked 'faithful'; if this option is 2, the corresponding triple is
    marked unfaithful. Second argument: 1 do not consider the initial sepset,
    2 also consider the initial sepset.
  }
}

\value{
  \item{unfTripl}{Vector with triples (coded as number using
    \code{triple2numb}) that were marked as unfaithful.}

  \item{vers}{Vector containing the version (1 or 2) of the
    corresponding triple saved in unfTripl (1=normal unfaithful triple that
    is, B is in some sepsets but not all or none; 2=triple coming from version.unf[1]==2 that is
    a and c are indep given the initial sepset but there doesn't exist a
    subset of the neighbours that d-separates them.)}
}

\seealso{\code{\link{skeleton}} for estimating a skeleton
  using the PC algorithm; \code{\link{pc}} for estimating a CPDAG using
  the PC algorithm; \code{\link{fci}} for computing
  the FCI algorithm.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Diego Colombo.
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
