\name{fciPlus}
\alias{fciPlus}
\title{Estimate a PAG, using the FCI+ algorithm}
\description{
  Estimate a Partial Ancestral Graph (PAG) from observational data, using
  the FCI+ Algorithm.
}
\usage{
fciPlus(suffStat, indepTest, alpha, labels, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{suffStat}{
    Sufficient statistics: List containing all necessary
    elements for the conditional independence decisions in the
    function \code{indepTest}.
  }
  \item{indepTest}{
    Predefined function for testing conditional independence. The
    function is internally called as \code{indepTest(x,y,S,suffStat)}, and
    tests conditional independence of \code{x} and \code{y} given
    \code{S}. Here, \code{x} and \code{y} are variables, and \code{S} is
    a (possibly empty) vector of variables (all variables are denoted
    by their column numbers
    in the adjacency matrix). \code{suffStat} is a list containing
    all relevant elements for the conditional independence
    decisions. The return value of \code{indepTest} is the p-value of
    the test for conditional independence.
  }
  \item{alpha}{
    Significance level for the individual conditional
    independence tests.
  }
  \item{labels}{
    Character vector of variable (or
    \dQuote{node}) names.  Typically preferred to specifying \code{p}.}
  \item{p}{(optional) number of variables (or nodes).  May be specified
    if \code{labels} are not, in which case \code{labels} is set to
    \code{1:p}.
  }
}
\details{
  A variation of FCI. For details, please see the references.
}
\value{
An object of \code{\link{class}} \code{fciAlgo} (see
 \code{\linkS4class{fciAlgo}}) containing the estimated graph
 (in the form of an adjacency matrix with various possible edge marks),
 the conditioning sets that lead to edge removals (sepset) and several other
 parameters.
}
\references{
  T. Claassen, J. Mooij, and T. Heskes (2013). Learning Sparse Causal Models is not NP-hard. In \emph{UAI 2013, Proceedings of the 29th Conference on Uncertainty in Artificial Intelligence}
}
\author{
  Emilija Perkovic and Markus Kalisch (\email{kalisch@stat.math.ethz.ch}).
}
\seealso{
  \code{\link{fci}} for estimating a PAG using the FCI algorithm.
}
\examples{
##################################################
## Example without latent variables
##################################################

set.seed(42)
p <- 7
## generate and draw random DAG :
myDAG <- randomDAG(p, prob = 0.4)

## find skeleton and PAG using the FCI algorithm
suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
res <- fciPlus(suffStat, indepTest=gaussCItest,
           alpha = 0.9999, p=p)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

