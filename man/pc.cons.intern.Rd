\name{pc.cons.intern}
\alias{pc.cons.intern}
\alias{triple2numb}
\title{Utility for conservative and majority rule in PC and FCI}
\description{
  The \code{pc.cons.intern()} function is used in \code{\link{pc}} and 
  \code{\link{fci}}, notably when
  \code{conservative = TRUE} (conservative orientation of v-structures) or
  \code{maj.rule = TRUE} (majority rule orientation of v-structures). 
}
\usage{
pc.cons.intern(sk, suffStat, indepTest, alpha, version.unf = c(NA, NA),
               maj.rule = FALSE, verbose = FALSE)
%% triple2numb(p,i,j,k)
}
\arguments{
  \item{sk}{A skeleton object as returned from \code{\link{skeleton}()}.}
  \item{suffStat}{Sufficient statistic: List containing all necessary
    elements for the conditional independence decisions in the
    function \code{indepTest}.}
  \item{indepTest}{Pre-defined function for testing conditional independence. The
    function is internally called as \code{indepTest(x,y,S,suffStat)}, and
    tests conditional independence of \code{x} and \code{y} given
    \code{S}. Here, \code{x} and \code{y} are variables, and \code{S} is
    a (possibly empty) vector of variables (all variables are denoted
    by their column numbers
    in the adjacency matrix). \code{suffStat} is a list containing
    all relevant elements for the conditional independence
    decisions. The return value of \code{indepTest} is the p-value of
    the test for conditional independence.}
  \item{alpha}{Significance level for the individual conditional
    independence tests.}
  \item{version.unf}{Vector of length two.  If \code{version.unf}[2]==1,
    the intitial separating set found by the PC/FCI algorithm is
    added to the set of separating sets; if \code{version.unf[2]==}2, it
    is not added.  In the latter case, if the set of separating sets is
    empty, the triple is marked as unambiguous if
    \code{version.unf}[1]==1, and as ambiguous if
    \code{version.unf}[1]==2.}
  \item{maj.rule}{Logical indicatin if the triples are checked for
    ambiguity using the majority rule idea, which is less strict than
    the standard conservative method.}
  \item{verbose}{Logical asking for detailed output.}
}
\details{
  For any unshielded triple A-B-C, consider all subsets of the neighbors
  of A and of the neighbors of C, and record all such sets D for which A and C
  are conditionally independent given D.  We call such sets
  \dQuote{separating sets}. 

  If \code{version.unf}[2]==1, the initial separating set found in the
  PC/FCI algorithm is added to this set of separating sets.  
  If \code{version.unf}[2]==2, the initial separating set is not added (as in Tetrad). 

  In the latter case, if the set of separating sets is empty, then the
  triple is marked as \sQuote{ambiguous} if \code{version.unf}[1]==2, for
  example in \code{\link{pc}}, or as \sQuote{unambiguous} if
  \code{version.unf}[1]==1, for example in \code{\link{fci}}.% and \code{\link{rfci}}.  
  Otherwise, there is at least one separating set.
  If \code{maj.rule=FALSE}, the conservative PC algorithm is used
  (Ramsey et al., 2006): If B is in some but not all separating sets,
  the triple is marked as ambiguous.  Otherwise it is treated as in the
  standard PC algorithm.  If \code{maj.rule=TRUE}, the majority rule is
  applied (Colombo and Maathuis, 2014): The triple is marked as
  \sQuote{ambiguous} if B is in exactly 50 percent of the separating sets.  If
  it is in less than 50 percent it is marked as a v-structure, and if it
  is in more than 50 percent it is marked as a non v-structure.

  Note: This function modifies the separating sets for unambiguous
  triples in the skeleton object (adding or removing B) to ensure that
  the usual orientations rules later on lead to the correct
  v-structures/non v-structures. 
}
\value{
  \item{unfTripl}{numeric vector of triples coded as numbers (via
    \code{triple2numb()}) that were marked as ambiguous.}
  \item{vers}{
    Vector containing the version (1 or 2) of the corresponding triple
    saved in unfTripl (1=normal ambiguous triple, i.e., B is in some
    sepsets but not all or none; 2=triple coming from version.unf[1]==2,
    i.e., a and c are indep given the initial sepset but there does not
    exist a subset of the neighbours of a or of c that d-separates
    them.)}
  \item{sk}{The updated skeleton-object (separating sets might have been updated).}
}
\references{
  D. Colombo and M.H. Maathuis (2014).Order-independent constraint-based
  causal structure learning. \emph{Journal of Machine Learning Research}
  \bold{15} 3741-3782.

  J. Ramsey, J. Zhang and P. Spirtes (2006).
  Adjacency-faithfulness and conservative causal inference. In
  \emph{Proceedings of the 22nd Annual Conference on Uncertainty in
  Artificial Intelligence}, Arlington, VA. AUAI Press.
}
\seealso{\code{\link{skeleton}}, \code{\link{pc}}, \code{\link{fci}}
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Diego Colombo.
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
