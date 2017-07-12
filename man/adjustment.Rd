\name{adjustment}
\alias{adjustment}
\title{
Compute adjustment sets for covariate adjustment.
}
\description{
This function is a wrapper for convenience to the function \code{adjustmentSet} from package \pkg{dagitty}. 
}
\usage{
adjustment(amat, amat.type, x, y, set.type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{
  adjacency matrix of type \code{amat.cpdag} or \code{amat.pag}.
}
  \item{amat.type}{
  string specifying the type of graph of the adjacency matrix amat. It can be a DAG (type="dag"), a CPDAG (type="cpdag") or a maximally oriented PDAG (type="pdag") from Meek (1995); then the type of adjacency matrix is assumed to be amat.cpdag. It can also be a MAG (type = "mag") or a PAG (type="pag"); then the type of the adjacency matrix is assumed to be amat.pag.
}
  \item{x}{
  (integer) position of variable x in the adjacency matrix.
}
  \item{y}{
  (integer) position of variable y in the adjacency matrix.
}
  \item{set.type}{
  string specifying the type of adjustment set that should be computed. It can be "minimal" ,"all" and "canonical". See Details explanations.
}
}
\details{
If \code{set.type} is "minimal", then only minimal sufficient adjustment sets are returned. If \code{set.type} is "all", all valid
adjustment sets are returned. If \code{set.type} is "canonical", a single adjustment set
is returned that consists of all (possible) ancestors of \code{x} and {y},
minus (possible) descendants of nodes on proper causal paths.  This canonical adjustment set is always valid if any valid set exists at all.
}
\value{
List of adjustment sets.
}
\references{
  E. Perkovic, J. Textor, M. Kalisch and M.H. Maathuis (2015). A
  Complete Generalized Adjustment Criterion. In \emph{Proceedings of UAI
    2015.}

    E. Perkovic, J. Textor, M. Kalisch and M.H. Maathuis (2017). Complete graphical characterization and construction of adjustment sets in Markov equivalence classes of ancestral graphs. To appear in \emph{Journal of Machine Learning Research.}

    B. van der Zander, M. Liskiewicz and J. Textor (2014). Constructing
  separators and adjustment sets in ancestral graphs. In
  \emph{Proceedings of UAI 2014.}
}
\author{
Emilija Perkovic and Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\seealso{
\code{\link{gac}} for testing if a set satisfies the Generalized Adjustment Criterion.
}
\examples{
## Example 4.1 in Perkovic et. al (2015), Example 2 in Perkovic et. al (2017)
mFig1 <- matrix(c(0,1,1,0,0,0, 1,0,1,1,1,0, 0,0,0,0,0,1,
                  0,1,1,0,1,1, 0,1,0,1,0,1, 0,0,0,0,0,0), 6,6)
type <- "cpdag"
x <- 3; y <- 6
## plot(as(t(mFig1), "graphNEL"))

## all
adjustment(amat = mFig1, amat.type = type, x = x, y = y, set.type = "all")
## finds adjustment sets: (2,4), (1,2,4), (4,5), (1,4,5), (2,4,5), (1,2,4,5)

## minimal
adjustment(amat = mFig1, amat.type = type, x = x, y = y, set.type = "minimal")
## finds adjustment sets: (2,4), (4,5), i.e., the valid sets with the fewest elements

## canonical
adjustment(amat = mFig1, amat.type = type, x = x, y = y, set.type = "canonical")
## finds adjustment set: (1,2,4,5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{graphs}% __ONLY ONE__ keyword per line
