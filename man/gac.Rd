\name{gac}
\alias{gac}
\title{Test If Set Satisfies Generalized Adjustment Criterion (GAC)}
\description{
  This function tests if \code{z} satisfies  the Generalized Adjustment
  Criterion (GAC) relative to \code{(x,y)} in the graph represented by
  adjacency matrix \code{amat} and interpreted as \code{type} (DAG,
  maximal PDAG, CPDAG, MAG, PAG). If yes, \code{z} can be used in covariate adjustment
  for estimating causal effects of \code{x} on \code{y}.
}
\usage{
gac(amat, x, y, z, type = "pag")
}
\arguments{
  \item{amat}{adjacency matrix of type \link{amat.cpdag} or \link{amat.pag}}
  \item{x,y,z}{(integer) positions of variables in \code{x},
    \code{y} or \code{z} in the adjacency matrix. \code{x}, \code{y} and
  \code{z} can be vectors representing several nodes.}
  \item{type}{string specifying the type of graph of the adjacency matrix
    \code{amat}.  It can be a DAG (type="dag"), a PDAG (type="pdag") or a CPDAG (type="cpdag");
    then the type of the adjacency matrix is assumed to be
    \link{amat.cpdag}. It can also be a MAG (type="mag"), or a PAG
    (type="pag"); then the type of the adjacency matrix is assumed to be
  \link{amat.pag}.}
}
\details{
  This work is a generalization of the work of Shpitser et al. (2012)
  (necessary and sufficient criterion in DAGs/ADMGs)
  and van der Zander et al. (2014) (necessary and sufficient criterion
  in MAGs). Moreover, it is a generalization of the
  Generalized Backdoor Criterion (GBC) of Maathuis and Colombo (2013):
  While GBC is sufficient but not necessary, GAC is both sufficient and
  necessary for DAGs, CPDAGs, MAGs and PAGs. For more details see
  Perkovic et al. (2015, 2017a, 2017b).

  The motivation to find a set \code{z} that satisfies the GAC with
  respect to \code{(x,y)} is the following:

  \emph{A set of variables \code{z} satisfies the GAC relative to
    \code{(x,y)} in the given graph, if and only if
    the causal effect of \code{x} on \code{y} is identifiable by
    covariate adjustment and is given by}
  \deqn{%
    P(Y|do(X = x)) = \sum_Z P(Y|X,Z) \cdot P(Z),}{%
    P(Y|do(X = x)) =  sum_Z P(Y|X,Z) * P(Z),}
  (for any joint distribution \dQuote{compatible} with the graph; the formula
    is for discrete variables with straightforward modifications for
    continuous variables).
  This result allows to write post-intervention densities (the one
  written using Pearl's do-calculus) using only observational densities
  estimated from the data.

  For \code{z} to satisfy the GAC relative to \code{(x,y)} and the graph, the
  following three conditions must hold:
  \describe{
    \item{(0)}{ The graph is adjustment amenable relative to \code{(x,y)}.}
    \item{(1)}{ The intersection of \code{z} and the forbidden set
      (explained in Perkovic et al. (2015, 2017b) is empty.}
    \item{(2)}{ All proper definite status non-causal paths in the graph from
      \code{x} to \code{y} are blocked by \code{z}.}
  }

  It is important to note that there can be \code{x} and
  \code{y} for which there is no set Z that satisfies the GAC, but the
  total causal effect might be identifiable via some technique other
  than covariate adjustment.
  
  For details on the GAC for DAGs, CPDAGs, PAGs see Perkovic et. al (2015,2017a). For details on the GAC for MAGs see van der Zander et. al (2014) and for details on the GAC for maximal PDAGs see Perkovic et. al (2017b).

  For the coding of the adjacency matrix see \link{amatType}.
}
\value{
  A \code{\link{list}} with three components:
  \item{gac}{logical; TRUE if \code{z} satisfies the GAC relative to
    \code{(x,y)} in the graph represented by \code{amat} and \code{type}}
  \item{res}{logical vector of length three indicating if each of the three
    conditions (0), (1) and (2) are true}
  \item{f}{node positions of nodes in the forbidden set (see Perkovic
    et al. (2015, 2017b)}
}
\references{
  E. Perkovic, J. Textor, M. Kalisch and M.H. Maathuis (2015). A
  Complete Generalized Adjustment Criterion. In \emph{Proceedings of UAI
    2015.}

    E. Perkovic, J. Textor, M. Kalisch and M.H. Maathuis (2017a). Complete graphical characterization and construction of adjustment sets in Markov equivalence classes of ancestral graphs. To appear in \emph{Journal of Machine Learning Research.}

  E. Perkovic, M. Kalisch and M.H. Maathuis (2017b). Interpreting and using CPDAGs with background knowledge. In \emph{Proceedings of UAI 2017.}

  I. Shpitser, T. VanderWeele and J.M. Robins (2012). On the validity of
  covariate adjustment for estimating causal effects. In
  \emph{Proceedings of UAI 2010.}

  B. van der Zander, M. Liskiewicz and J. Textor (2014). Constructing
  separators and adjustment sets in ancestral graphs. In
  \emph{Proceedings of UAI 2014.}

  M.H. Maathuis and D. Colombo (2013). A generalized backdoor
  criterion. \emph{Annals of Statistics} 43 1060-1088.
}
\author{Emilija Perkovic and Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{backdoor}} for the Generalized Backdoor Criterion,
  \code{\link{pc}} for estimating a CPDAG and \code{\link{fci}} and
  \code{\link{fciPlus}} for estimating a PAG.
}
\examples{%% note:  Tests in  ../tests/test_gac.R
## We reproduce the four examples in Perkovic et. al (2015, 2017a)

##############################
## Example 4.1 in Perkovic et. al (2015), Example 2 in Perkovic et. al (2017a)
##############################
mFig1 <- matrix(c(0,1,1,0,0,0, 1,0,1,1,1,0, 0,0,0,0,0,1,
                  0,1,1,0,1,1, 0,1,0,1,0,1, 0,0,0,0,0,0), 6,6)
type <- "cpdag"
x <- 3; y <- 6
## Z satisfies GAC :
gac(mFig1, x,y, z=c(2,4),    type)
gac(mFig1, x,y, z=c(4,5),    type)
gac(mFig1, x,y, z=c(4,2,1),  type)
gac(mFig1, x,y, z=c(4,5,1),  type)
gac(mFig1, x,y, z=c(4,2,5),  type)
gac(mFig1, x,y, z=c(4,2,5,1),type)
## Z does not satisfy GAC :
gac(mFig1,x,y, z=2,    type)
gac(mFig1,x,y, z=NULL, type)

##############################
## Example 4.2 in Perkovic et. al (2015), Example 3 in Perkovic et. al (2017a)
##############################
mFig3a <- matrix(c(0,1,0,0, 1,0,1,1, 0,1,0,1, 0,1,1,0), 4,4)
mFig3b <- matrix(c(0,2,0,0, 3,0,3,3, 0,2,0,3, 0,2,2,0), 4,4)
mFig3c <- matrix(c(0,3,0,0, 2,0,3,3, 0,2,0,3, 0,2,2,0), 4,4)
type <- "pag"
x <- 2; y <- 4
## Z does not satisfy GAC
gac(mFig3a,x,y, z=NULL, type) ## not amenable rel. to (X,Y)
gac(mFig3b,x,y, z=NULL, type) ## not amenable rel. to (X,Y)
## Z satisfies GAC
gac(mFig3c,x,y, z=NULL, type) ## amenable rel. to (X,Y)

##############################
## Example 4.3 in Perkovic et. al (2015), Example 4 in Perkovic et. al (2017a)
##############################
mFig4a <- matrix(c(0,0,1,0,0,0, 0,0,1,0,0,0, 2,2,0,3,3,2,
                   0,0,2,0,2,2, 0,0,2,1,0,2, 0,0,1,3,3,0), 6,6)
mFig4b <- matrix(c(0,0,1,0,0,0, 0,0,1,0,0,0, 2,2,0,0,3,2,
                   0,0,0,0,2,2, 0,0,2,3,0,2, 0,0,2,3,2,0), 6,6)
type <- "pag"
x <- 3; y <- 4
## both PAGs are amenable rel. to (X,Y)
## Z satisfies GAC in Fig. 4a
gac(mFig4a,x,y, z=6, type)
gac(mFig4a,x,y, z=c(1,6), type)
gac(mFig4a,x,y, z=c(2,6), type)
gac(mFig4a,x,y, z=c(1,2,6), type)
## no Z satisfies GAC in Fig. 4b
gac(mFig4b,x,y, z=NULL, type)
gac(mFig4b,x,y, z=6, type)
gac(mFig4b,x,y, z=c(5,6), type)

##############################
## Example 4.4 in Perkovic et. al (2015), Example 8 in Perkovic et. al (2017a)
##############################
mFig5a <- matrix(c(0,1,0,0,0, 1,0,1,0,0, 0,0,0,0,1, 0,0,1,0,0, 0,0,0,0,0), 5,5)
type <- "cpdag"
x <- c(1,5); y <- 4
## Z satisfies GAC
gac(mFig5a,x,y, z=c(2,3), type)
## Z does not satisfy GAC
gac(mFig5a,x,y, z=2, type)

mFig5b <- matrix(c(0,1,0,0,0,0,0, 2,0,2,3,0,3,0, 0,1,0,0,0,0,0,
0,2,0,0,3,0,0, 0,0,0,2,0,2,3, 0,2,0,0,2,0,0, 0,0,0,0,2,0,0), 7,7)
type <- "pag"
x<-c(2,7); y<-6
## Z satisfies GAC
gac(mFig5b,x,y, z=c(4,5), type)
gac(mFig5b,x,y, z=c(4,5,1), type)
gac(mFig5b,x,y, z=c(4,5,3), type)
gac(mFig5b,x,y, z=c(1,3,4,5), type)
## Z does not satisfy GAC
gac(mFig5b,x,y, z=NULL, type)

##############################
## Example 4.7 in Perkovic et. al (2017b)
##############################
mFig3a <- matrix(c(0,1,0,0, 1,0,1,1, 0,1,0,1, 0,1,1,0), 4,4)
mFig3b <- matrix(c(0,1,0,0, 0,0,1,1, 0,0,0,1, 0,0,1,0), 4,4)
mFig3c <- matrix(c(0,0,0,0, 1,0,1,0, 0,1,0,1, 0,1,1,0), 4,4)
type <- "pdag"
x <- 2; y <- 4
## Z does not satisfy GAC
gac(mFig3a,x,y, z=NULL, type) ## not amenable rel. to (X,Y)
gac(mFig3c,x,y, z=NULL, type) ## amenable rel. to (X,Y), but no set can block X <- Y
## Z satisfies GAC
gac(mFig3b,x,y, z=NULL, type) ## amenable rel. to (X,Y)


}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
