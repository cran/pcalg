\name{searchAM}
\alias{searchAM}
\title{Search for certain nodes in a DAG/CPDAG/MAG/PAG}
\description{
  Searches for all ancestors, descendants, anteriors, spouses, neighbors,
  parents, children or possible descendants of a (set of) node(s) 
  in a DAG, CPDAG, MAG or PAG.
}
\usage{
searchAM(amat,x,
	 type = c("an", "de", "ant", "sp", "nb", "pa", "ch", "pde"))
}
\arguments{
	\item{amat}{Adjacency matrix of type \link{amat.pag}.}
	\item{x}{Target node(s), given as (a vector of) column number(s) 
		of the node(s) in the adjacency matrix.}
	\item{type}{Character string specifying which relation to the
		target nodes in \code{x} is asked for. It can be one of:
    \describe{
      \item{"an"}{ancestors: nodes \code{y} s.t. \code{y-->...-->x};}
      \item{"de"}{descendants: nodes \code{y} s.t. \code{y<--...<--x};}
      \item{"ant"}{anteriors: nodes \code{y} s.t. \code{y---...---z-->...-->x},
                   i.e. there is an undirected path from \code{y} to a node \code{z}
                   followed by a directed path from \code{z} to \code{x};}
      \item{"sp"}{spouses: nodes \code{y} s.t. \code{y<->x};}
      \item{"nb"}{neighbors: nodes \code{y} s.t. \code{y---x};}
      \item{"pa"}{parents: nodes \code{y} s.t. \code{y-->x};}
      \item{"ch"}{children: nodes \code{y} s.t. \code{y<--x};}
      \item{"pde"}{possible descendants: nodes \code{y} s.t. there is a possibly
                   directed path from \code{y} to \code{x}: \code{y {o,-}--{o,>} ... {o,-}--{o,>} x}.}
    }
	For the precise definitions of these concepts, see the references.}
}
\value{Vector of column numbers of the nodes related to \code{x} as specified by \code{type}.
}
\details{
  This function performs a search for nodes related to the set of target nodes \code{x} 
  in the way specified by \code{type} in adjacency matrix \code{amat} of type \link{amat.pag}.
}
\references{
  T.S. Richardson and P. Spirtes (2002). Ancestral graph Markov models.
  \emph{Annals of Statistics} \bold{30} 962-1030.

  J. Zhang (2008). Causal Reasoning with Ancestral Graphs.
  \emph{Journal of Machine Learning Research} \bold{9} 1437-1474.
}
\author{
  Joris Mooij.
}
\examples{
# Y-structure MAG
# Encode as adjacency matrix
p <- 10 # total number of variables
V <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10") # variable labels
# amat[i,j] = 0 iff no edge btw i,j
# amat[i,j] = 1 iff i *-o j
# amat[i,j] = 2 iff i *-> j
# amat[i,j] = 3 iff i *-- j
amat <- rbind(c(0,3,0,0,0,0,0,0,0,0),
              c(3,0,3,0,0,0,0,0,0,0),
              c(0,3,0,2,0,0,0,0,0,0),
              c(0,0,3,0,2,0,0,0,0,0),
              c(0,0,0,3,0,2,0,2,2,1),
              c(0,0,0,0,3,0,2,0,0,0),
              c(0,0,0,0,0,3,0,0,0,0),
              c(0,0,0,0,2,0,0,0,0,0),
              c(0,0,0,0,1,0,0,0,0,0),
              c(0,0,0,0,1,0,0,0,0,0))
rownames(amat)<-V
colnames(amat)<-V

stopifnot(all.equal(searchAM(amat,5,type = "an"), c(3,4,5))) # ancestors of X5
stopifnot(all.equal(searchAM(amat,5,type = "de"), c(5,6,7))) # descendants of X5
stopifnot(all.equal(searchAM(amat,5,type = "ant"), c(1,2,3,4,5))) # anteriors of X5
stopifnot(all.equal(searchAM(amat,5,type = "sp"), c(8))) # spouses of X5
stopifnot(all.equal(searchAM(amat,2,type = "nb"), c(1,3))) # neighbors of X2
stopifnot(all.equal(searchAM(amat,c(4,6),type = "pa"), c(3,5))) # parents of {X4,X6}
stopifnot(all.equal(searchAM(amat,c(3,5),type = "ch"), c(4,6))) # children of {X3,X5}
stopifnot(all.equal(searchAM(amat,5,type = "pde"), c(5,6,7,9,10))) # possible descendants of X5
}
\keyword{graphs}
