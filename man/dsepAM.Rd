\name{dsepAM}
\alias{dsepAM}
\title{Test for d-separation in a MAG}
\description{
  This function tests for d-separation (also known as m-separation) of nodes \code{X} and nodes \code{Y} given nodes \code{S} in a MAG.
}
\usage{
dsepAM(X, Y, S = NULL, amat, verbose=FALSE)
}
\arguments{
  \item{X}{Vector of column numbers of nodes \code{X} in the adjacency matrix}
  \item{Y}{Vector of column numbers of nodes \code{Y} in the adjacency matrix}
  \item{S}{Vector of column numbers of nodes \code{S} in the adjacency matrix (may be empty)}
  \item{amat}{The Maximal Ancestral Graph encoded as adjacency matrix 
  	  of type \link{amatType}}
  \item{verbose}{If true, more detailed output is provided.}
}
\value{TRUE if \code{X} and \code{Y} are d-separated by \code{S} in the 
       MAG encoded by \code{amat}, otherwise FALSE.
}
\details{This function checks separation in the moralized graph as
  explained in Richardson and Spirtes (2002).}
\references{
  T.S. Richardson and P. Spirtes (2002). Ancestral graph Markov models.
  \emph{Annals of Statistics} \bold{30} 962-1030.
}
\seealso{\code{\link{dsepAMTest}} for a wrapper of this function that can
  easily be included into \code{\link{skeleton}}, \code{\link{fci}} or
  \code{\link{fciPlus}}. \code{\link{dsep}} for a similar function for DAGs.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}), Joris Mooij
}
\keyword{graphs}
\examples{
# Y-structure MAG
# Encode as adjacency matrix
p <- 4 # total number of variables
V <- c("X1","X2","X3","X4") # variable labels
# amat[i,j] = 0 iff no edge btw i,j
# amat[i,j] = 1 iff i *-o j
# amat[i,j] = 2 iff i *-> j
# amat[i,j] = 3 iff i *-- j
amat <- rbind(c(0,0,2,0),
              c(0,0,2,0),
              c(3,3,0,2),
              c(0,0,3,0))
rownames(amat)<-V
colnames(amat)<-V

## d-separated
cat('X1 d-separated from X2? ', dsepAM(1,2,S=NULL,amat),'\n')
## not d-separated given node 3
cat('X1 d-separated from X2 given X4? ', dsepAM(1,2,S=4,amat),'\n')
## not d-separated by node 3 and 4
cat('X1 d-separated from X2 given X3 and X4? ', dsepAM(1,2,S=c(3,4),amat),'\n') 
}
