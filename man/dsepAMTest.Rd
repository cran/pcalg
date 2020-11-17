\name{dsepAMTest}
\alias{dsepAMTest}
\title{Test for d-separation in a MAG}
\description{
  This function tests for d-separation (also known as m-separation) of node \code{x} and node \code{y} given nodes \code{S} in a MAG.
  \code{dsepAMTest()} is written to be easily used in \code{\link{skeleton}}, \code{\link{fci}},
  \code{\link{fciPlus}}.
}
\usage{
dsepAMTest(x, y, S = NULL, suffStat)
}
\arguments{
  \item{x}{Column number of node \code{x} in the adjacency matrix}
  \item{y}{Column number of node \code{y} in the adjacency matrix}
  \item{S}{Vector of column numbers of nodes \code{S} in the adjacency
 	   matrix, may be empty}
  \item{suffStat}{a \code{\link{list}} with two elements,
    \describe{
      \item{amat}{The Maximal Ancestral Graph encoded as adjacency matrix 
                  of type \link{amatType}}
      \item{verbose}{If true, more detailed output is provided.}
      }}
}
\value{Returns 1 if \code{x} and \code{y} are d-separated by \code{S} in the 
       MAG encoded by \code{amat}, otherwise 0.
	
  This is analogous to the p-value of an ideal
  (without sampling error) conditional independence test on any
  distribution that is faithful to the MAG.
}
\details{The function is a wrapper for \code{\link{dsepAM}}, which checks
	separation in the moralized graph as explained in Richardson and Spirtes (2002).
}
\references{
  T.S. Richardson and P. Spirtes (2002). Ancestral graph Markov models.
  \emph{Annals of Statistics} \bold{30} 962-1030.
}
\seealso{\code{\link{dsepTest}} for a similar function for DAGs.
  \code{\link{gaussCItest}}, \code{\link{disCItest}} and
  \code{\link{binCItest}} for similar functions for a conditional
  independence test for gaussian, discrete and
  binary variables, respectively.
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

suffStat<-list(g=amat,verbose=FALSE)
## d-separated
cat('X1 d-separated from X2? ', dsepAMTest(1,2,S=NULL,suffStat),'\n')
## not d-separated given node 3
cat('X1 d-separated from X2 given X4? ',
dsepAMTest(1,2,S=4,suffStat),'\n')
## not d-separated by node 3 and 4
cat('X1 d-separated from X2 given X3 and X4? ', dsepAMTest(1,2,S=c(3,4),
suffStat),'\n')

# Derive PAG that represents the Markov equivalence class of the MAG with the FCI algorithm
# Make use of d-separation oracle as "independence test"
indepTest <- dsepAMTest
fci.pag <- fci(suffStat,indepTest,alpha = 0.5,labels = V,verbose=FALSE)
cat('True MAG:\n')
print(amat)
cat('PAG output by FCI:\n')
print(fci.pag@amat)
}
