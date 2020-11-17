\name{pag2edge}
\alias{pag2edge}
\title{Reads off identifiable parents and non-parents from a directed PAG}
\description{
  Constructs a matrix which contains identifiable parental and
  non-parental relations in the Markov equivalence class represented
  by a directed partial ancestral graph.
}
\usage{
pag2edge(P)
}
\arguments{
	\item{P}{Adjacency matrix of type \link{amat.pag}, which should encode
	a directed PAG (i.e., it should not contain any undirected edges of the
	form \code{x --- y} or any circle-tail edges of the form \code{x o-- y}).}
}
\value{Matrix \code{A}, where entry \code{A[i,j]} equals \describe{
      \item{1}{if node \code{i} is an identifiable parent of node \code{j},}
      \item{-1}{if node \code{i} is an identifiable non-parent of node \code{j},}
      \item{0}{in case the parental relationship between nodes \code{i} and \code{j} is unknown.}}
}
\details{
  We say that node \code{i} is parent of node \code{j} in a directed mixed graph (DMG) iff there
  exists a directed edge \code{i-->j} in that graph. If the directed mixed graph
  has a causal interpretation (for example, if it is the graph of a simple SCM) then parental
  relations coincide (generically) with direct causal relations (relative to the variables in the graph).

  This function implements the sufficient conditions (Propositions 7 and 8) in Mooij and Claassen (2020) 
  for concluding whether a parental relation between two nodes must be present or absent
  in all directed mixed graphs in the Markov equivalence class represented by the directed PAG \code{P}.
  It applies to both the
  acyclic case as well as the cyclic (simple SCM) case, assuming the d-separation resp. \eqn{\sigma}-separation
  Markov property.

  The output is a matrix containing for each ordered pair of nodes whether
  the presence of a parental relation was identified, or the absence, or neither.
  
  It is not known whether these sufficient conditions for identifiability are complete. 
  Hence, zero entries in the result indicate that the sufficient condition gives no 
  conclusion, rather than that the Markov equivalence class represented by the directed PAG 
  necessarily contains DMGs where a parental relation is present as well as DMGs where it is absent.

  \code{P} should be an adjacency matrix of type \link{amat.pag} that contains no undirected
  and circle-tail edges.
}
\references{
  J. M. Mooij and T. Claassen (2020).
  Constraint-Based Causal Discovery using Partial Ancestral Graphs in the presence of Cycles.
  \emph{In Proc. of the 36th Conference on Uncertainty in Artificial Intelligence (UAI-20)}, 1159-1168.
}
\author{
  Joris Mooij.
}
\examples{
##################################################
## Mooij et al. (2020), Fig. 43(a), p. 97
##################################################

# Encode ADMG as adjacency matrix
p <- 8 # total number of variables
V <- c("Ca","Cb","Cc","X0","X1","X2","X3","X4") # 3 context variables, 5 system variables
# amat[i,j] = 0 iff no edge btw i,j
# amat[i,j] = 1 iff i *-o j
# amat[i,j] = 2 iff i *-> j
# amat[i,j] = 3 iff i *-- j
amat <- rbind(c(0,2,2,2,0,0,0,0),
              c(2,0,2,0,2,0,0,0),
              c(2,2,0,0,2,2,0,0),
              c(3,0,0,0,0,0,2,0),
              c(0,3,3,0,0,3,0,2),
              c(0,0,3,0,2,0,0,0),
              c(0,0,0,3,0,0,0,2),
              c(0,0,0,0,2,0,3,0))
rownames(amat)<-V
colnames(amat)<-V

# Make use of d-separation oracle as "independence test"
indepTest <- dsepAMTest
suffStat<-list(g=amat,verbose=FALSE)

\donttest{
# Derive PAG that represents the Markov equivalence class of the ADMG with the FCI algorithm
# (assuming no selection bias)
fci.pag <- fci(suffStat,indepTest,alpha = 0.5,labels = V,verbose=TRUE,selectionBias=FALSE)

# Read off causal features from the FCI PAG
cat('Identified absence (-1) and presence (+1) of direct causal relations from FCI PAG:\n')
print(pag2edge(fci.pag@amat))
}
}
\keyword{graphs}
