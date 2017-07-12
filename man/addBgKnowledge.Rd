\name{addBgKnowledge}
\alias{addBgKnowledge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Add background knowledge to a CPDAG or PDAG
}
\description{
Add background knowledge x -> y to an adjacency matrix and complete the orientation rules from Meek (1995).
}
\usage{
addBgKnowledge(gInput, x = c(), y = c(), verbose = FALSE, checkInput = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gInput}{\code{graphNEL} object or adjacency matrix of type \code{amat.cpdag} (see \code{\link{amatType}})}

  \item{x,y}{node labels of \code{x} or \code{y} in the adjacency matrix. \code{x} and \code{y} can be vectors representing several nodes (see details below).}
  \item{verbose}{If TRUE, detailed output is provided.}
  \item{checkInput}{If TRUE, the input adjacency matrix is carefully
  checked to see if it is a valid graph using function \code{\link{isValidGraph}}}
}

\details{
If the input is a \code{graphNEL} object, it will be converted into an adjacency matrix of type \code{amat.cpdag}.
If \code{x} and \code{y} are given and if \code{amat[y,x] != 0}, this function adds orientation x -> y to the adjacency matrix \code{amat} and completes the orientation rules from Meek (1995).
If \code{x} and \code{y} are not specified (or empty vectors) this function simply completes the orientation rules from Meek (1995). If \code{x} and \code{y} are vectors of length k, k>1, this function tries to add \code{x[i] -> y[i]} to the adjacency matrix amat and complete the orientation rules from Meek (1995) for every i in {1,...,k} (see Algorithm 1 in Perkovic et. al, 2017). 
}
\value{
An adjacency matrix of type \code{amat.cpdag} of the maximally oriented pdag with added background knowledge \code{x -> y} or \code{NULL}, if the backgound knowledge is not consistent with any DAG represented by the PDAG with the adjacency matrix \code{amat}.
}
\references{
C. Meek (1995). Causal inference and causal explanation with background knowledge, In Proceedings of UAI 1995, 403-410.
  
E. Perkovic, M. Kalisch and M.H. Maathuis (2017). Interpreting and using CPDAGs with background knowledge. In Proceedings of UAI 2017. 
}
\author{
Emilija Perkovic and Markus Kalisch
}
\examples{
## a -- b -- c
amat <- matrix(c(0,1,0, 1,0,1, 0,1,0), 3,3)
colnames(amat) <- rownames(amat) <- letters[1:3]
## plot(as(t(amat), "graphNEL"))             
addBgKnowledge(gInput = amat) ## amat is a valid CPDAG
## b -> c is directed; a -- b is not directed by applying
## Meek's orientation rules
bg1 <- addBgKnowledge(gInput = amat, x = "b", y = "c") 
## plot(as(t(bg1), "graphNEL"))
## b -> c and b -> a are directed
bg2 <- addBgKnowledge(gInput = amat, x = c("b","b"), y = c("c","a")) 
## plot(as(t(bg2), "graphNEL")) 

## c -> b is directed; as a consequence of Meek's orientation rules,
## b -> a is directed as well
bg3 <- addBgKnowledge(gInput = amat, x = "c", y = "b") 
## plot(as(t(bg3), "graphNEL")) 

amat2 <- matrix(c(0,1,0, 1,0,1, 0,1,0), 3,3)
colnames(amat2) <- rownames(amat2) <- letters[1:3]
## new collider is inconsistent with original CPDAG; thus, NULL is returned
addBgKnowledge(gInput = amat2, x = c("c", "a"), y = c("b", "b"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")

