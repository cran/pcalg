\name{isValidGraph}
\alias{isValidGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check for a DAG, CPDAG or a maximally oriented PDAG
  %%  ~~function to do ... ~~
}
\description{
Check whether the adjacency matrix \code{amat} matches the specified \code{type}.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
isValidGraph(amat, type = c("pdag", "cpdag", "dag"), verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{
    adjacency matrix of type \code{amat.cpdag} (see \code{\link{amatType}})
}
  \item{type}{
    string specifying the type of graph of the adjacency matrix amat. It
    can be a DAG (\code{type="dag"}), a CPDAG (\code{type="cpdag"}) or a
    maximally oriented PDAG (\code{type="pdag"}) from Meek (1995). 
}
  \item{verbose}{
    If TRUE, detailed output on why the graph might not be valid is provided.
}
}
\details{
For a given adjacency matrix \code{amat} and graph \code{type}, this
function checks whether the two match. 

For \code{type = "dag"} we require that \code{amat} does NOT contain
directed cycles. 

For \code{type = "cpdag"} we require that \code{amat} does NOT contain
directed or partially directed cycles. We also require that the
undirected part of the CPDAG (represented by \code{amat}) is made up of
chordal components and that our graph is maximally oriented according to
rules from Meek (1995). 

For \code{type = "pdag"} we require that \code{amat} does NOT contain
directed cycles. We also require that the PDAG is maximally oriented
according to rules from Meek (1995). Additionally, we require that the
adjacency matrix \code{amat1} of the CPDAG corresponding to our PDAG
(represented by \code{amat}), satisfies \code{isValidGraph(amat =
  amat1,type = "cpdag") == TRUE} and that there is no mismatch in the
orientations implied by \code{amat} and \code{amat1}. We obtain
\code{amat1} by extracting the skeleton and v-structures from
\code{amat} and then closing the orientation rules from Meek (1995). 
}
\value{
TRUE, if the adjacency matrix \code{amat} is of the \code{type} specified and FALSE, otherwise.
}
\references{
 C. Meek (1995). Causal inference and causal explanation with background knowledge, In Proceedings of UAI 1995, 403-410.
}
\author{
  Emilija Perkovic and Markus Kalisch
}
\examples{
## a -> b -> c
amat <- matrix(c(0,1,0, 0,0,1, 0,0,0), 3,3)
colnames(amat) <- rownames(amat) <- letters[1:3]
## graph::plot(as(t(amat), "graphNEL"))             
isValidGraph(amat = amat, type = "dag") ## is a valid DAG
isValidGraph(amat = amat, type = "cpdag") ## not a valid CPDAG 
isValidGraph(amat = amat, type = "pdag") ## is a valid PDAG

## a -- b -- c
amat <- matrix(c(0,1,0, 1,0,1, 0,1,0), 3,3)
colnames(amat) <- rownames(amat) <- letters[1:3]
## plot(as(t(amat), "graphNEL"))             
isValidGraph(amat = amat, type = "dag") ## not a valid DAG
isValidGraph(amat = amat, type = "cpdag") ## is a valid CPDAG
isValidGraph(amat = amat, type = "pdag") ## is a valid PDAG

## a -- b -- c -- d -- a
amat <- matrix(c(0,1,0,1, 1,0,1,0, 0,1,0,1, 1,0,1,0), 4,4)
colnames(amat) <- rownames(amat) <- letters[1:4]
## plot(as(t(amat), "graphNEL"))             
isValidGraph(amat = amat, type = "dag") ## not a valid DAG
isValidGraph(amat = amat, type = "cpdag") ## not a valid CPDAG
isValidGraph(amat = amat, type = "pdag") ## not a valid PDAG
}
\keyword{graphs}% use one of  RShowDoc("KEYWORDS")

