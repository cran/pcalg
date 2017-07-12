\name{pcalg2dagitty}
\alias{pcalg2dagitty}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Transform the adjacency matrix from \pkg{pcalg} into a \pkg{dagitty} object
}
\description{
Transform the adjacency matrix of type \code{amat.cpdag} or
\code{amat.pag} (for details on coding see \code{\link{amatType}}).
}
\usage{pcalg2dagitty(amat, labels, type = "cpdag")}
\arguments{
  \item{amat}{adjacency matrix of type \code{amat.cpdag} or \code{amat.pag}}

  \item{labels}{\code{character} vector of variable (or “node”) names.}
  \item{type}{string specifying the type of graph of the adjacency matrix \code{amat}. It can be a DAG (\code{type="dag"}), a CPDAG (\code{type="cpdag"}) or a maximally oriented PDAG (\code{type="pdag"}) from Meek (1995); then the type of adjacency matrix is assumed to be \code{amat.cpdag}. It can also be a MAG (\code{type = "mag"}) or a PAG (\code{type="pag"}); then the type of the adjacency matrix is assumed to be \code{amat.pag}.}
}
\details{
For a given adjacency matrix \code{amat} the form \code{amat.cpdag} or \code{amat.pag} and a specified graph \code{type}, this function returns a dagitty object corresponding to the graph structure specified by \code{amat}, \code{labels} and \code{type}.
The resulting object is compatible with the \pkg{dagitty} package.
}
\value{
A dagitty graph (see the \pkg{dagitty} package).
}
\author{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Emilija Perkovic and Markus Kalisch
}
\examples{
data(gmG)
n <- nrow    (gmG8$x)
V <- colnames(gmG8$x) # labels aka node names

amat <- wgtMatrix(gmG8$g)
amat[amat != 0] <- 1
dagitty_dag1 <- pcalg2dagitty(amat,V,type="dag")
}