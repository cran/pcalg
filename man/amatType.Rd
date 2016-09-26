\name{amatType}
\title{Types and Display of Adjacency Matrices in Package 'pcalg'}
\alias{amatType}
\alias{amat.cpdag}
\alias{amat.pag}
\alias{coerce,LINGAM,amat-method}%  as(<LINGAM>,  "amat")
\alias{coerce,pcAlgo,amat-method}%  as(<pcAlgo>,  "amat")
\alias{coerce,fciAlgo,amat-method}% as(<fciAlgo>, "amat")
%% the as(*, "matrix") methods where introduced even earlier.
\alias{coerce,pcAlgo,matrix-method}
\alias{coerce,fciAlgo,matrix-method}
%% hidden, but for us:
\alias{show.pc.amat}
\alias{show.fci.amat}
\description{
  Two types of adjacency matrices are used in package \pkg{pcalg}: Type
  \code{amat.cpdag} for DAGs and CPDAGs and type \code{amat.pag} for
  MAGs and PAGs.  The required type of adjacency matrix is documented
  in the help files of the respective functions or classes.  If in some
  functions more detailed information on the graph type is needed
  (i.e. DAG or CPDAG; MAG or PAG) this information will be passed in a
  separate argument (see e.g. \code{\link{gac}} and the examples below).

  Note that you get (\sQuote{extract}) such adjacency matrices as (S3)
  objects of \code{\link{class}} \code{"amat"} via the usual
  \code{\link{as}(., "<class>")} coercion, \preformatted{  as(from, "amat")
}%pre
}
\arguments{
  \item{from}{an \R object of

    class \code{\linkS4class{pcAlgo}}, as returned from
    \code{\link{skeleton}()} or \code{\link{pc}()} or an object of

    class \code{\linkS4class{fciAlgo}}, as from \code{\link{fci}()}
    (or \code{\link{rfci}}, \code{\link{fciPlus}}, and
    \code{\link{dag2pag}}), or an object of

    class \code{"LINGAM"} as returned from \code{\link{lingam}()}.}
}
\details{
  Adjacency matrices are integer valued square matrices with zeros on the
  diagonal.  They can have row- and columnnames; however, most functions
  will work on the (integer) column positions in the adjacency matrix.

  \bold{Coding for type \code{amat.cpdag}:}
  \describe{
    \item{\code{0}:}{No edge or tail}
    \item{\code{1}:}{Arrowhead}
  }
  Note that the edgemark-code refers to the \emph{row} index (as opposed
  adjacency matrices of type mag or pag).  E.g.:\preformatted{%
    amat[a,b] = 0  and  amat[b,a] = 1   implies a --> b.
    amat[a,b] = 1  and  amat[b,a] = 0   implies a <-- b.
    amat[a,b] = 0  and  amat[b,a] = 0   implies a     b.
    amat[a,b] = 1  and  amat[b,a] = 1   implies a --- b.}

  \bold{Coding for type \code{amat.pag}:}
  \describe{
    \item{\code{0}:}{No edge}
    \item{\code{1}:}{Circle}
    \item{\code{2}:}{Arrowhead}
    \item{\code{3}:}{Tail}
  }
  Note that the edgemark-code refers to the \emph{column} index (as
  opposed adjacency matrices of type dag or cpdag).  E.g.:\preformatted{%
  amat[a,b] = 2  and  amat[b,a] = 3   implies   a --> b.
  amat[a,b] = 3  and  amat[b,a] = 2   implies   a <-- b.
  amat[a,b] = 2  and  amat[b,a] = 2   implies   a <-> b.
  amat[a,b] = 1  and  amat[b,a] = 3   implies   a --o b.
  amat[a,b] = 0  and  amat[b,a] = 0   implies   a     b.}
}
\seealso{E.g. \code{\link{gac}} for a function which takes an
  adjacency matrix as input; \code{\linkS4class{fciAlgo}} for a class
  which has an adjacency matrix in one slot.

  \code{\link{getGraph}(x)} extracts the \code{\linkS4class{graph}}
  object from \code{x}, whereas \code{as(*, "amat")} gets the
  corresponding adjacency matrix.
}
\examples{
##################################################
## Function gac() takes an adjecency matrix of
## any kind as input.  In addition to that, the
## precise type of graph (DAG/CPDAG/MAG/PAG) needs
## to be passed as a different argument
##################################################
## Adjacency matrix of type 'amat.cpdag'
m1 <- matrix(c(0,1,0,1,0,0, 0,0,1,0,1,0, 0,0,0,0,0,1,
               0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0), 6,6)
## more detailed information on the graph type needed by gac()
gac(m1, x=1,y=3, z=NULL, type = "dag")

## Adjacency matrix of type 'amat.cpdag'
m2 <- matrix(c(0,1,1,0,0,0, 1,0,1,1,1,0, 0,0,0,0,0,1,
               0,1,1,0,1,1, 0,1,0,1,0,1, 0,0,0,0,0,0), 6,6)
## more detailed information on the graph type needed by gac()
gac(m2, x=3, y=6, z=c(2,4), type = "cpdag")

## Adjacency matrix of type 'amat.pag'
m3 <- matrix(c(0,2,0,0, 3,0,3,3, 0,2,0,3, 0,2,2,0), 4,4)
## more detailed information on the graph type needed by gac()
mg3 <- gac(m3, x=2, y=4, z=NULL, type = "mag")
pg3 <- gac(m3, x=2, y=4, z=NULL, type = "pag")
%% FIXME:  mg3 is identical to pg3! -- need an example where 'type' matters!

############################################################
## as(*, "amat") returns an adjacency matrix incl. its type
############################################################
## Load predefined data
data(gmG)
n <- nrow    (gmG8$x)
V <- colnames(gmG8$x)

## define sufficient statistics
suffStat <- list(C = cor(gmG8$x), n = n)
## estimate CPDAG
skel.fit <- skeleton(suffStat, indepTest = gaussCItest,
                     alpha = 0.01, labels = V)
## Extract the "amat" [and show nicely via  'print()' method]:
as(skel.fit, "amat")

##################################################
## Function fci() returns an adjacency matrix
## of type amat.pag as one slot.
##################################################
set.seed(42)
p <- 7
## generate and draw random DAG :
myDAG <- randomDAG(p, prob = 0.4)

## find skeleton and PAG using the FCI algorithm
suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
res <- fci(suffStat, indepTest=gaussCItest,
           alpha = 0.9999, p=p, doPdsep = FALSE)
str(res)
## get the a(djacency) mat(rix)  and nicely print() it:
as(res, "amat")

##################################################
## pcAlgo object
##################################################
## Load predefined data
data(gmG)
n <- nrow    (gmG8$x)
V <- colnames(gmG8$x)

## define sufficient statistics
suffStat <- list(C = cor(gmG8$x), n = n)
## estimate CPDAG
skel.fit <- skeleton(suffStat, indepTest = gaussCItest,
                     alpha = 0.01, labels = V)
## Extract Adjacency Matrix - and print (via method 'print.amat'):
as(skel.fit, "amat")

pc.fit <- pc(suffStat, indepTest = gaussCItest,
             alpha = 0.01, labels = V)
pc.fit # (using its own print() method 'print.pcAlgo')

as(pc.fit, "amat")
}
\keyword{graphs}
