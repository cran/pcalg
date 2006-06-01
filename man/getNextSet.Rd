\name{getNextSet}
\alias{getNextSet}
\title{Iteration through a list of all combinations of choose(n,k)}
\description{
  Given a combination of \eqn{k} elements out of the elements
  \eqn{1,\dots,n}, the next 
  set of size \code{k} in a specified sequence is computed.
}
\usage{
getNextSet(n,k,set)
}
\arguments{
  \item{n}{number of elements to choose from (integer)}
  \item{k}{size of chosen set (integer)}
  \item{set}{previous set in list (numeric vector)}
}
\value{List with two elements:
  \item{nextSet}{next set in list (numeric vector)}
  \item{wasLast}{logical indicating whether the end of the specified sequence is reached.}
}
\details{
  The initial set is \code{1:k}.  Last index varies quickest.  Using the
  dynamic creation of sets reduces the memory demands dramatically for
  large sets.  If complete lists of combination sets have to be produced
  and memory is no problem, the function \code{\link[combinat]{combn}}
  from package \pkg{combinat} is an alternative.
}
\author{
  Markus Kalisch \email{kalisch@stat.math.ethz.ch} and Martin Maechler
}
\examples{
## start from first set (1,2) and get the next set of size 2 out of 1:5
## notice that res$wasLast is FALSE :
str(r <- getNextSet(5,2,c(1,2)))

## input is the last set; notice that res$wasLast now is TRUE:
str(r2 <- getNextSet(5,2,c(4,5)))

## Show all sets of size k out of 1:n :
## {if you really want this in practice, use something like combn() !}
n <- 5
k <- 3
currentSet <- 1:k
(res <- rbind(currentSet, deparse.level = 0))
repeat {
  newEl <- getNextSet(n,k,currentSet)
  if (newEl$wasLast)
     break
  ## otherwise continue:
  currentSet <- newEl$nextSet
  res <- rbind(res, currentSet, deparse.level = 0)
}
res
stopifnot(choose(n,k) == nrow(res)) ## must be identical
}
\keyword{arith}
\keyword{utilities}
