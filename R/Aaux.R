`%nin%` <- function (x, table) is.na(match(x, table))

## TODO 1) Fast version of  !all(x %in% tab) == any(x %nin% tab)
##      2) Fast version of  !any(x %in% tab) == all(x %nin% tab)
##     3a) Fast version of  !(x %in% A  |  x %in  B)
##     3b) Fast version of  !(x %in% A ||  x %in  B) = 3a) where length(x)==1

##' Compute \eqn{\log( (1+r)/(1-r) )} in a numerically stable way.  Note that
##' \eqn{(1+r)/(1-r) = 1 + 2r/(1-r)};  hence the logarithm is
##' \eqn{log(*) = log(1 + 2r/(1-r)) = log1p(2r/(1-r))}.
##'
##' @title  log( (1+r)/(1-r) )  in a Numerically Stable Way
##' @param r numeric (vector) in \eqn{[0,1]}.
##' @return numeric vector \dQuote{as} \code{r}.
##' @author Martin Maechler
logQ1pm <- function(r) log1p(2*r/(1-r))


## MM: "Lifted" from Matrix package ~/R/Pkgs/Matrix/R/Auxiliaries.R
## "Theory" behind this: /u/maechler/R/MM/MISC/lower-tri-w.o-matrix.R
indTri <- function(n, upper = TRUE, diag = FALSE) {
    ## Indices of (strict) upper/lower triangular part
    ## == which(upper.tri(diag(n), diag=diag) or
    ##	  which(lower.tri(diag(n), diag=diag) -- but
    ## more efficiently for largish 'n'
    stopifnot(length(n) == 1, n == (n. <- as.integer(n)), (n <- n.) >= 0)
    if(n <= 2) {
        if(n == 0) return(integer(0))
        if(n == 1) return(if(diag) 1L else integer(0))
        ## else n == 2
        v <- if(upper) 3L else 2L
	return(if(diag) c(1L, v, 4L) else v)
    }

    ## n >= 3 [also for n == 2 && diag (==TRUE)] :

    ## First, compute the 'diff(.)' of the result [fast, using integers]
    n. <- if(diag) n else n - 1L
    n1 <- n. - 1L
    ## all '1' but a few
    r <- rep.int(1L, choose(n.+1, 2) - 1)
    tt <- if(diag) 2L else 3L
    r[cumsum(if(upper) 1:n1 else n.:2)] <- if(upper) n:tt else tt:n
    ## now have differences; revert to "original":
    cumsum(c(if(diag) 1L else if(upper) n+1L else 2L, r))
}

##' Count Edges, for directed graphs, *not* as numEdges(g)
##' but such that  "<-->" or "<--o" ... counts as 1
##' (for undirected graphs, this is the same as numEdges(.))
numGedges <- function(amat) {
    dimnames(amat) <- NULL # speed
    A <- (amat + t(amat)) != 0
    n <- nrow(A)
    if(n <= 40) ## faster
	sum(A[lower.tri(A)])
    else
	sum(A[indTri(n)])
}

check.Rgraphviz <- function() {
  if(!requireNamespace("Rgraphviz"))
    stop("Package 'Rgraphviz' (from Bioconductor) must be installed for plotting graphs!")
}

