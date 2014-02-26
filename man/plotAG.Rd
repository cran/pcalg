\name{plotAG}
\alias{plotAG}
\title{Plot partial ancestral graphs (PAG)}
\description{
This function is DEPRECATED! Use the \code{plot} method of the
\code{\linkS4class{fciAlgo}} class instead.
}
\usage{
plotAG(amat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{Adjacency matrix (coding 0,1,2,3 for no edge, circle,
    arrowhead, tail; e.g., \code{amat[a,b] = 2} and \code{amat[b,a] = 3}
    implies a -> b)} 
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}

\seealso{\code{\link{fci}}}

\keyword{multivariate}
\keyword{models}
\keyword{graphs}


