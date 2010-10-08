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
  \item{amat}{adjacency matrix M with edge marks. The edge marks are
  coded in the following way: M[i,j]=M[j,i]=0: no edge; M[i,j]=1, M[j,i]
  != 0: i *-o j; M[i,j]=2, M[j,i]!=0: i*->j; M[i,j]=3, M[j,i]!=0: i*-j.}
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}

\seealso{\code{\link{fci}}}

\keyword{multivariate}
\keyword{models}
\keyword{graphs}


