\name{pcalg-internal}
%
\alias{adj.check}
\alias{allDags}
\alias{amat2dag}
\alias{check.new.coll}
\alias{checkTriple}
\alias{ci.test}
\alias{discr.path}
\alias{faith.check}
\alias{find.sink}
\alias{find.upd}
\alias{flipEdges}
\alias{gAlgo-class}%-> *do* document! (FIXME)
\alias{labelEdges}
\alias{legal.psep}
\alias{lm.cov}
\alias{make.edge.df}
\alias{orderEdges}
\alias{pcAlgo.Perfect}
\alias{reach}
\alias{triple2numb}% maybe document ?
\alias{trueCov}
\alias{ucp}
\alias{upd}
\alias{wgtMatrix}
%
\title{Internal Pcalg Functions}
\description{
  Internal functions of package \pkg{pcalg}
}
\usage{
wgtMatrix(g, transpose=TRUE)
adj.check(gm,x)
amat2dag(amat)
find.sink(gm)
labelEdges(amat)
make.edge.df(amat)
orderEdges(amat)
}
\details{
  These are not to be called by the user.
}
\keyword{ internal }
