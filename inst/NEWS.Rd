% Check from R:
% news(db= tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/pcalg/inst/NEWS.Rd"))
\name{NEWS}
\title{pcalg News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 2.3-0 (2015-..-.., svn ...)}{
  \itemize{
    \item \code{print()} and \code{summary()} methods for \code{pc()},
    \code{fci()}, etc.

    \item new \code{as(*, "amat")} methods to get the adjacency matrix
    of a model, including a \code{print()} method for "amat"s.

    \code{?amatType} documents the different kinds of adjacency matrices
    in \pkg{pcalg}.

    \item New functions ...... % FIXME -- only if *not* mentioned below

    \item Deprecated \code{LINGAM()} in favor of new \code{lingam()}.
    \item Deprecated \code{allDags()} in favor of \code{pdag2allDags()}.
    %% more deprecated ?
  }
}

\section{CHANGES IN VERSION 2.2-4 (2015-07-22, svn r344)}{
  \itemize{
    \item Bug fix in \code{simulate(<gies>)}.
    \item Bug fix in \code{LINGAM(., verbose=FALSE)} (the default, now tested).
    \item Fix \code{plot(pc(..), zvalue.lwd=TRUE)}, and better document
    that plot method.
    \item Much faster \code{gSquareDis()} and \code{gSquareBin()} for
     \eqn{|S|} in 1:4 or so.
    \item More test coverage, notably for \code{pcSelect()} and \code{mcor()}.
    \item \code{r.gauss.pardag()} gets an optional argument \code{labels}.
  }
}
\section{CHANGES IN VERSION 2.2-3 (2015-07-12, svn r322)}{
  \itemize{
    \item Adapt to API change in \pkg{igraph} package.
    \item More test coverage and first speed tweaks.
  }
}
\section{CHANGES IN VERSION 2.2-2 (2015-06-24, svn r304)}{
  \itemize{
    \item Function \code{gac} for testing if a set of nodes satisfies
    the Generalized Adjustment Criterion (GAC). GAC is a sufficient and
    necessary criterion for estimating causal effects via covariate
    adjustment in DAGs, CPDAGs, MAGs and PAGs.
  }
}

\section{CHANGES IN VERSION 2.2-0}{
  \itemize{
    \item Function \code{randDAG} for generating random directed acyclic
    graphs (DAGs) with fixed expected number of neighbours. Several
    different methods are provided, each intentionally biased towards
    certain properties. The methods are based on the analogue functions
    in the \pkg{igraph} package.
    \item Function \code{unifDAG} and \code{unifDAG.approx} for uniform
    sampling of a labelled directed acyclic graph (DAG) with combinatorial
    enumeration.
  }
}

\section{CHANGES IN VERSION 2.1-0}{
  \itemize{
    \item Function \code{LINGAM} for estimating the structure of Linear
    non-Gaussian Additive Models
  }
}

\section{CHANGES IN VERSION 2.0-4}{
  \itemize{
    \item Function \code{fciPlus} for a fast (but accurate) version of
    \code{fci} (note that the already existing function \code{rfci} is a
    fast \emph{approximation} of \code{fci})
  }
}

\section{CHANGES IN VERSION 2.0-1}{
  \itemize{
    \item Functions \code{gds}, \code{gies}, \code{ges} for score based learning
    \item Function \code{backdoor} for generalized backdoor criterion
    \item Bugfix in function \code{wgtMatrix}; this affects functions
    \code{trueCov} and
    \code{rmvDAG}; to reproduce old version, use option 'back.compatible =
    TRUE'
    \item \code{skeleton}, \code{pc}, \code{fci} and \code{rfci}  now
    have 'alpha' as third argument,
    then 'labels, p' where 'p' is optional; namely, derived from 'labels'
    when these are specified.
    \item The default for extending the skeleton to a CPDAG is now
    \code{udag2pdagRelaxed}, which is completely deterministic (was
    \code{udag2pdag}
    before, which can use random orientations if sampling errors occur).
  }
}


