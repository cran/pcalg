\name{NEWS}
\title{pcalg News}
\encoding{UTF-8}
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


