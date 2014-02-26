\name{NEWS}
\title{pcalg News}
\encoding{UTF-8}
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


