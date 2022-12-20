% Check from R:
% news(db= tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/pcalg/inst/NEWS.Rd"))
\name{NEWS}
\title{pcalg News}
\encoding{UTF-8}


\section{CHANGES IN VERSION 2.7-8 (2022-12-17, svn r528)}{
  \itemize{
    \item Fix \code{idaFast()} regression tests, using \code{all.equal()}
    instead of \code{round(*) == round(*)}, notably for ATLAS-using
    BLAS/Lapack.
  }
}

\section{CHANGES IN VERSION 2.7-7 (2022-10-04, svn r526)}{
  \itemize{
    \item Tiny bug fix for fci.
  }
}

\section{CHANGES IN VERSION 2.7-6 (2022-04-04, svn r525)}{
  \itemize{
    \item \sQuote{Proper} registration of \code{.Call()}ed C routines.
    \item Not using \code{if(class(.) = "...")}
    \item Adapt to upcoming deprecation of C++ \code{str::binary_function()}
    \item \pkg{dagitty} masking \code{randomDAG()}.
  }
}

\section{CHANGES IN VERSION 2.7-5 (2022-02-21, svn r518)}{
  \itemize{
    \item fixing bug. % no tests, no doc (hmm)
  }
}

\section{CHANGES IN VERSION 2.7-4 (2021-11-24, svn r517)}{
  \itemize{
    \item updated links, removed old vignette.
  }
}

\section{CHANGES IN VERSION 2.7-3 (2021-05-20, svn r516)}{
  \itemize{
    \item adapt to R-devel's new/changed warning of
    \code{matrix(<non-matching-length>, n,m)}.
  }
}

%% Martin's shell script (run at SfS where there's a CRAN mirror):

%% R-pkg-archive-show pcalg | tac    {plus manually sorting 2.6-12 etc}

%% 2.7-6 &          2022-4-4 &           2022-04-04 11:50:02 UTC
%% 2.7-5 &          2022-2-21 &          2022-02-22 12:40:02 UTC
%% 2.7-4 &          2021-11-24 &         2021-11-25 10:10:02 UTC
%% 2.7-3 &          2021-5-20 &          2021-06-01 08:30:02 UTC
%% 2.7-2 &          2021-4-20 &          2021-04-22 04:50:13 UTC
%% 2.7-1 &          2021-1-8 &           2021-01-09 06:30:08 UTC
%% 2.7-0 &          2020-11-15 &         2020-11-17 08:20:10 UTC
%% 2.6-12 &         2020-09-22 &         2020-09-22 15:30:03 UTC
%% 2.6-11 &         2020-08-12 &         2020-08-16 16:00:02 UTC
%% 2.6-10 &         2020-03-3 &          2020-03-03 14:40:07 UTC
%% 2.6-9 &          2020-02-11 &         2020-02-11 13:00:06 UTC
%% 2.6-8 &          2019-11-18 &         2019-11-21 12:00:06 UTC
%% 2.6-7 &          2019-10-20 &         2019-10-23 14:40:02 UTC
%% 2.6-6 &          2019-09-16 &         2019-09-16 15:50:05 UTC
%% 2.6-5 &          2019-08-26 &         2019-08-27 12:00:05 UTC
%% 2.6-2 &          2019-04-08 &         2019-04-25 22:40:03 UTC
%% 2.6-0 &          2018-05-14 &         2018-06-04 18:23:47 UTC
%% 2.5-0 &          2017-07-11 &         2017-07-12 10:41:16 UTC
%% 2.4-5 &          2017-02-17 &         2017-02-22 09:57:28
%% 2.4-3 &          2016-09-21 &         2016-09-28 17:10:41

%% 2.4-2 &          2016-09-21 &         2016-09-26 16:47:02
\section{CHANGES IN VERSION 2.4-2 (2016-09-22, svn r403)}{

%% FIXME: The comments below are *ONLY* up to 2015-09-09, r365 :

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

%%==> ../ChangeLog <<<  for older versions
%%    ~~~~~~~~~~~~

%% Version          Date  Date CRAN Publication
%%   1.1-6    2013-03-22    2013-03-22 15:46:01
%%   1.1-5    2012-04-26    2012-04-26 17:34:25
%%   1.1-4    2010-05-06    2011-05-09 15:30:33
%%   1.1-2    2010-10-08    2010-10-27 05:09:28
%%   1.0-2    2010-03-26    2010-04-16 09:51:38
%%   1.0-0    2010-03-09    2010-03-14 07:50:36
%%   0.1-9    2009-3-03     2009-09-23 16:45:42
%%   0.1-8    2009-3-03
%%   0.1-7    2008-4-16
%%   0.1-6    2008-3-31
%%   0.1-5    2007-3-23
%%   0.1-4    2006-11-1
%%   0.1-3    2006-06-1

