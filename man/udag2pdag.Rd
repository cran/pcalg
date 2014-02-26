\name{udag2pdag}
\alias{udag2pdag}
\alias{udag2pdagRelaxed}
\alias{udag2pdagSpecial}
\title{Last PC Algorithm Step: Extend Object with Skeleton to Completed PDAG}
\description{
  These functions perform the last step in the PC algorithm:
  Transform an object of the class
  \code{"\linkS4class{pcAlgo}"} containing a skeleton and corresponding
  conditional independence information into a completed partially
  directed acyclic graph (CPDAG).  The functions first determine the
  v-structures, and then apply the three orientation rules as described
  in Sprirtes et al (2000) and Meek (1995) to orient as many of the
  remaining edges as possible.

  In the oracle version and when all assumptions hold, all three
  functions yield the same CPDAG.  In the sample version, however, the
  resulting CPDAG may be invalid in the sense that one cannot extend it
  a DAG without additional unshielded colliders by orienting the
  undirecting edges.  This can for example happen due to errors in the
  conditional indepedence tests or violations of the faithfulness
  assumption.  The three functions deal with such conflicts in different
  ways, as described in Details.
}
\usage{
udag2pdag       (gInput, verbose)
udag2pdagRelaxed(gInput, verbose, unfVect=NULL, solve.confl=FALSE)
udag2pdagSpecial(gInput, verbose, n.max=100)
}
\arguments{
  \item{gInput}{\code{"pcAlgo"}-object containing skeleton and
    conditional indepedence information.}
  \item{verbose}{0: No output; 1: Details}
  \item{unfVect}{vector containing numbers that encode ambiguous
    triples (as returned by \code{\link{pc.cons.intern}}).  This is
    needed in the conservative and majority rule PC algorithms.}
  \item{solve.confl}{if \code{TRUE}, the orientation of the v-structures and the
    orientation rules work with lists for candidate sets and allow
    bi-directed edges to resolve conflicting edge orientations. Note that
    therefore the resulting object is order-independent but might not be
    a PDAG because bi-directed edges can be present.}
  \item{n.max}{maximum number of tries for re-orienting doubly visited
    edges in \code{udag2pdagSpecial}.}
}
\details{
  \describe{
    \item{for \code{\link{udag2pdag}}:}{%
      If there are edges that are part of more than one v-structure
      (i.e., the edge b - c in the v-structures a -> b <- c and b -> c
      <- d), earlier edge orientations are simply overwritten by later
      ones.  Thus, if a -> b <- c is considered first, the edge b - c is
      first oriented as b <- c and later overwritten by b -> c.  The
      v-structures are considered in lexicographical ordering.

      If the resulting graph is extendable to a DAG without additional
      v-structures, then the rules of Meek (1995) and Spirtes et al
      (2000) are applied to obtain the corresponding CPDAG.  Otherwise,
      the edges are oriented randomly to obtain a DAG that fits on the
      skeleton, discarding all information about the v-structures.  The
      resulting DAG is then transformed into its CPDAG.  Note that the
      output of \code{udag2pdag} is random whenever the initial graph
      was not extendable. %%TODO: please check if this is correct

      Although the output of \code{\link{udag2pdag}} is always
      extendable, it is not necessarily a valid CPDAG in the sense that
      it describes a Markov equivalence class of DAGs.  For example, two
      v-structures a -> b <- c and b -> c <- d (considered in this
      order) would yield the output a -> b -> c <- d.  This is
      extendable to a DAG (it already \emph{is} a DAG), but it does not
      describe a Markov equivalence class of DAGs, since the DAG a <- b
      -> c <- d describes the same conditional independencies.
    }

    \item{for \code{\link{udag2pdagSpecial}}:}{%

      If the graph after orienting the v-structures as in
      \code{udag2pdag} is extendable to a DAG without additional
      v-structures, then the rules of Meek (1995) and Spirtes et al
      (2000) are applied to obtain the corresponding CPDAG.  Otherwise,
      the algorithm tries at most \code{n.max} different random
      orderings of the v-structures (hence overwriting orientations in
      different orders), until it finds one that yields an extendable
      CPDAG.  If this fails, the edges are oriented randomly to obtain a
      DAG that fits on the skeleton, discarding all information about
      the v-structures.  The resulting DAG is then transformed into its
      CPDAG.  Note that the output of \code{udag2pdagSpecial} is random
      whenever the initial graph was not extendable.

      Although the output of \code{\link{udag2pdag}} is always
      extendable, it is not necessarily a valid CPDAG in the sense that
      it describes a Markov equivalence class of DAGs.  For example, two
      v-structures a -> b <- c and b -> c <- d (considered in this
      order) would yield the output a -> b -> c <- d.  This is
      extendable to a DAG (it already IS a DAG), but it does not
      describe a Markov equivalence class of DAGs, since the DAG a <- b
      -> c <- d describes the same conditional independencies.
    }

    \item{for \code{\link{udag2pdagRelaxed}}:}{%

      This is the default version in the PC/RFCI/FCI algorithm.  It does
      \bold{not} test whether the output is extendable to a DAG without
      additional v-structures.

      If \code{unfVect = NULL} (no ambiguous triples), the three
      orientation rules are applied to each eligible structure until no
      more edges can be oriented.  Otherwise, \code{unfVect} contains
      the numbers of all ambiguous triples in the graph as determined by
      \code{\link{pc.cons.intern}}.  Then the orientation rules take
      this information into account.  For example, if a -> b - c and
      <a,b,c> is an unambigous triple and a non-v-structure, then rule 1
      implies b -> c.  On the other hand, if a -> b - c but <a,b,c> is
      an ambiguous triple, then the edge b - c is not oriented.

      If \code{solve.confl = FALSE}, earlier edge orientations are
      overwritten by later ones as in \code{udag2pdag} and
      \code{udag2pdagSpecial}.

      If \code{solv.confl = TRUE}, both the v-structures and the
      orientation rules work with lists for the candidate edges and
      allow bi-directed edges if there are conflicting orientations. For
      example, two v-structures a -> b <- c and b -> c <- d then yield a
      -> b <-> c <-d.  This option can be used to get an
      order-independent version of the PC algorithm (see Colombo and
      Maathuis (2013)).

      We denote bi-directed edges, for example between two variables i
      and j, in the adjacency matrix M of the graph as M[i,j]=2 and
      M[j,i]=2.  Such edges should be interpreted as indications of
      conflicts in the algorithm, for example due to errors in the
      conditional independence tests or violations of the faithfulness
      assumption.
    }
  }
}
\value{
  \describe{
    \item{for \code{udag2pdag()} and \code{udag2pdagRelaxed()}:}{%
      oriented \code{"pcAlgo"}-object.}

    \item{for \code{udag2pdagSpecial}:}{a \code{\link{list}} with
      components \describe{
	\item{pcObj}{An oriented \code{"pcAlgo"}-object.}
	\item{evisit}{Matrix counting the number of orientation attempts per edge}
	%% TODO both xtbl.orig and xtbl:  what is different from 'status'  <<<<<,
	\item{xtbl.orig}{Logical indicating whether the original graph
	  with v-structures is extendable.}
	\item{xtbl}{Logical indicating whether the final graph with
	  v-structures is extendable}
	\item{amat0}{Adjacency matrix of original graph with
	  v-structures (coding 0/1 for no edge or tail / arrowhead;
	  e.g. \code{amat[a,b] = 0} and \code{amat[b,a] = 1} implies a -> b) .}
	\item{amat1}{Adjacency matrix of final graph with v-structures
	  after changing the ordering in which the v-structures are
	  considered (coding 0/1 for no edge or tail / arrowhead;
	  e.g. \code{amat[a,b] = 0} and \code{amat[b,a] = 1} implies a -> b) .}
	\item{status}{Integer code with values \describe{
	    \item{0:}{Original try is extendable;}
	    \item{1:}{Reorienting double edge visits helps;}
	    \item{2:}{Original try is not extendable; reorienting double
	      visits does not help; result is acyclic, has original
	      v-structures, but perhaps additional v-structures.}
	  }
	}
	\item{counter}{Number of orderings of the v-structures until
	  success or \code{n.max}.}
      }}
  }
}
\references{
  C. Meek (1995). Causal inference and causal explanation with
  background knowledge. In \emph{Proceedings of the Eleventh Conference
    on Uncertainty in Artificial Intelligence (UAI-95)},
  pp. 403-411. Morgan Kaufmann Publishers, Inc.

  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.

  J. Pearl (2000), \emph{Causality}, Cambridge University Press.

  D. Colombo and M.H. Maathuis (2013).
  Order-independent constraint-based causal structure learning.
  arXiv:1211.3295v2.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{pc}}, \code{\link{pdag2dag}},
  \code{\link{dag2cpdag}}, \code{\link{udag2pag}},
  \code{\link{udag2apag}}, \code{\link{dag2pag}}.
}
\examples{
## simulate data
set.seed(123)
p <- 10
myDAG <- randomDAG(p, prob = 0.2)
trueCPDAG <- dag2cpdag(myDAG)
n <- 1000
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

## estimate skeleton
resU <- skeleton(suffStat = list(C = cor(d.mat), n = n),
                 indepTest = gaussCItest, ## (partial correlations)
                 alpha = 0.05, p=p)

## orient edges using three different methods
resD1 <- udag2pdagRelaxed(resU, verbose=0)
resD2 <- udag2pdagSpecial(resU, verbose=0, n.max=100)
resD3 <- udag2pdag       (resU, verbose=0)

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
