\name{pc}
\alias{pc}
\title{Estimate the Equivalence Class of a DAG using the PC Algorithm}
\description{
  Estimate the equivalence class of a directed acyclic
  graph (DAG) from observational data, using the PC-algorithm.
}
\usage{
pc(suffStat, indepTest, alpha, labels, p,
   fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,
   u2pd = c("relaxed", "rand", "retry"),
   skel.method = c("stable", "original", "stable.fast"),
   conservative = FALSE, maj.rule = FALSE, solve.confl = FALSE, 
   numCores = 1, verbose = FALSE)
}
\arguments{
  \item{suffStat}{A \code{\link{list}} of sufficient statistics,
    containing all necessary elements for the conditional independence
    decisions in the function \code{indepTest}.}
  \item{indepTest}{A \code{\link{function}} for testing conditional
    independence.  It is internally called as
    \code{indepTest(x,y,S,suffStat)}, and tests conditional independence
    of \code{x} and \code{y} given \code{S}.  Here, \code{x} and
    \code{y} are variables, and \code{S} is a (possibly empty) vector of
    variables (all variables are denoted by their (integer) column positions
    in the adjacency matrix).  \code{suffStat} is a list, see the
    argument above.  The return value of \code{indepTest} is the p-value
    of the test for conditional independence.}
  \item{alpha}{significance level (number in \eqn{(0,1)} for the
    individual conditional independence tests.}
  \item{labels}{(optional) character vector of variable (or
    \dQuote{node}) names.  Typically preferred to specifying \code{p}.}
  \item{p}{(optional) number of variables (or nodes).  May be specified
    if \code{labels} are not, in which case \code{labels} is set to
    \code{1:p}.}
  \item{numCores}{Specifies the number of cores to be used for parallel
    estimation of \code{\link{skeleton}}.}
  \item{verbose}{If \code{TRUE}, detailed output is provided.}
  \item{fixedGaps}{A logical matrix of dimension p*p. If entry
    \code{[i,j]} or \code{[j,i]} (or both) are TRUE, the edge i-j is
    removed before starting the algorithm. Therefore, this edge is
    guaranteed to be absent in the resulting graph.}
  \item{fixedEdges}{A logical matrix of dimension p*p. If entry
    \code{[i,j]} or \code{[j,i]} (or both) are TRUE, the edge i-j is
    never considered for removal. Therefore, this edge is
    guaranteed to be present in the resulting graph.}
  \item{NAdelete}{If indepTest returns \code{NA} and this option is
    \code{TRUE}, the corresponding edge is deleted. If this option is
    \code{FALSE}, the edge is not deleted.}
  \item{m.max}{Maximal size of the conditioning sets that are considered in the
    conditional independence tests.}
  \item{u2pd}{String specifying the method for dealing with conflicting
    information when trying to orient edges (see details below).}
  \item{skel.method}{Character string specifying method; the default,
    \code{"stable"} provides an \emph{order-independent} skeleton, see
    \code{\link{skeleton}}.}
  \item{conservative}{Logical indicating if the conservative PC is used.
    In this case, only option \code{u2pd = "relaxed"} is supported.
    Note that therefore the resulting object might not be extendable to
    a DAG.  See details for more information.}
  \item{maj.rule}{Logical indicating that the triples shall be checked
    for ambiguity using a majority rule idea, which is less strict than the
    conservative PC algorithm.  For more information, see details.}
  \item{solve.confl}{If \code{TRUE}, the orientation of the v-structures and the
    orientation rules work with lists for candidate sets and allow
    bi-directed edges to resolve conflicting edge orientations. In this
    case, only option \code{u2pd = relaxed} is supported. Note, that
    therefore the resulting object might not be a CPDAG because
    bi-directed edges might be present. See details
    for more information.}
}
\value{An object of \code{\link{class}} \code{"pcAlgo"} (see
  \code{\linkS4class{pcAlgo}}) containing an estimate of the equivalence
  class of the underlying DAG.
}
\details{
  Under the assumption that the distribution of the observed variables
  is faithful to a DAG, this function estimates the Markov equivalence class of
  the DAG. We do not estimate the DAG itself, because this is typically
  impossible (even with an infinite amount of data), since different
  DAGs can describe the same conditional independence relationships.
  Since all DAGs in an equivalence class describe the same conditional
  independence relationships, they are equally valid ways to
  describe the conditional dependence structure that was given as
  input.

  All DAGs in a Markov equivalence class have the same skeleton (i.e.,
  the same adjacency information) and the same v-structures (see
  definition below). However, the direction of some edges may be
  undetermined, in the sense that they point one way in one DAG in the
  equivalence class, while they point the other way in another DAG in
  the equivalence class.

  A Markov equivalence class can be uniquely represented by a completed
  partially directed acyclic graph (CPDAG). A CPDAG
  contains undirected and directed edges. The edges have the following
  interpretation: (i) there is a (directed or undirected) edge between i
  and j if and only if variables i and j are conditionally dependent
  given S for all possible subsets S of the remaining nodes; (ii) a directed
  edge \eqn{i \longrightarrow j}{i → j} means that this directed edge is 
  present in all DAGs in the Markov equivalence class; (iii) an undirected 
  edge \eqn{i - j} means that there is at least one DAG in the Markov 
  equivalence class with edge \eqn{i \longrightarrow j}{i → j} and
  there is at least one DAG in the Markov equivalence class with edge 
  \eqn{i \longleftarrow j}{i ← j}.

  The CPDAG is estimated using the PC algorithm (named after its inventors
  \bold{P}eter Spirtes and \bold{C}lark Glymour). The skeleton is
  estimated by the function \code{\link{skeleton}} which uses a modified
  version of the original PC algorithm (see Colombo and Maathuis (2014) for
  details). The original PC algorithm is known to be
  order-dependent, in the sense that the output depends on the order in
  which the variables are given. Therefore, Colombo and Maathuis (2014)
  proposed a simple modification, called PC-stable, that yields
  order-independent adjacencies in the skeleton (see the help file
  of this function for details). Subsequently, as many edges as possible
  are oriented. This is done in two steps. It is important to note that
  if no further actions are taken (see below) these two steps still
  remain order-dependent.

  The edges are oriented as follows. First, the algorithm considers all
  triples \code{(a,b,c)}, where \eqn{a} and \eqn{b} are adjacent, \eqn{b} and 
  \eqn{c} are adjacent, but \eqn{a} and \eqn{c} are not adjacent. For all such 
  triples, we direct both edges towards \eqn{b} 
  (\eqn{a \longrightarrow b \longleftarrow c}{a → b ← c}) if and only if 
  \eqn{b} was not part of the conditioning set that made the edge between 
  \eqn{a} and \eqn{c} drop out. These conditioning sets were saved in
  \code{sepset}. The structure 
  \eqn{a \longrightarrow b \longleftarrow c}{a → b ← c} is called a 
  v-structure.

  After determining all v-structures, there may still
  be undirected edges. It may be possible to direct some of these edges, since
  one can deduce that one of the two possible directions of the edge is
  invalid because it introduces
  a new v-structure or a directed cycle. Such edges are found by
  repeatedly applying rules R1-R3 of the PC algorithm as given in
  Algorithm 2 of Kalisch and Bühlmann (2007). The algorithm stops if
  none of the rules is applicable to the graph.

  The conservative PC algorithm (\code{conservative = TRUE}) is a
  slight variation of the PC algorithm (see Ramsey et al. 2006). After
  the skeleton is computed, all potential v-structures \eqn{a - b - c} are 
  checked in the following way. We test whether a and c are independent
  conditioning on all subsets of the neighbors of \eqn{a} and all subsets of 
  the neighbors of \eqn{c}. When a subset makes \eqn{a} and \eqn{c} 
  conditionally independent, we call it a separating set. If \eqn{b} is in no 
  such separating set or in all such separating sets, no further action is 
  taken and the usual PC is continued. If, however, \eqn{b} is in only some 
  separating sets, the triple \eqn{a - b - c} is marked as 'ambiguous'.
  Moreover, if no separating set is found among the neighbors, the triple is 
  also marked as 'ambiguous'. An ambiguous triple is not oriented as a
  v-structure. Furthermore, no further orientation rule that needs to
  know whether \eqn{a - b - c} is a v-structure or not is applied. Instead of
  using the conservative version, which is quite strict towards the
  v-structures, Colombo and Maathuis (2014) introduced a less strict
  version for the v-structures called majority rule. This adaptation can
  be called using \code{maj.rule = TRUE}. In this case, the triple 
  \eqn{a - b - c} is marked as 'ambiguous' if and only if \eqn{b} is in 
  exactly 50 percent of such separating sets or no separating set was found. 
  If \eqn{b} is in less than 50 percent of the separating sets it is set as a 
  v-structure, and if in more than 50 percent it is set as a non v-structure 
  (for more details see Colombo and Maathuis, 2014). The usage of both the
  conservative and the majority rule versions resolve the
  order-dependence issues of the determination of the v-structures.

  Sampling errors (or hidden variables) can lead to conflicting
  information about edge directions. For example, one may find that
  \eqn{a - b - c} and \eqn{b - c - d} should both be directed as v-structures. 
  This gives conflicting information about the edge \eqn{b - c}, since it should 
  be directed as \eqn{b \longleftarrow c}{b ← c} in v-structure 
  \eqn{a \longrightarrow b \longleftarrow c}{a → b ← c}, while it should be 
  directed as \eqn{b \longrightarrow c}{b → c} in v-structure 
  \eqn{b \longrightarrow c \longleftarrow d}{b → c ← d}. With the option 
  \code{solve.confl = FALSE}, in such cases, we simply overwrite the
  directions of the conflicting edge. In the example above this means
  that we obtain 
  \eqn{a \longrightarrow b \longrightarrow c \longleftarrow d}{a → b → c ← d} 
  if \eqn{a - b - c} was visited first, and 
  \eqn{a \longrightarrow b \longleftarrow c \longleftarrow d}{a → b ← c ← d}
  if \eqn{b - c - d} was visited first, meaning that the final orientation on 
  the edge depends on the ordering in which the v-structures were
  considered. With the option \code{solve.confl = TRUE} (which is only
  supported with option \code{u2pd = "relaxed"}), we first generate a list
  of all (unambiguous) v-structures (in the example above \eqn{a - b - c} and
  \eqn{b - c - d}), and then we simply orient them allowing both directions on 
  the edge \eqn{b - c}, namely we allow the bi-directed edge 
  \eqn{b \leftrightarrow c}{b ↔ c} resolving the order-dependence issues on the 
  edge orientations. We denote bi-directed edges in the adjacency matrix 
  \eqn{M} of the graph as \code{M[b,c] = 2} and \code{M[c,b] = 2}. In a similar 
  way, using lists for the candidate edges for each orientation rule and 
  allowing bi-directed edges, the order-dependence issues in the orientation 
  rules can be resolved. Note that bi-directed edges merely represent a 
  conflicting orientation and they should not to be interpreted causally. The 
  useage of these lists for the candidate edges and allowing bi-directed edges 
  resolves the order-dependence issues on the orientation of the v-structures 
  and on the orientation rules, see Colombo and Maathuis (2014) for
  more details.

  Note that calling (\code{conservative = TRUE}), or \code{maj.rule =
  TRUE}, together with \code{solve.confl = TRUE} produces a fully
  order-independent output, see Colombo and Maathuis (2014).

  Sampling errors, non faithfulness, or hidden variables can also lead
  to non-extendable CPDAGs, meaning that there does not exist a DAG that
  has the same skeleton and v-structures as the graph found by the
  algorithm. An example of this is an undirected cycle consisting of the
  edges \eqn{a - b - c - d} and \eqn{d - a}. In this case it is impossible to 
  direct the edges without creating a cycle or a new v-structure. The option
  \code{u2pd} specifies what should be done in such a situation. If the
  option is set to \code{"relaxed"}, the algorithm simply outputs the
  invalid CPDAG. If the option is set to \code{"rand"}, all direction
  information is discarded and a random DAG is generated on the
  skeleton, which is then converted into its CPDAG. If the option is set
  to \code{"retry"}, up to 100 combinations of possible directions of
  the ambiguous edges are tried, and the first combination that results
  in an extendable CPDAG is chosen. If no valid combination is found, an
  arbitrary DAG is generated on the skeleton as in the option "rand",
  and then converted into its CPDAG. Note that the output can also be an
  invalid CPDAG, in the sense that it cannot arise from the oracle PC
  algorithm, but be extendible to a DAG, for example 
  \eqn{a \longrightarrow b \longleftarrow c \longleftarrow d}{a → b ← c ← d}. 
  In this case, \code{u2pd} is not used.

  Using the function \code{\link{isValidGraph}} one can check if the final output is indeed a valid CPDAG.
  
  Notes: (1) Throughout, the algorithm works with the column positions
  of the variables in the adjacency matrix, and not with the names of
  the variables. (2) When plotting the object, undirected and bidirected
  edges are equivalent.}

\references{
  D. Colombo and M.H. Maathuis (2014).Order-independent constraint-based
  causal structure learning. \emph{Journal of Machine Learning Research}
  \bold{15} 3741-3782. 

  M. Kalisch, M. Maechler, D. Colombo, M.H. Maathuis and P. Buehlmann
  (2012). Causal Inference Using Graphical Models with the R Package
  pcalg. \emph{Journal of Statistical Software} \bold{47(11)} 1--26,
  \url{http://www.jstatsoft.org/v47/i11/}.

  M. Kalisch and P. Buehlmann (2007).
  Estimating high-dimensional directed acyclic graphs with the PC-algorithm.
  \emph{JMLR} \bold{8} 613-636.

  J. Ramsey, J. Zhang and P. Spirtes (2006).
  Adjacency-faithfulness and conservative causal inference. In
  \emph{Proceedings of the 22nd Annual Conference on Uncertainty in
  Artificial Intelligence}. AUAI Press, Arlington, VA.

  P. Spirtes, C. Glymour and R. Scheines (2000).
  \emph{Causation, Prediction, and Search}, 2nd edition. The MIT Press.
}
\seealso{\code{\link{skeleton}} for estimating a skeleton of a DAG;
  \code{\link{udag2pdag}} for converting the
  skeleton to a CPDAG; \code{\link{gaussCItest}},
  \code{\link{disCItest}}, \code{\link{binCItest}} and
  \code{\link{dsepTest}} as examples for \code{indepTest}. \code{\link{isValidGraph}} for testing whether the output is a valid CPDAG.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}), Martin Maechler,
  and Diego Colombo.
}
\examples{
##################################################
## Using Gaussian Data
##################################################
## Load predefined data
data(gmG)
n <- nrow    (gmG8$ x)
V <- colnames(gmG8$ x) # labels aka node names

## estimate CPDAG
pc.fit <- pc(suffStat = list(C = cor(gmG8$x), n = n),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha=0.01, labels = V, verbose = TRUE)
if (require(Rgraphviz)) {
  ## show estimated CPDAG
  par(mfrow=c(1,2))
  plot(pc.fit, main = "Estimated CPDAG")
  plot(gmG8$g, main = "True DAG")
}
##################################################
## Using d-separation oracle
##################################################
## define sufficient statistics (d-separation oracle)
suffStat <- list(g = gmG8$g, jp = RBGL::johnson.all.pairs.sp(gmG8$g))
## estimate CPDAG
fit <- pc(suffStat, indepTest = dsepTest, labels = V,
          alpha= 0.01) ## value is irrelevant as dsepTest returns either 0 or 1
if (require(Rgraphviz)) {
  ## show estimated CPDAG
  plot(fit, main = "Estimated CPDAG")
  plot(gmG8$g, main = "True DAG")
}

##################################################
## Using discrete data
##################################################
## Load data
data(gmD)
V <- colnames(gmD$x)
## define sufficient statistics
suffStat <- list(dm = gmD$x, nlev = c(3,2,3,4,2), adaptDF = FALSE)
## estimate CPDAG
pc.D <- pc(suffStat,
           ## independence test: G^2 statistic
           indepTest = disCItest, alpha = 0.01, labels = V, verbose = TRUE)
if (require(Rgraphviz)) {
  ## show estimated CPDAG
  par(mfrow = c(1,2))
  plot(pc.D, main = "Estimated CPDAG")
  plot(gmD$g, main = "True DAG")
}

##################################################
## Using binary data
##################################################
## Load binary data
data(gmB)
V <- colnames(gmB$x)
## estimate CPDAG
pc.B <- pc(suffStat = list(dm = gmB$x, adaptDF = FALSE),
           indepTest = binCItest, alpha = 0.01, labels = V, verbose = TRUE)
pc.B
if (require(Rgraphviz)) {
  ## show estimated CPDAG
  plot(pc.B, main = "Estimated CPDAG")
  plot(gmB$g, main = "True DAG")
}

##################################################
## Detecting ambiguities due to sampling error
##################################################
## Load predefined data
data(gmG)
n <- nrow    (gmG8$ x)
V <- colnames(gmG8$ x) # labels aka node names

## estimate CPDAG
pc.fit <- pc(suffStat = list(C = cor(gmG8$x), n = n),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha=0.01, labels = V, verbose = TRUE)

## due to sampling error, some edges were overwritten:
isValidGraph(as(pc.fit, "amat"), type = "cpdag")

## re-fit with solve.confl = TRUE
pc.fit2 <- pc(suffStat = list(C = cor(gmG8$x), n = n),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha=0.01, labels = V, verbose = TRUE,
             solve.confl = TRUE)

## conflicting edge is V5 - V6
as(pc.fit2, "amat")
}%{examples}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
