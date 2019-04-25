\name{ages}
\alias{ages}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate an APDAG within the Markov equivalence class of a DAG using AGES
}
\description{
Estimate an APDAG (a particular PDAG) using the aggregative greedy equivalence search (AGES) algorithm, which uses the solution path of the greedy equivalence search (GES) algorithm of Chickering (2002).
}
\usage{
ages(data, lambda_min = 0.5 * log(nrow(data)), labels = NULL,
     fixedGaps = NULL, adaptive = c("none", "vstructures", "triples"),
     maxDegree = integer(0), verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A \eqn{n*p} matrix (or data frame) containing the observational data.
}
  \item{lambda_min}{
The smallest penalty parameter value used when computing the solution path of GES.
}
  \item{labels}{
Node labels; if NULL the names of the columns of the data matrix (or the names in the data frame) are used. If these are not specified the sequence 1 to p is used.
}
  \item{fixedGaps}{
{logical \emph{symmetric} matrix of dimension p*p.  If entry
    \code{[i, j]} is \code{TRUE}, the result is guaranteed to have no edge
    between nodes \eqn{i} and \eqn{j}.}
}
  \item{adaptive}{
{indicating whether constraints should be adapted to
    newly detected v-structures or unshielded triples (cf. details).}
}
  \item{maxDegree}{
{Parameter used to limit the vertex degree of the estimated
    graph.  Valid arguments:
    \enumerate{
      \item Vector of length 0 (default): vertex degree is not limited.
      \item Real number \eqn{r}, \eqn{0 < r < 1}: degree of vertex \eqn{v} is
        limited to \eqn{r \cdot n_v}, where \eqn{n_v} denotes the number of
        data points where \eqn{v} was not intervened.
      \item Single integer: uniform bound of vertex degree for all vertices
        of the graph.
      \item Integer vector of length \code{p}: vector of individual bounds
        for the vertex degrees.
    }
    }
}
  \item{verbose}{
{If \code{TRUE}, detailed output is provided.}
}
  \item{\dots}{
{Additional arguments for debugging purposes and fine tuning.}
}
}
\details{
This function tries to add orientations to the essential graph (CPDAG) found by \code{\link{ges}} (ran with lambda=lambda_min). It does it aggregating several CPDAGs present in the solution path of GES. Conceptually, AGES starts with the essential graph found by GES ran with \code{lambda = lambda_min}. Then, it checks for further (compatible) orientation information in other essential graphs present in the solution path of GES, i.e., in essential graphs outputted by GES for larger penalty parameters. With compatible we mean that the aggregation process is done such that the final APDAG is still within the Markov equivalence graph represented by the essential graph found by GES in the following sense: an APDAG can always be extended to a DAG without creating new v-structures. This DAG lies in the Markov equivalence class represented by the essential graph found by GES. The algorithm is explained in detail in Eigenmann, Nandy, and Maathuis (2017).

The arguments \code{fixedgaps} and \code{adaptive} work also with AGES. However, they have not been studied in Eigenmann, Nandy, and Maathuis (2017).
Using the argument \code{fixedGaps}, one can make sure that certain edges
  will \emph{not} be present in the resulting essential graph: if the entry
  \code{[i, j]} of the matrix passed to \code{fixedGaps} is \code{TRUE}, there
  will be no edge between nodes \eqn{i} and \eqn{j}. The argument \code{adaptive} can be
  used to relax the constraints encoded by \code{fixedGaps} according to a 
  modification of GES called ARGES (adaptively restricted greedy 
  equivalence search) which has been presented in Nandy, Hauser and Maathuis
  (2018):
  \itemize{
    \item When \code{adaptive = "vstructures"} and the algorithm introduces a 
    new v-structure \eqn{a \longrightarrow b \longleftarrow c}{a → b ← c} in the 
    forward phase, then the edge \eqn{a - c} is removed from the list of fixed 
    gaps, meaning that the insertion of an edge between \eqn{a} and \eqn{c} 
    becomes possible even if it was forbidden by the initial matrix passed to 
    \code{fixedGaps}.
    
    \item When \code{adaptive = "triples"} and the algorithm introduces a new
    unshielded triple in the forward phase (i.e., a subgraph of three nodes
    \eqn{a}, \eqn{b} and \eqn{c} where \eqn{a} and \eqn{b} as well as \eqn{b}
    and \eqn{c} are adjacent, but \eqn{a} and \eqn{c} are not), then the edge
    \eqn{a - c} is removed from the list of fixed gaps.
  }
  With one of the adaptive modifications, the successive application of a 
  skeleton estimation method and GES restricted to an estimated skeleton still
  gives a \emph{consistent} estimator of the DAG, which is not the case without
  the adaptive modification.

For a detailed explanation of the GES function as well as its related object like essential graphs, we refer to the \code{\link{ges}} function.

Differences in the arguments with respect to GES: AGES uses \code{data} to initialize several scores taken as argument by GES. AGES modifies the forward and backward phases of GES performing single steps in either directions. For this reason, \code{phase}, \code{iterate}, and \code{turning} are not available arguments.
}
\value{
  \code{ages} returns a list with the following four components:
  \item{essgraph}{An object of class \code{\linkS4class{EssGraph}} containing an
    estimate of the equivalence class of the underlying DAG.}
  \item{repr}{An object of a class derived from \code{\linkS4class{ParDAG}}
    containing a (random) representative of the estimated equivalence class.}
  \item{CPDAGsList}{A list of p*p matrices containing all CPDAGs considered by AGES in the aggregation processes}
  \item{lambda}{A vector containing the penalty parameter used to obtain the list of CPDAGs mentioned above. GES returns the list of CPDAGs when used with this vector of penalty parameters if used with phases = c("forward", "backward") and iterate = FALSE.}
}
\references{
  D.M. Chickering (2002).  Optimal structure identification with greedy search.
  \emph{Journal of Machine Learning Research} \bold{3}, 507--554

  M.F. Eigenmann, P. Nandy, and M.H. Maathuis (2017). Structure learning of linear Gaussian structural equation models with weak edges. In \emph{Proceedings of UAI 2017}
  
  P. Nandy, A. Hauser and M.H. Maathuis (2018). High-dimensional consistency in score-based and hybrid structure learning. \emph{Annals of Statistics}, to appear.
}
\author{
  Marco Eigenmann (\email{eigenmann@stat.math.ethz.ch})
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{ges}}, \code{\linkS4class{EssGraph}}
}
\examples{

## Example 1: ages adds correct orientations: Bar --> V6 and Bar --> V8

set.seed(77)

p <- 8
n <- 5000
## true DAG:
vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
gGtrue <- randomDAG(p, prob = 0.3, V = vars)
data = rmvDAG(n, gGtrue)


## Estimate the aggregated PDAG with ages
ages.fit <- ages(data = data)


## Estimate the essential graph with ges
## We specify the phases in order to have a fair comparison of the algorithms
## Without the phases specified it would be easy to find examples
## where each algorithm outperforms the other
score <- new("GaussL0penObsScore", data)
ges.fit <- ges(score, phase = c("forward","backward"), iterate = FALSE)

## Plots
par(mfrow=c(1,3))
plot(ges.fit$essgraph, main="Estimated CPDAG with GES")
plot(ages.fit$essgraph, main="Estimated APDAG with AGES")
plot(gGtrue, main="TrueDAG")



## Example 2: ages adds correct orientations: Author --> Goal and Author --> V5

set.seed(50)

p <- 9
n <- 5000
## true DAG:
vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:9))
gGtrue <- randomDAG(p, prob = 0.5, V = vars)
data = rmvDAG(n, gGtrue)


## Estimate the aggregated PDAG with ages
ages.fit <- ages(data = data)


## Estimate the essential graph with ges
## We specify the phases in order to have a fair comparison of the algorithms
## Without the phases specified it would be easy to find examples
## where each algorithm outperforms the other
score <- new("GaussL0penObsScore", data)
ges.fit <- ges(score, phase = c("forward","backward"), iterate = FALSE)

## Plots
par(mfrow=c(1,3))
plot(ges.fit$essgraph, main="Estimated CPDAG with GES")
plot(ages.fit$essgraph, main="Estimated APDAG with AGES")
plot(gGtrue, main="TrueDAG")


## Example 3: ges and ages return the same graph

data(gmG)

data <- gmG8$x

## Estimate the aggregated PDAG with ages
ages.fit <- ages(data = data)


## Estimate the essential graph with ges
score <- new("GaussL0penObsScore", data)
ges.fit <- ges(score)


## Plots
par(mfrow=c(1,3))
plot(ges.fit$essgraph, main="Estimated CPDAG with GES")
plot(ages.fit$essgraph, main="Estimated APDAG with AGES")
plot(gmG8$g, main="TrueDAG")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{graphs}% __ONLY ONE__ keyword per line
