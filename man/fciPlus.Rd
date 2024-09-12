\name{fciPlus}
\alias{fciPlus}
\title{Estimate a PAG with the FCI+ Algorithm}
\description{
  Estimate a Partial Ancestral Graph (PAG) from observational data, using
  the FCI+ (Fast Causal Inference) algorithm, or from a combination of
  data from different (e.g., observational and interventional) contexts, 
  using the FCI+-JCI (Joint Causal Inference) extension.
}
\usage{
fciPlus(suffStat, indepTest, alpha, labels, p, verbose=TRUE,
	selectionBias = TRUE, jci = c("0","1","12","123"), contextVars = NULL)
}
\arguments{
  \item{suffStat}{sufficient statistics: A named \code{\link{list}}
    containing all necessary elements for the conditional independence
    decisions in the function \code{indepTest}.}
  \item{indepTest}{a \code{\link{function}} for testing conditional
    independence.  The function is internally called as
    \code{indepTest(x,y, S, suffStat)}, and tests conditional independence
    of \code{x} and \code{y} given \code{S}.  Here, \code{x} and
    \code{y} are variables, and \code{S} is a (possibly empty) vector of
    variables (all variables are denoted by their column numbers
    in the adjacency matrix). \code{suffStat} is a list with all
    relevant information, see above.  The return value of
    \code{indepTest()} is the p-value of the test for conditional
    independence.}
  \item{alpha}{numeric significance level (in \eqn{(0, 1)}) for the
    individual conditional independence tests.}
  \item{labels}{(optional) \code{\link{character}} vector of variable (or
    \dQuote{node}) names.  Typically preferred to specifying \code{p}.}
  \item{p}{(optional) number of variables (or nodes).  May be specified
    if \code{labels} are not, in which case \code{labels} is set to
    \code{1:p}.}% FIXME: can get *both* from  suffStat$C  in the typical case
  \item{selectionBias}{If \code{TRUE}, allow for selection bias. 
    If \code{FALSE}, selection bias is excluded by assumption and hence 
    rules R5-R7, as in Zhang (2008), are disabled.}
  \item{jci}{String specifying the JCI assumptions that are used. 
    It can be one of:
    \describe{
      \item{"0"}{No JCI assumption is made (default),}
      \item{"1"}{JCI assumption 1 (no system variable causes any context 
	variable),}
      \item{"12"}{JCI assumptions 1 and 2 (no system variable causes any 
        context variable, and no system variable is confounded with any 
	context variable),}
      \item{"123"}{JCI assumptions 1, 2 and 3 (no system variable causes
	any context variable, no system variable is confounded with any 
	context variable, and all context variables are confounded but 
	are not direct causes of each other).}
    }
    For more information, see Mooij et al. (2020).
    }
  \item{contextVars}{Subset of variable indices \{1,...,p\} that will be 
    treated as context variables in the JCI extension of FCI+.}
  \item{verbose}{logical indicating if progress of the algorithm should
    be printed.  The default is true, which used to be hard coded previously.}
}
\details{
  A (possibly much faster) variation of FCI (Fast Causal Inference).  
  For details, please see the references, and also \code{\link{fci}}.
}
\value{An object of \code{\link{class}} \code{fciAlgo} (see
 \code{\linkS4class{fciAlgo}}) containing the estimated graph
 (in the form of an adjacency matrix with various possible edge marks),
 the conditioning sets that lead to edge removals (sepset) and several other
 parameters.
}
\references{
  T. Claassen, J. Mooij, and T. Heskes (2013).
  Learning Sparse Causal Models is not NP-hard.
  In \emph{UAI 2013, Proceedings of the 29th Conference on Uncertainty in Artificial Intelligence}
}
\author{
  Emilija Perkovic, Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
  and Joris Mooij.
}
\seealso{
  \code{\link{fci}} for estimating a PAG using the FCI algorithm.
}
\examples{
##################################################
## Example without latent variables
##################################################

## generate a random DAG ( p = 7 )
set.seed(42)
p <- 7
myDAG <- randomDAG(p, prob = 0.4)

## find PAG using the FCI+ algorithm on "Oracle"
suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
m.fci <- fciPlus(suffStat, indepTest=gaussCItest,
                 alpha = 0.9999, p=p)
summary(m.fci)

if (require(Rgraphviz)) {
  par(mfrow=c(1,2))
plot(myDAG)
plot(m.fci)
}

##################################################
## Joint Causal Inference Example
## Mooij et al. (2020), Fig. 43(a), p. 97
##################################################

# Encode MAG as adjacency matrix
p <- 8 # total number of variables
V <- c("Ca","Cb","Cc","X0","X1","X2","X3","X4") # 3 context variables, 5 system variables
# amat[i,j] = 0 iff no edge btw i,j
# amat[i,j] = 1 iff i *-o j
# amat[i,j] = 2 iff i *-> j
# amat[i,j] = 3 iff i *-- j
amat <- rbind(c(0,2,2,2,0,0,0,0),
              c(2,0,2,0,2,0,0,0),
              c(2,2,0,0,2,2,0,0),
              c(3,0,0,0,0,0,2,0),
              c(0,3,3,0,0,3,0,2),
              c(0,0,3,0,2,0,0,0),
              c(0,0,0,3,0,0,0,2),
              c(0,0,0,0,2,0,3,0))
rownames(amat)<-V
colnames(amat)<-V

# Make use of d-separation oracle as "independence test"
indepTest <- dsepAMTest
suffStat<-list(g=amat,verbose=FALSE)

# Derive PAG that represents the Markov equivalence class of the MAG with the FCI+ algorithm
# (assuming no selection bias)
# fci.pag <- fciPlus(suffStat, indepTest, alpha = 0.5, labels = V,
# selectionBias=FALSE,verbose=TRUE)

# Derive PAG with FCI+-JCI, the Joint Causal Inference extension of FCI
# (assuming no selection bias, and all three JCI assumptions)
# fcijci.pag <- fciPlus(suffStat, indepTest, alpha = 0.5, labels = V,
# selectionBias=FALSE, contextVars=c(1,2,3), jci="123", verbose=TRUE)

# Report results
# cat('True MAG:\n')
# print(amat)
# cat('PAG output by FCI+:\n')
# print(fci.pag@amat)
# cat('PAG output by FCI+-JCI:\n')
# print(fcijci.pag@amat)

# Read off causal features from the FCI PAG
#cat('Identified absence (-1) and presence (+1) of ancestral causal relations from FCI+ PAG:\n')
#print(pag2anc(fci.pag@amat))
#cat('Identified absence (-1) and presence (+1) of direct causal relations from FCI+ PAG:\n')
#print(pag2edge(fci.pag@amat))
#cat('Identified absence (-1) and presence (+1) of pairwise latent confounding from FCI+ PAG:\n')
#print(pag2conf(fci.pag@amat))

# Read off causal features from the FCI-JCI PAG
#cat('Identified absence (-1) and presence (+1) of ancestral causal relations from FCI+-JCI PAG:\n')
#print(pag2anc(fcijci.pag@amat))
#cat('Identified absence (-1) and presence (+1) of direct causal relations from FCI+-JCI PAG:\n')
#print(pag2edge(fcijci.pag@amat))
#cat('Ident. absence (-1) and presence (+1) of pairwise latent confounding from FCI+-JCI PAG:\n')
#print(pag2conf(fcijci.pag@amat))
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
