\name{pcSelect.presel}
\alias{pcSelect.presel}
\title{PC-Select preselection: Estimate subgraph around a response variable using preselection}
\description{
  This function uses pcSelect to preselect some covariates and then runs
  pcSelect again on the reduced data set.
}
\usage{
pcSelect.presel(y, dm, alpha, alphapre, corMethod = "standard", verbose = 0, directed=FALSE)
}
\arguments{
  \item{y}{Response Vector (length(y)=nrow(dm))}
  \item{dm}{Data matrix (rows: samples, cols: nodes)}
  \item{alpha}{Significance level of individual partial correlation
    tests}
  \item{alphapre}{Significance level for pcSelect in preselection}
  \item{corMethod}{"standard" or "Qn" for standard or robust correlation
  estimation}
  \item{verbose}{0-no output, 1-small output, 2-details (using 1 and 2
  makes the function very much slower)}
  \item{directed}{Boolean; should the output graph be directed?}
}
\value{
  \item{pcs}{A boolean vector indicating which column of \code{dm} is
  associated with \code{y}}
  \item{zMin}{The minimal z-values when testing partial correlations
  between \code{y} and each column of \code{dm}. The larger the number,
  the more consistent is the edge with the data.}
  \item{Xnew}{Preselected Variables.}
}
\details{
  First, \code{pcSelect} is run using \code{alphapre}. Then,
  only the important variables are kept and \code{pcSelect} is run on
  them again.
}
\seealso{\code{\link{pcSelect}}
}
\author{
  Philipp Ruetimann.
}
\examples{
p <- 10
## generate and draw random DAG :
set.seed(101)
myDAG <- randomDAG(p, prob = 0.2)
plot(myDAG, main = "randomDAG(10, prob = 0.2)")

## generate 1000 samples of DAG using standard normal error distribution
n <- 1000
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

## let's pretend that the 10th column is the response and the first 9
## columns are explanatory variable. Which of the first 9 variables
## "cause" the tenth variable?
y <- d.mat[,10]
dm <- d.mat[,-10]
res <- pcSelect.presel(d.mat[,10],d.mat[,-10],alpha=0.05,alphapre=0.6)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
