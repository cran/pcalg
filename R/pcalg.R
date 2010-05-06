###  MM: (ess-set-style 'DEFAULT)  -- for indenting

check.Rgraphviz <- function() {
  if(!require("Rgraphviz"))
    stop("Package 'Rgraphviz' (from Bioconductor) must be installed for plotting graphs!")
}

##' .. content for \description{} (no empty lines) ..
##'
##' @title
##' @param g a "graph" or adjacency matrix
##' @return
trueCov <- function(g) {
  if (is(g, "graphNEL")) {
    w <- wgtMatrix(g)
  } else {
    w <- g
  }
  ## find all parents
  p <- ncol(w)
  pa <- vector("list",p)
  for (i in 1:p) pa[[i]] <- which(w[i,]!=0)
  l.pa <- sapply(pa,length)

  ## make EX-Matrix
  ecov <- diag(1,p)
  for (j in 1:(p-1)) {
    for (i in (j+1):p) {
      if (l.pa[i]>0) {
        ecov[j,i]  <- sum(w[i,pa[[i]]]*ecov[j,pa[[i]]])
      } else {
        ecov[j,i] <- 0
      }
    }
  }

  ## make Covariance Matrix
  xcov <- matrix(0,p,p)

  for (i in 1:p) {
    for (j in 1:i) {
      case.id <- as.character((l.pa[i]!=0)*10+(l.pa[j]!=0))
      switch(case.id,
             "0" = {
               ## l.pa[i]=0 & l.pa[j]=0
               if (i!=j) {
                 xcov[i,j] <- xcov[j,i] <- 0
               } else {
                 xcov[i,j] <- 1
               }},
             "1" = {
               ## l.pa[i]=0 & l.pa[j] != 0
               xcov[i,j] <- xcov[j,i] <- ecov[i,j]
             },
             "10" = {
               ## l.pa[i]!=0 & l.pa[j] = 0
               xcov[i,j] <- xcov[j,i] <- ecov[j,i]
             },
             "11" = {
               ## l.pa[i]!=0 & l.pa[j] = 0
               xcov[i,j] <- xcov[j,i] <-
                 w[i,pa[[i]],drop=FALSE] %*% xcov[pa[[i]],pa[[j]]] %*% t(w[j,pa[[j]],drop=FALSE]) +
                   ecov[j,i]
             })
    }
  }
  xcov
}

randomDAG <- function (n, prob, lB = 0.1, uB = 1)
{
    stopifnot(n >= 2, is.numeric(prob), length(prob) == 1, 0 <= prob, prob <= 1,
              is.numeric(lB), is.numeric(uB), lB <= uB)
    V <- as.character(1:n)
    edL <- as.list(V)
    names(edL) <- V
    nmbEdges <- 0
    for (i in seq(length = n - 2)) {
        listSize <- rbinom(1, n - i, prob)
        nmbEdges <- nmbEdges + listSize
        edgeList <- sample(seq(i + 1, n), size = listSize)
        weightList <- runif(length(edgeList), min = lB, max = uB)
        edL[[i]] <- list(edges = edgeList, weights = weightList)
    }
    if (rbinom(1, 1, prob) == 1) {
        edL[[n - 1]] <- list(edges = n, weights = runif(1, min = lB,
                                        max = uB))
    }
    else {
        edL[[n - 1]] <- list(edges = integer(0), weights = numeric(0))
    }
    edL[[n]] <- list(edges = integer(0), weights = numeric(0))
    if (nmbEdges > 0) {
        res <- new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
    } else {
        res <- new("graphNEL", nodes = V, edgemode = "directed")
    }
    res
}

## A version of this is also in	/u/maechler/R/MM/Pkg-ex/graph/weightmatrix.R
## another on in  Matrix/R/sparseMatrix.R  function graph.wgtMatrix() :

wgtMatrix <- function(g, transpose = TRUE)
{
  ## Purpose: work around "graph" package's  as(g, "matrix") bug
  ## ----------------------------------------------------------------------
  ## ACHTUNG: mat_[i,j]==1 iff j->i,
  ## whereas with as(g,"matrix") mat_[i,j]==1 iff i->j
  ## ----------------------------------------------------------------------
  ## Arguments: g: an object inheriting from (S4) class "graph"
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, based on Seth Falcon's code;  Date: 12 May 2006

  ## MM: another buglet for the case of  "no edges":
  if(numEdges(g) == 0) {
    p <- length(nd <- nodes(g))
    return( matrix(0, p,p, dimnames = list(nd, nd)) )
  }
  ## Usual case, when there are edges:
  if(!("weight" %in% names(edgeDataDefaults(g))))
    edgeDataDefaults(g, "weight") <- 1L
  w <- unlist(edgeData(g, attr = "weight"))
  ## we need the *transposed* matrix typically:
  tm <- if(transpose) t(as(g, "matrix")) else as(g, "matrix")
  ## now is a 0/1 - matrix (instead of 0/wgts) with the 'graph' bug
  if(any(w != 1)) ## fix it
    tm[tm != 0] <- w
  ## tm_[i,j]==1 iff i->j
  tm
}

rmvDAG <- function(n, dag, errDist = c("normal", "cauchy", "mix", "mixt3", "mixN100","t4"),
                   mix = 0.1, errMat = NULL)
{
  ## Purpose: Generate data according to a given DAG (with weights) and
  ## given node distribution (rows: number of samples; cols: node values in
  ## topological order)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - n     : Number of samples
  ## - dag   : Graph object containing the DAG and weights
  ## - errDist: "normal" or "mix" for pure standard normal node distribution
  ##           or mixing with standard cauchy
  ## - mix   : Percentage of points sampled from standard cauchy
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006;  Martin Maechler

  ## check input &  initialize variables
  stopifnot(is(dag, "graph"),
            (p <- length(nodes(dag))) >= 2)

  ##  as(.,"matrix") now {for some versions of 'graph' pkg} is 0/1
  ## weightMatrix <- t(as(dag,"matrix"))
  weightMatrix <- wgtMatrix(dag)

  ## check if top. sorted
  nonZeros <- which(weightMatrix != 0, arr.ind = TRUE)
  if (nrow(nonZeros)>0) {
    if (any((nonZeros[,1] - nonZeros[,2])<0) ||
        any(diag(weightMatrix) != 0) )
      stop("Input DAG must be topologically ordered!")
  }

  errDist <- match.arg(errDist)

  if(is.null(errMat)) {
    ## generate errors e_i
    errMat <-
      switch(errDist,
             "normal" = matrix(rnorm  (n*p),  nrow = n),
             "cauchy" = matrix(rcauchy(n*p),  nrow = n),
             "t4" =     matrix(rt(n*p, df=4), nrow = n),
             "mix" = {
               cauchySamples <- rcauchy(round(mix*n*p))
               normSamples <- rnorm(n*p-length(cauchySamples))

               samples <- c(normSamples,cauchySamples)
               matrix(sample(samples,length(samples)), nrow = n)
             },
             "mixt3" = {
               t3Samples <- rt(round(mix*n*p),df=3)
               normSamples <- rnorm(n*p-length(t3Samples))

               samples <- c(normSamples,t3Samples)
               matrix(sample(samples,length(samples)), nrow = n)
             },
             "mixN100" = {
               outlierSamples <- rnorm(round(mix*n*p),sd=10)
               normSamples <- rnorm(n*p-length(outlierSamples))

               samples <- c(normSamples,outlierSamples)
               matrix(sample(samples,length(samples)), nrow = n)
             })
  }
  else { ## check & use 'errMat' argument:
    stopifnot(!is.null(dim.eM <- dim(errMat)),
              dim.eM == c(n,p), is.numeric(errMat))
  }
  X <- matrix(0, n, p) # will contain the result

  ## compute X matrix X_i
  if (sum(weightMatrix) > 0) {
    X[,1] <- errMat[,1] # for (one of) the root node(s)
    for (j in 2:p) { ## uses X[*, 1:(j-1)] -- "recursively" !
      ij <- 1:(j-1)
      X[,j] <- X[, ij, drop = FALSE] %*% weightMatrix[j, ij] + errMat[,j]
    }
    X
  }
  else { ## return
    errMat
  }
}


pcSelect <- function(y,dm, alpha, corMethod = "standard", verbose = FALSE, directed=FALSE)
{
  ## Purpose: Find columns in dm, that have nonzero parcor with y given
  ## any other set of columns in dm
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - y: Response Vector (length(y)=nrow(dm))
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - verbose: 0-no output, 1-small output, 2-details
  ## ----------------------------------------------------------------------
  ## Value: List
  ## - G: boolean vector with connected nodes
  ## - zMin: Minimal z values
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 27.4.07

  stopifnot((n <- nrow(dm)) >= 1,
            (p <- ncol(dm)) >= 1)
  vNms <- colnames(dm)
  cl <- match.call()

  zMin <- c(0,rep.int(Inf,p))
  C <- mcor(cbind(y,dm), method = corMethod)
  cutoff <- qnorm(1 - alpha/2)
  n.edgetests <- numeric(1)# final length = max { ord}
  ## G := complete graph :
  G <- c(FALSE,rep.int(TRUE,p))
  seq_p <- 1:(p+1)

  done <- FALSE
  ord <- 0
  while (!done && any(G)) {
    n.edgetests[ord+1] <- 0
    done <- TRUE
    ind <- which(G)
    remainingEdgeTests <- length(ind)
    if(verbose>=1)
      cat("Order=",ord,"; remaining edges:",remainingEdgeTests,"\n", sep='')
    for (i in 1:remainingEdgeTests) {
      if(verbose && i%%100==0) cat("|i=",i,"|iMax=",nrow(ind),"\n")
      y <- 1
      x <- ind[i]

      if (G[x]) {
        nbrsBool <- G
        nbrsBool[x] <- FALSE
        nbrs <- seq_p[nbrsBool]
        ## neighbors of y without itself and x
        length_nbrs <- length(nbrs)

        if (length_nbrs >= ord) {
          if (length_nbrs > ord) done <- FALSE
          S <- seq(length = ord)

          ## now includes special cases (ord == 0) or (length_nbrs == 1):
          repeat {
            n.edgetests[ord+1] <- n.edgetests[ord+1]+1
            z <- zStat(x,y, nbrs[S], C,n)
            if(abs(z)<zMin[x]) zMin[x] <- abs(z)
            if (verbose >= 2)
              cat(paste("x:",vNms[x-1],"y:",(ytmp <- round((p+1)/2)),"S:"),
                  c(ytmp,vNms)[nbrs[S]],paste("z:",z,"\n"))
            if (abs(z) <= cutoff) {
              G[x] <- FALSE
              break
            }
            else {
              nextSet <- getNextSet(length_nbrs, ord, S)
              if(nextSet$wasLast)
                break
              S <- nextSet$nextSet
            }
          } ## {repeat}
        }
      } ## end if( G )
    } ## end for(i ..)
    ord <- ord+1
  } ## end while
  Gres <- G[-1]
  names(Gres) <- vNms
  list(G = Gres, zMin = zMin[-1])
}## pcSelect

zStat <- function(x,y, S, C, n)
{
  ## Purpose: Fisher's z-transform statistic of partial corr.(x,y | S)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x,y,S: Are x,y cond. indep. given S?
  ## - C: Correlation matrix among nodes
  ## - n: Samples used to estimate correlation matrix
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, 22 May 2006; Markus Kalisch

  ## for C use:
  ## dyn.load("/u/kalisch/cCode/pcAlgo/parcorC.so")
  ##  res <- 0

  ##   if (length(S) < 4) {
  r <- pcorOrder(x,y, S, C)
  ##  } else {
  ##    k <- solve(C[c(x,y,S),c(x,y,S)])
  ##    r <- -k[1,2]/sqrt(k[1,1]*k[2,2])
  ##      r <- .C("parcorC",as.double(res),as.integer(x-1),as.integer(y-1),as.integer(S-1),as.integer(length(S)),as.integer(dim(C)[1]),as.double(as.vector(C)))[[1]]
  ##  }

  res <- sqrt(n- length(S) - 3) * ( 0.5*log( (1+r)/(1-r) ) )
  if (is.na(res)) res <- 0
  res
}

condIndFisherZ <- function(x,y,S,C,n, cutoff,
                           verbose= isTRUE(getOption("verbose.pcalg.condIFz")))
{
  ## Purpose: Return boolean result on conditional independence using
  ## Fisher's z-transform
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x,y,S: Are x,y cond. indep. given S?
  ## - C: Correlation matrix among nodes
  ## - n: Samples used to estimate correlation matrix
  ## - cutoff: Cutoff for significance level for individual
  ##           partial correlation tests
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006, 17:32

  ## R Variante
  r <- pcorOrder(x,y, S, C)

  ## C Variante
  ## C res <- 0
  ## C p <- dim(C)[1]
  ## C r <- .C("parcor",as.double(res),as.integer(x-1),as.integer(y-1),as.integer(S-1),as.integer(length(S)),as.integer(p),as.double(C))[[1]]

  z <- 0.5*log( (1+r)/(1-r) )
  T <- sqrt(n-length(S)-3)*z
  ## cat(" (",x,",",y,") | ",S," : T = ",T,"\n", sep='')

  ## MM: T is only NA when 'r' already is (r = +- 1  <==>  T = +- Inf) -- better check there (FIXME?)
  ## is.na(T) <==>  T <- 0 # if problem, delete edge: be conservative
  is.na(T) || abs(T) <= cutoff
}

pcorOrder <- function(i,j, k, C, cut.at = 0.9999999) {
  ## Purpose: Compute partial correlation
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - i,j,k: Partial correlation of i and j given k
  ## - C: Correlation matrix among nodes
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler
  if (length(k)==0) {
    r <- C[i,j]
  } else if (length(k)==1) {
    r <- (C[i, j] - C[i, k] * C[j, k])/sqrt((1 - C[j, k]^2) * (1 - C[i, k]^2))
  } else { ## length(k) >= 2
    ## If corpcor was decent and had a name space, we'd just use
    ## PM <- corpcor::pseudoinverse(C[c(i,j,k), c(i,j,k)])
    stopifnot(require("corpcor", quietly=TRUE))
    PM <- pseudoinverse(C[c(i,j,k), c(i,j,k)])
    r <- -PM[1, 2]/sqrt(PM[1, 1] * PM[2, 2])
  }
  if(is.na(r)) r <- 0
  min(cut.at, max(-cut.at, r))
}


compareGraphs <- function(gl,gt) {
  ## Purpose: Return TPR, FPR and TDR of comparison of two undirected graphs
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gl: Estimated graph (may be directed, but the direction will
  ##       be dropped internally)
  ## - gt: True graph (may be directed, but the direction will
  ##       be dropped internally)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006, 17:35

  ## When 'graph' returns the 'weight matrix' again:
  ##   ml <- as(ugraph(gl), "matrix")
  ##   mt <- as(ugraph(gt), "matrix")
  ##   p <- dim(ml)[1]
  ml <- wgtMatrix(ugraph(gl))
  mt <- wgtMatrix(ugraph(gt))
  p <- dim(ml)[2]

  mt[mt != 0] <- rep(1,sum(mt != 0))
  ml[ml != 0] <- rep(1,sum(ml != 0)) ## inserted to fix bug

  ## FPR: #misplaced edges/#true gaps
  diffm <- ml-mt
  nmbTrueGaps <- (sum(mt == 0)-p)/2
  ##  print(list(p=p,sum=sum(mt==0),mt=mt,ml=ml,nmbTrueGaps=nmbTrueGaps,diffm=diffm))
  fpr <- if (nmbTrueGaps == 0) 1 else (sum(diffm > 0)/2)/nmbTrueGaps

  ## TPR: #correctly found edges/#true edges
  diffm2 <- mt-ml
  nmbTrueEdges <- (sum(mt == 1)/2)
  tpr <- if (nmbTrueEdges == 0) 0 else 1 - (sum(diffm2 > 0)/2)/nmbTrueEdges

  ## TDR: #correctly found edges/#found edges
  trueEstEdges <- (nmbTrueEdges-sum(diffm2 > 0)/2) ## #true edges-#not detected
  tdr <-
    if (sum(ml == 1) == 0) { ## no edges detected
      if (trueEstEdges == 0) 1 ## no edges in true graph
      else 0
    } else trueEstEdges/(sum(ml == 1)/2)

  ## return named vector:
  c(tpr = tpr, fpr = fpr, tdr = tdr)
}

getNextSet <- function(n,k,set) {
  ## Purpose: Generate the next set in a list of all possible sets of size
  ##          k out of 1:n;
  ##  Also returns a boolean whether this set was the last in the list.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - n,k: Choose a set of size k out of numbers 1:n
  ## - set: previous set in list
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006, 17:37

  ## chInd := changing Index
  chInd <- k - (zeros <- sum((seq(n-k+1,n)-set) == 0))
  wasLast <- (chInd == 0)
  if (!wasLast) {
    set[chInd] <- set[chInd] + 1
    if (chInd < k)
      set[(chInd+1):k] <- seq(set[chInd]+1, set[chInd]+zeros)
  }
  list(nextSet = set, wasLast = wasLast)
}

mcor <- function(dm, method =
                 c("standard", "Qn", "QnStable", "ogkScaleTau2", "ogkQn","shrink"))
{
  ## Purpose: Compute correlation matrix (perhaps elementwise)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix; rows: samples, cols: variables
  ## - method: "Qn" or "standard" (default) envokes robust (based on Qn
  ##           scale estimator) or standard correlation estimator, respectively.
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 27 Jan 2006, 10:07

  p <- ncol(dm)
  method <- match.arg(method)
  switch(method,
	 "Qn" = {
	   res <- matrix(0, p,p)
	   qnVec <- apply(dm, 2, Qn)

	   for (i in 1:(p-1)) {
	     for (j in i:p) {
	       qnSum <- Qn(dm[,i] + dm[,j])
	       qnDiff <- Qn(dm[,i] - dm[,j])
	       res[i,j] <- max(-1,
                               min(1,
                                   (qnSum^2 - qnDiff^2) / (4*qnVec[i]*qnVec[j])))
	     }
	   }
	   res <- res + t(res)
	   diag(res) <- rep(1,p)
	   res
	 },
         "QnStable" = {
	   res <- matrix(0, p,p)
	   qnVec <- apply(dm, 2, Qn)

	   for (i in 1:(p-1)) {
	     for (j in i:p) {
	       qnSum <- Qn(dm[,i]/Qn(dm[,i]) + dm[,j]/Qn(dm[,j]))
	       qnDiff <- Qn(dm[,i]/Qn(dm[,i]) - dm[,j]/Qn(dm[,j]))
	       res[i,j] <- max(-1,
                               min(1,
                                   (qnSum^2 - qnDiff^2) / (qnSum^2 + qnDiff^2)))
	     }
	   }
	   res <- res + t(res)
	   diag(res) <- rep(1,p)
	   res
	 },
	 "ogkScaleTau2" = {
	   cov2cor(covOGK(dm, n.iter = 2, sigmamu = scaleTau2,
			  weight.fn = hard.rejection)$cov)
	 },
	 "ogkQn" = {
	   cov2cor(covOGK(dm, n.iter = 2, sigmamu = s_Qn,
			  weight.fn = hard.rejection)$cov)
	 },
	 "standard" = cor(dm),
         "shrink"={
           X <- dm
           CM <- cor(X)
           n <- nrow(X)
           p <- ncol(X)
           g <- 0
           S1 <- 0
           S2 <- 0
           scor3 <- CM

           for(k in 2:p){
             S1 <- S1+(CM[1,k])^2
             S2 <- S2+(1-(CM[1,k])^2)^2
           }

           g <-2*S1/(2*S1+(2/n)*S2)

           for(i in 2:p){
             scor3[1,i] <- g*CM[1,i]
             scor3[i,1] <- g*CM[i,1]
           }

           scor3
         }
         )
}

pcSelect.presel <- function(y,dm, alpha, alphapre, corMethod = "standard", verbose = 0, directed=FALSE)
{
  ## Purpose: Find columns in dm, that have nonzero parcor with y given
  ## any other set of columns in dm with use of some kind of preselection
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - y: Response Vector (length(y)=nrow(dm))
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - alphapre: Significance level in preselective use of pcSelect
  ## ----------------------------------------------------------------------
  ## Author: Philipp Ruetimann, Date: 5 Mar 2008

  tmp <- pcSelect(y,dm,alphapre,corMethod,verbose,directed)
  n <- nrow(dm)
  pcs <- tmp$G
  zmi <- tmp$zMin
  ppcs <- which(pcs==1)
  Xnew <- dm[,ppcs]
  if(length(ppcs)==1){Xnew <- matrix(Xnew,n,1)}
  lang <- length(pcs)
  tmp2 <- pcSelect(y,Xnew,alpha,corMethod,verbose,directed)
  pcSnew <- tmp2$G
  zminew <- tmp2$zMin
  zmi[ppcs] <- zminew
  ## MM FIXME -- do without for() :
  k <- 1
  for (i in 1:lang){
    if(pcs[i]==1){
      pcs[i] <- pcSnew[k]
      k <- k+1
    }
  }

  list(pcs=pcs,Xnew=Xnew,zMin=zmi)
}


corGraph <- function(dm, alpha = 0.05, Cmethod = "pearson")
{
  ## Purpose: Computes a correlation graph. Significant correlations are
  ## shown using the given correlation method.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix with rows as samples and cols as variables
  ## - alpha: Significance level for correlation test
  ## - Cmethod: a character string indicating which correlation coefficient
  ##          is to be  used for the test.  One of '"pearson"',
  ##          '"kendall"', or '"spearman"', can be abbreviated.
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 30 Jan 2006, 12:08

  stopifnot(is.numeric(p <- ncol(dm)), p >= 2)
  mat <- matrix(0, p,p)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      mat[i,j] <- cor.test(dm[,i], dm[,j], alternative = "two.sided",
                           method = Cmethod)$p.value < alpha
    }
  }
  mat <- mat + t(mat)
  rownames(mat) <- colnames(mat) <- 1:p
  as(mat, "graphNEL")
}

##################################################
## CPDAG
##################################################

##################################################
## dag2cpdag
##################################################

make.edge.df <- function(amat) {
  ## Purpose: Generate a data frame describing some properties of a DAG
  ## (for extending to a CPDAG)
  ## The output contains xmin,xmax,head,tail,order (NA or number),
  ## type (1="d",0="u") in lexikographic order
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: Adjacency matrix of DAG [x_ij=1 means i->j]
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:43

  ## INPUT: Adjacency matrix
  stopifnot(sum(amat)>0)
  e <- which(amat==1,arr.ind=TRUE)
  e.dup <- duplicated(t(apply(e,1,sort)))
  nmb.edges <- sum(!e.dup)
  res <- data.frame(xmin=rep(NA,nmb.edges),xmax=rep(NA,nmb.edges),
                    tail=rep(NA,nmb.edges),head=rep(NA,nmb.edges),
                    order=rep(NA,nmb.edges),type=rep(1,nmb.edges))
  pure.edges <- e[!e.dup,]
  if(length(pure.edges)==2) dim(pure.edges) <- c(1,2)
  for (i in 1:dim(pure.edges)[1]) {
    if (all(amat[pure.edges[i,1],pure.edges[i,2]]==
            amat[pure.edges[i,2],pure.edges[i,1]])) {
      res$type[i] <- 0
      res$head[i] <- NA
      res$tail[i] <- NA
    } else {
      res$head[i] <- pure.edges[i,2]
      res$tail[i] <- pure.edges[i,1]
    }
  }
  s.pure.edges <- t(apply(pure.edges,1,sort))
  ii <- order(s.pure.edges[,1],s.pure.edges[,2])
  res <- res[ii,]
  res$xmin <- s.pure.edges[ii,1]
  res$xmax <- s.pure.edges[ii,2]
  res
}

orderEdges <- function(amat) {
  ## Purpose: Order the edges of a DAG according to Chickering
  ## (for extension to CPDAG)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: Adjacency matrix of DAG
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:42

  stopifnot(isAcyclic(amat))
  ordered.nodes <- topOrder(amat) ## parents before children
  edge.df <- make.edge.df(amat)

  eOrder <- 0
  while(any(unOrdered <- is.na(edge.df$order))) {
    counter <- 0
    ## find y
    y <- NA
    found <- FALSE
    while(!found) {
      counter <- counter+1
      node <- ordered.nodes[counter]
      ## which edges are incident to node?
      nbr.nodes <- which(amat[,node]==1)
      if(length(nbr.nodes)>0) {
        unlabeled <- rep.int(FALSE, length(nbr.nodes))
        for(i in seq_along(nbr.nodes)) {
          x <- nbr.nodes[i]
          ## is edge edge x-y unlabeled?
	  unlabeled[i] <- length(intersect(which(edge.df$xmin==min(node,x) &
						 edge.df$xmax==max(node,x)),
					   which(unOrdered))) > 0
        }
        ## choose unlabeled edge with highest order node
        if(any(unlabeled)) {
          nbr.unlab <- nbr.nodes[unlabeled] #nbrnodes w. unlabeled edges
          tmp <- ordered.nodes[ordered.nodes %in% nbr.unlab]
          y <- tmp[length(tmp)]
          ## y <- last(ordered.nodes[which(ordered.nodes %in% nbr.unlab)])
          edge.df$order[edge.df$xmin==min(node,y) &
                        edge.df$xmax==max(node,y)] <- eOrder
          eOrder <- eOrder+1
          found <- TRUE
        }
      }

    } ## while !found

  } ## while any(unOrdered)
  edge.df
}


labelEdges <- function(amat) {
  ## Purpose: Label the edges in a DAG with "compelled" and "reversible"
  ## (for extension to a CPDAG)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: Adjacency matrix of DAG
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006;  prettified: MMaechler

  ## label=TRUE -> compelled
  ## label=FALSE -> reversible
  edge.df <- orderEdges(amat)
  lab <- rep(NA,dim(edge.df)[1])
  edge.df <- edge.df[order(edge.df$order),]
  Head <- edge.df$head
  Tail <- edge.df$tail

  while(any(ina <- is.na(lab))) {
    x.y <- which(ina)[1]
    x <- Tail[x.y]
    y <- Head[x.y]
    y.is.head <- Head == y
    e1 <- which(Head == x & lab)
    for(ee in e1) {
      w <- Tail[ee]
      if (any(wt.yh <- w == Tail & y.is.head))
        lab[wt.yh] <- TRUE
      else {
        lab[y.is.head] <- TRUE
        break
      }
    }
    ## edges going to y not starting from x
    cand <- which(y.is.head  &  Tail != x)
    if (length(cand) > 0) {
      valid.cand <- rep(FALSE,length(cand))
      for (iz in seq_along(cand)) {
        z <- Tail[cand[iz]]
        if (!any(Tail==z & Head==x)) ## NOT.parent.of.x :
          valid.cand[iz] <- TRUE
      }
      cand <- cand[valid.cand]
    }
    lab[which(y.is.head & is.na(lab))] <- (length(cand) > 0)
  }
  edge.df$label <- lab
  edge.df
}


dag2cpdag <- function(dag) {
  ## Purpose: Compute the (unique) completed partially directed graph (CPDAG)
  ## that corresponds to the input DAG; result is a graph object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dag: input DAG (graph object)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:30

  p <- numNodes(dag)
  ## transform DAG to adjacency matrix if any edges are present
  if (numEdges(dag)==0) {
    cpdag.res <- dag
  } else {
    dag <- as(dag,"matrix")
    dag[dag!=0] <- 1

    ## dag is adjacency matrix
    e.df <- labelEdges(dag)
    cpdag <- matrix(0, p,p)
    for (i in 1:dim(e.df)[1]) {
      if (e.df$label[i]) {
        cpdag[e.df$tail[i],e.df$head[i]] <- 1
      } else {
        cpdag[e.df$tail[i],e.df$head[i]] <- cpdag[e.df$head[i],e.df$tail[i]] <- 1
      }
    }
    rownames(cpdag) <- colnames(cpdag) <- as.character(seq(1,p))
    cpdag.res <- as(cpdag,"graphNEL")
  }
  cpdag.res
}

##################################################
## pdag2dag
##################################################
find.sink <- function(gm) {
  ## Purpose: Find sink of an adj matrix; return numeric(0) if there is none;
  ## a sink may have incident undirected edges, but no directed ones
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: Adjacency matrix (gm_i_j is edge from j to i)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:28

  ## treat undirected edges
  undir.nbrs <- which(gm == t(gm) & gm == 1, arr.ind=TRUE)
  gm[undir.nbrs] <- 0
  ## treat directed edges
  which(colSums(gm) == 0)
}

adj.check <- function(gm,x) {
  ## Purpose:  Return "TRUE", if:
  ## For every vertex y, adj to x, with (x,y) undirected, y is adjacent to
  ## all the other vertices which are adjacent to x
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: adjacency matrix of graph
  ## - x: node number (number)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:27

  res <- TRUE
  ## find undirected neighbors of x
  un <- which(gm[x,]==1 & gm[,x]==1)
  if (length(un)>0) {
    for (i in seq_along(un)) {
      y <- un[i]
      adj.x <- setdiff(which(gm[x,]==1 | gm[,x]==1),y)
      adj.y <- setdiff(which(gm[y,]==1 | gm[,y]==1),x)
      axINay <- all(adj.x %in% adj.y)
      res <- (res & axINay)
    }
  }
  res
}


amat2dag <- function(amat) {
  ## Purpose: Transform the adjacency matrix of an PDAG to the adjacency
  ## matrix of a DAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: adjacency matrix; x -> y if amat[x,y]=1,amat[y,x]=0
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:23

  p <- dim(amat)[1]
  ## amat to skel
  skel <- amat+t(amat)
  skel[which(skel>1)] <- 1

  ## permute skel
  ord <- sample(1:p)
  skel <- skel[ord,ord]

  ## skel to dag
  for (i in 2:p) {
    for (j in 1:(i-1)) {
      if(skel[i,j]==1) skel[i,j] <- 0
    }
  }
  ## inverse permutation
  i.ord <- order(ord)
  skel[i.ord,i.ord]
}

##################################################
## udag2pdag
##################################################
udag2pdag <- function(gInput, verbose=FALSE) {
  ## Purpose: Transform the Skeleton of a pcAlgo-object to a PDAG using
  ## the rules of Pearl. The output is again a pcAlgo-object.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gInput: pcAlgo object
  ## - verbose: 0 - no output, 1 - detailed output
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Sep 2006, 15:03

  res <- gInput
  if (numEdges(gInput@graph)>0) {
    g <- as(gInput@graph,"matrix") ## g_ij if i->j
    p <- dim(g)[1]
    pdag <- g
    ind <- which(g==1,arr.ind=TRUE)

    ## Create minimal pattern
    for (i in 1:dim(ind)[1]) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(g[y,]==1),x) ## x-y-z
      for (z in allZ) {
        if (g[x,z]==0  &&
            !(y %in% gInput@sepset[[x]][[z]] ||
              y %in% gInput@sepset[[z]][[x]])) {
          if (verbose) {
            cat("\n",x,"->",y,"<-",z,"\n")
            cat("Sxz=",gInput@sepset[[z]][[x]],"Szx=",gInput@sepset[[x]][[z]])
          }
          pdag[x,y] <- pdag[z,y] <- 1
          pdag[y,x] <- pdag[y,z] <- 0
        }
      }
    }

    ## Test whether this pdag allows a consistent extension
    res2 <- pdag2dag(as(pdag,"graphNEL"))

    if (res2$success) {
      ## Convert to complete pattern: use rules by Pearl
      old_pdag <- matrix(0, p,p)
      while (!all(old_pdag == pdag)) {
        old_pdag <- pdag
        ## rule 1
        ind <- which((pdag==1 & t(pdag)==0), arr.ind=TRUE) ## a -> b
        for (i in seq_len(nrow(ind))) {
          a <- ind[i,1]
          b <- ind[i,2]
          indC <- which( (pdag[b,]==1 & pdag[,b]==1) & (pdag[a,]==0 & pdag[,a]==0))
          if (length(indC)>0) {
            pdag[b,indC] <- 1
            pdag[indC,b] <- 0
            if (verbose)
              cat("\nRule 1:",a,"->",b," and ",b,"-",indC,
                  " where ",a," and ",indC," not connected: ",b,"->",indC,"\n")
          }
        }
        ## x11()
        ## plot(as(pdag,"graphNEL"), main="After Rule1")

        ## rule 2
        ind <- which((pdag==1 & t(pdag)==1), arr.ind=TRUE) ## a -> b
        for (i in seq_len(nrow(ind))) {
          a <- ind[i,1]
          b <- ind[i,2]
          indC <- which( (pdag[a,]==1 & pdag[,a]==0) & (pdag[,b]==1 & pdag[b,]==0))
          if (length(indC)>0) {
            pdag[a,b] <- 1
            pdag[b,a] <- 0
            if (verbose) cat("\nRule 2: Kette ",a,"->",indC,"->",
                  b,":",a,"->",b,"\n")
          }
        }
        ## x11()
        ## plot(as(pdag,"graphNEL"), main="After Rule2")

        ## rule 3
        ind <- which((pdag==1 & t(pdag)==1), arr.ind=TRUE) ## a - b
        for (i in seq_len(nrow(ind))) {
          a <- ind[i,1]
          b <- ind[i,2]
          indC <- which( (pdag[a,]==1 & pdag[,a]==1) & (pdag[,b]==1 & pdag[b,]==0))
          if (length(indC)>=2) {
            ## cat("R3: indC = ",indC,"\n")
            g2 <- pdag[indC,indC]
            ## print(g2)
            if (length(g2)<=1) {
              g2 <- 0
            } else {
              diag(g2) <- rep(1,length(indC)) ## no self reference
            }
            if (any(g2==0)) { ## if two nodes in g2 are not connected
              pdag[a,b] <- 1
              pdag[b,a] <- 0
              if (verbose) cat("\nRule 3:",a,"->",b,"\n")
            }
          }
        }
        ## x11()
        ## plot(as(pdag,"graphNEL"), main="After Rule3")

        ## rule 4
        ##-         ind <- which((pdag==1 & t(pdag)==1), arr.ind=TRUE) ## a - b
        ##-         if (length(ind)>0) {
        ##-           for (i in 1:dim(ind)[1]) {
        ##-             a <- ind[i,1]
        ##-             b <- ind[i,2]
        ##-             indC <- which( (pdag[a,]==1 & pdag[,a]==1) & (pdag[,b]==0 & pdag[b,]==0))
        ##-             l.indC <- length(indC)
        ##-             if (l.indC>0) {
        ##-               found <- FALSE
        ##-               ic <- 0
        ##-               while(!found & (ic < l.indC)) {
        ##-                 ic <- ic + 1
        ##-                 c <- indC[ic]
        ##-                 indD <- which( (pdag[c,]==1 & pdag[,c]==0) & (pdag[,b]==1 & pdag[b,]==0))
        ##-                 if (length(indD)>0) {
        ##-                   found <- TRUE
        ##-                   pdag[b,a] = 0
        ##-                   if (verbose) cat("Rule 4 applied \n")
        ##-                 }
        ##-               }
        ##-             }
        ##-           }
        ##-         }

      }
      res@graph <- as(pdag,"graphNEL")
    } else {
      ## was not extendable; random DAG chosen
      res@graph <- res2$graph
      ## convert to CPDAG
      res@graph <- dag2cpdag(res@graph)
    }
  }
  return(res)
} ## udag2pdag

shd <- function(g1,g2)
{
  ## Purpose: Compute Structural Hamming Distance between graphs g1 and g2
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - g1, g2: Input graphs
  ## (graph objects;connectivity matrix where m[x,y]=1 iff x->1
  ## and m[x,y]=m[y,x]=1 iff x-y; pcAlgo-objects)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  1 Dec 2006, 17:21

                                        # Idea: Transform g1 into g2
                                        # Transform g1 and g2 into adjacency matrices
  if (is(g1, "pcAlgo")) g1 <- g1@graph
  if (is(g2, "pcAlgo")) g2 <- g2@graph

  if (is(g1, "graphNEL")) {
    m1 <- wgtMatrix(g1, transp=FALSE)
    m1[m1 != 0] <- 1
  }
  if (is(g2, "graphNEL")) {
    m2 <- wgtMatrix(g2, transp=FALSE)
    m2[m2 != 0] <- 1
  }

  p <- dim(m1)[2]
  shd <- 0
                                        # Remove superfluous edges from g1
  s1 <- m1 + t(m1)
  s2 <- m2 + t(m2)
  s1[s1==2] <- 1
  s2[s2==2] <- 1
  ds <- s1 - s2
  ind <- which(ds > 0)
  m1[ind] <- 0
  shd <- shd + length(ind)/2
                                        # Add missing edges to g1
  ind <- which(ds < 0)
  m1[ind] <- m2[ind]
  shd <- shd + length(ind)/2
                                        # Compare Orientation
  d <- abs(m1-m2)
  ## return
  shd + sum((d + t(d)) > 0)/2
}

################################################################################
## New in V8 ; uses  vcd  package
################################################################################
ci.test <- function(x,y, S=NULL, dm.df) {
  stopifnot(is.data.frame(dm.df), ncol(dm.df) > 1)
  tab <- table(dm.df[,c(x,y,S)])
  if (any(dim(tab) < 2))
    1
  else if (length(S)==0)
    fisher.test(tab, simulate.p.value=TRUE)$p.value
  else
    vcd::coindep_test(tab,3:(length(S)+2))$p.value
}

pcAlgo <- function(dm = NA, C = NA, n=NA, alpha, corMethod = "standard",
                   verbose = FALSE, directed=FALSE,
                   G=NULL, datatype='continuous', NAdelete=TRUE,
                   m.max=Inf, u2pd="rand", psepset=FALSE) {
  ## Purpose: Perform PC-Algorithm, i.e., estimate skeleton of DAG given data
  ## Output is an unoriented graph object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - C: correlation matrix (only for continuous)
  ## - n: sample size
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - G: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical)
  ## - datatype: distinguish between discrete and continuous data
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - u2pd: Function for converting udag to pdag
  ##   "rand": udag2pdagu
  ##   "relaxed": udag2pdagRelaxed
  ##   "retry": udag2pdagSpecial
  ## - psepset: Also check possible sep sets.
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler
  ## Modifications: Sarah Gerster, Date: July 2007
  ## Modifications: Diego Colombo, Date: Sept 2009

  .Deprecated(msg = "pcAlgo() is deprecated and only kept for backward compatibility.
 Please use skeleton, pc, or fci instead\n")
  cl <- match.call()

  if (any(is.na(dm))) {
    stopifnot(all(!is.na(C)),!is.na(n), (p <- ncol(C))>0)
  } else {
    n <- nrow(dm)
    p <- ncol(dm)
  }
  n <- as.integer(n)

  if (is.null(G)) {
    ## G := complete graph :
    G <- matrix(TRUE, p,p)
    diag(G) <- FALSE
  } else if (!(identical(dim(G),c(p,p))))
      stop("Dimensions of the dataset and G do not agree.")

  seq_p <- 1:p
  sepset <- pl <- vector("list",p)
  for (i in seq_p) sepset[[i]] <- pl
  zMin <- matrix(Inf, p,p)
  n.edgetests <- numeric(1)# final length = max { ord}
  done <- FALSE
  ord <- 0

  if (datatype=='continuous') {
    diag(zMin) <- rep(0,p)
    if (any(is.na(C))) C <- mcor(dm, method = corMethod)
    cutoff <- qnorm(1 - alpha/2)
    while (!done && any(G) && ord<=m.max) {
      n.edgetests[ord+1] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[,1]) ,]
      remainingEdgeTests <- nrow(ind)
      if(verbose)
        cat("Order=",ord,"; remaining edges:",remainingEdgeTests,"\n", sep='')
      for (i in 1:remainingEdgeTests) {
        if(verbose && i%%100 == 0) cat("|i=",i,"|iMax=",nrow(ind),"\n")
        x <- ind[i,1]
        y <- ind[i,2]
        if (G[y,x]) {
          nbrsBool <- G[,x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord) done <- FALSE
            S <- seq(length = ord)
            repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
              n.edgetests[ord+1] <- n.edgetests[ord+1]+1
              z <- zStat(x,y, nbrs[S], C,n)
              if (verbose) cat(paste("x:",x,"y:",y,"S:"),nbrs[S],paste("z:",z,"\n"))
              if(abs(z)<zMin[x,y]) zMin[x,y] <- abs(z)
              if (abs(z) <= cutoff) {
                G[x,y] <- G[y,x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              } else {
                nextSet <- getNextSet(length_nbrs, ord, S)
                if(nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              }
            }
          }
        } ## end if(!done)

      } ## end for(i ..)
      ord <- ord+1
      ##    n.edgetests[ord] <- remainingEdgeTests
    } ## while

    for (i in 1:(p-1)) {
      for (j in 2:p) {
        zMin[i,j] <- zMin[j,i] <- min(zMin[i,j],zMin[j,i])
      }
    }
#########
#########
#########
######### DISCRETE DATA ######################################################
  } else {
    if (datatype=='discrete') {
      dm.df <- as.data.frame(dm)
      while (!done && any(G) && ord<=m.max) {
        n.edgetests[ord+1] <- 0
        done <- TRUE
        ind <- which(G, arr.ind = TRUE)
        ## For comparison with C++ sort according to first row
        ind <- ind[order(ind[,1]) ,]
        remainingEdgeTests <- nrow(ind)
        if(verbose)
          cat("Order=",ord,"; remaining edges:",remainingEdgeTests,"\n", sep='')
        for (i in 1:remainingEdgeTests) {
          if(verbose) { if(i%%100==0) cat("|i=",i,"|iMax=",nrow(ind),"\n") }
          x <- ind[i,1]
          y <- ind[i,2]
          if (G[y,x]) {
            nbrsBool <- G[,x]
            nbrsBool[y] <- FALSE
            nbrs <- seq_p[nbrsBool]
            length_nbrs <- length(nbrs)
            if (length_nbrs >= ord) {
              if (length_nbrs > ord) done <- FALSE
              S <- seq(length = ord)
              repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
                n.edgetests[ord+1] <- n.edgetests[ord+1]+1
                prob <- ci.test(x,y, nbrs[S], dm.df)
                if (verbose) cat("x=",x," y=",y," S=",nbrs[S],":",prob,"\n")
                if (is.na(prob)) prob <- if(NAdelete) 1 else 0
                if(prob >= alpha) { # independent
                  G[x,y] <- G[y,x] <- FALSE
                  sepset[[x]][[y]] <- nbrs[S]
                  break
                } else {
                  nextSet <- getNextSet(length_nbrs, ord, S)
                  if(nextSet$wasLast)
                    break
                  S <- nextSet$nextSet
                }
              }
            }
          } ## end if(!done)

        } ## end for(i ..)
        ord <- ord+1
        ##    n.edgetests[ord] <- remainingEdgeTests
      } ## while
    } else {
      stop("Datatype must be 'continuous' or 'discrete'.")
    }
  }

  if (psepset) {
    amat <- G
    amat[amat==TRUE] <- 1
    amat[amat==FALSE] <- 0
    ind <- which(amat==1, arr.ind=TRUE)
    ## Orient colliders
    for (i in 1:dim(ind)[1]) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(amat[y,]==1),x) ## x-y-z

      for (z in allZ) {
        if (amat[x,z] == 0 &&
            !((y %in% sepset[[x]][[z]]) ||
              (y %in% sepset[[z]][[x]]))) {
          if (verbose >= 2) {
            cat("\n",x,"*->",y,"<-*",z,"\n")
            cat("Sxz=",sepset[[z]][[x]],"and","Szx=",sepset[[x]][[z]],"\n")
          }

          ## x o-> y <-o z
          amat[x,y] <- amat[z,y] <- 2

        } ## for
      } ## if
    } ## for

    ## Compute poss. sepsets
    for (x in 1:p) {
      attr(x,'class') <- 'possibledsep'
      if (any(amat[x,]!=0)){
        tf1 <- setdiff(reach(x,-1,-1,amat), x)
        for (y in seq_p[amat[x,]!=0]) {
          ## tf = possible_d_sep(amat,x,y)
          tf <- setdiff(tf1,y)
          ## test
          if (length(tf)>0) {
            az <- abs(zStat(x,y,tf,C,n))
            if (az < zMin[x,y]) zMin[x,y] <- az
            if (az <= cutoff) {
              ## delete x-y
              amat[x, y] <- amat[y, x] <- 0
              ## save pos d-sepset in sepset
              sepset[[x]][[y]] <- tf
            }
            if (verbose >= 2)
              cat("Possible-D-Sep of", x, "and", y, "is", tf, " - |z| = ",az,"\n")
          }
        }
      }
    }
    G[amat==0] <- FALSE
    G[amat==1] <- TRUE
  } ## end if(psepset)

  if(verbose) { cat("Final graph adjacency matrix:\n"); print(symnum(G)) }

  ## transform matrix to graph object :
  if (sum(G) == 0) {
    Gobject <- new("graphNEL", nodes = as.character(seq_p))
  } else {
    colnames(G) <- rownames(G) <- as.character(seq_p)
    Gobject <- as(G,"graphNEL")
  }

  res <- new("pcAlgo",
             graph = Gobject,
             call = cl, n = n, max.ord = as.integer(ord-1),
             n.edgetests = n.edgetests, sepset = sepset,
             zMin = zMin)
  if (directed) {
    res <- switch (u2pd,
                   "rand" = udag2pdag(res),
                   "retry" = udag2pdagSpecial(res)$pcObj,
                   "relaxed" = udag2pdagRelaxed(res))
  }
  res
}

flipEdges <- function(amat,ind) {
  res <- amat
  for (i in seq_len(nrow(ind))) {
    x <- ind[i,]
    res[x[1],x[2]] <- amat[x[2],x[1]]
    res[x[2],x[1]] <- amat[x[1],x[2]]
  }
  res
}

pdag2dag <- function(g, keepVstruct=TRUE) {
  ## Purpose: Generate a consistent extension of a PDAG to a DAG; if this
  ## is not possible, a random extension of the skeleton is returned and
  ## a warning is issued.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - g: PDAG (graph object)
  ## - keepVstruct: TRUE - vStructures are kept
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Sep 2006, 15:21

  if (numEdges(g)==0) {
    success <- TRUE
    graph <- g
  } else {
    gm <- wgtMatrix(g) ## gm_i_j is edge from j to i
    gm[which(gm>0 & gm!=1)] <- 1
    p <- dim(gm)[1]

    gm2 <- gm
    a <- gm
    go.on <- TRUE
    go.on2 <- FALSE
    while(go.on && length(a) > 1 && sum(a) > 0) {
      go.on <- FALSE
      go.on2 <- TRUE
      sinks <- find.sink(a)
      if (length(sinks)>0) {
        counter <- 1
        while(go.on2 && counter <= length(sinks)) {
          x <- sinks[counter]
          if (!keepVstruct || adj.check(a,x)) {
            go.on2 <- FALSE
            ## orient edges
            inc.to.x <- which(a[,x]==1 & a[x,]==1) ## undirected
            if (length(inc.to.x)>0) {
              real.inc.to.x <- as.numeric(row.names(a)[inc.to.x])
              real.x <- as.numeric(row.names(a)[x])
              gm2[real.x,real.inc.to.x] <- rep(1,length(inc.to.x))
              gm2[real.inc.to.x,real.x] <- rep(0,length(inc.to.x))
            }
            ## remove x and all edges connected to it
            a <- a[-x,-x]
          }
          counter <- counter+1
        }
      }
      go.on <- !go.on2
    }## { while }

    if (go.on2) {
      graph <- as(amat2dag(gm),"graphNEL")
      ## warning("PDAG not extendible: Random DAG on skeleton drawn")
      success <- FALSE
    } else {
      graph <- as(t(gm2), "graphNEL")
      success <- TRUE
    }
  }
  list(graph=graph, success=success)
}

udag2pdagSpecial <- function(gInput, verbose=FALSE, n.max=100) {
  ## Purpose: Transform the Skeleton of a pcAlgo-object to a PDAG using
  ## the rules of Pearl. The output is again a pcAlgo-object. Ambiguous
  ## v-structures are reoriented until extendable or max number of tries
  ## is reached. If still not extendable, a DAG is produced starting from the
  ## current PDAG even if introducing new v-structures.
  ##
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gInput: pcAlgo object
  ## - verbose: 0 - no output, 1 - detailed output
  ## - n.max: Maximal number of tries to reorient v-strucutres
  ## ----------------------------------------------------------------------
  ## Values:
  ## - pcObj: Oriented pc-Object
  ## - evisit: Matrix counting the number of orientation attemps per edge
  ## - xtbl.orig: Is original graph with v-structure extendable
  ## - xtbl: Is final graph with v-structure extendable
  ## - amat0: Adj.matrix of original graph with v-structures
  ## - amat1: Adj.matrix of graph with v-structures after reorienting
  ##          edges from double edge visits
  ## - status:
  ##   0: original try is extendable
  ##   1: reorienting double edge visits helps
  ##   2: orig. try is not extendable; reorienting double visits don't help;
  ##      result is acyclic, has orig. v-structures, but perhaps
  ##      additional v-structures
  ## - counter: Number of reorientation tries until success or max.tries
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Sep 2006, 15:03
  counter <- 0
  res <- gInput
  status <- 0
  p <- length(nodes(res@graph))
  evisit <- amat0 <- amat1 <- matrix(0,p,p)
  xtbl <- xtbl.orig <- TRUE
  if (numEdges(gInput@graph)>0) {
    g <- as(gInput@graph,"matrix") ## g_ij if i->j
    p <- dim(g)[1]
    pdag <- g
    ind <- which(g==1,arr.ind=TRUE)
    ## ind <- unique(t(apply(ind,1,sort)))

    ## Create minimal pattern
    for (i in 1:dim(ind)[1]) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(g[y,]==1),x) ## x-y-z
      for(z in allZ) {
        if ((g[x,z]==0) &&
            !(y %in% gInput@sepset[[x]][[z]] ||
              y %in% gInput@sepset[[z]][[x]])) {
          if (verbose) {
            cat("\n",x,"->",y,"<-",z,"\n")
            cat("Sxz=",gInput@sepset[[z]][[x]],"Szx=",gInput@sepset[[x]][[z]])
          }
          ## check if already in other direction directed
          if (pdag[x,y]==0 && pdag[y,x]==1) {
            evisit[x,y] <- evisit[x,y] + 1
            evisit[y,x] <- evisit[y,x] + 1
          }
          if (pdag[z,y]==0 && pdag[y,z]==1) {
            evisit[z,y] <- evisit[z,y] + 1
            evisit[y,z] <- evisit[y,z] + 1
          }
          pdag[x,y] <- pdag[z,y] <- 1
          pdag[y,x] <- pdag[y,z] <- 0
        } ## if
      } ## for
    } ## for ( i )

    amat0 <- pdag
    ## Test whether this pdag allows a consistent extension
    res2 <- pdag2dag(as(pdag,"graphNEL"))
    xtbl <- res2$success
    xtbl.orig <- xtbl

    if (!xtbl && (max(evisit)>0)) {
      tmp.ind2 <- unique(which(evisit>0,arr.ind=TRUE))
      ind2 <- unique(t(apply(tmp.ind2,1,sort)))
      ## print(ind2)
      n <- nrow(ind2)
      n.max <- min(2^n-1,n.max)
      counter <- 0
      ## xtbl is FALSE because of if condition
      while((counter<n.max) & !xtbl) {
        ## if (counter%%100 == 0) cat("\n counter=",counter,"\n")
        counter <- counter + 1
        dgBase <- digitsBase(counter)
        dgBase <- dgBase[length(dgBase):1]
        ## print(dgBase)
        indBase <- matrix(0,1,n)
        indBase[1,seq_along(dgBase)] <- dgBase
        ## indTmp <- ind2[ss[[counter]],,drop=FALSE]
        indTmp <- ind2[(indBase==1),,drop=FALSE]
        ## print(indTmp)
        pdagTmp <- flipEdges(pdag,indTmp)
        resTmp <- pdag2dag(as(pdagTmp,"graphNEL"))
        xtbl <- resTmp$success
      }
      pdag <- pdagTmp
      status <- 1
    }
    amat1 <- pdag

    if (xtbl) {
      ## Convert to complete pattern: use rules by Pearl
      old_pdag <- matrix(0, p,p)
      while (any(old_pdag != pdag)) {
        old_pdag <- pdag
        ## rule 1
        ind <- which((pdag==1 & t(pdag)==0), arr.ind=TRUE) ## a -> b
        if (length(ind)>0) {
          for (i in 1:dim(ind)[1]) {
            a <- ind[i,1]
            b <- ind[i,2]
            indC <- which( (pdag[b,]==1 & pdag[,b]==1) & (pdag[a,]==0 & pdag[,a]==0))
            if (length(indC)>0) {
              pdag[b,indC] <- 1
              pdag[indC,b] <- 0
              if (verbose)
                cat("\nRule 1:",a,"->",b," and ",b,"-",indC,
                    " where ",a," and ",indC," not connected: ",b,"->",indC,"\n")
            }
          }
          ## x11()
          ## plot(as(pdag,"graphNEL"), main="After Rule1")
        }

        ## rule 2
        ind <- which((pdag==1 & t(pdag)==1), arr.ind=TRUE) ## a -> b
        if (length(ind)>0) {
          for (i in 1:dim(ind)[1]) {
            a <- ind[i,1]
            b <- ind[i,2]
            indC <- which( (pdag[a,]==1 & pdag[,a]==0) & (pdag[,b]==1 & pdag[b,]==0))
            if (length(indC)>0) {
              pdag[a,b] <- 1
              pdag[b,a] <- 0
              if (verbose) cat("\nRule 2: Kette ",a,"->",indC,"->",
                    b,":",a,"->",b,"\n")
            }
          }
        }
        ## x11()
        ## plot(as(pdag,"graphNEL"), main="After Rule2")

        ## rule 3
        ind <- which((pdag==1 & t(pdag)==1), arr.ind=TRUE) ## a - b
        if (length(ind)>0) {
          for (i in 1:dim(ind)[1]) {
            a <- ind[i,1]
            b <- ind[i,2]
            indC <- which( (pdag[a,]==1 & pdag[,a]==1) & (pdag[,b]==1 & pdag[b,]==0))
            if (length(indC)>=2) {
              ## cat("R3: indC = ",indC,"\n")
              g2 <- pdag[indC,indC]
              ## print(g2)
              if (length(g2)<=1) {
                g2 <- 0
              } else {
                diag(g2) <- rep(1,length(indC)) ## no self reference
              }
              if (any(g2==0)) { ## if two nodes in g2 are not connected
                pdag[a,b] <- 1
                pdag[b,a] <- 0
                if (verbose) cat("\nRule 3:",a,"->",b,"\n")
              }
            }
          }
        }
      }
      res@graph <- as(pdag,"graphNEL")
    } else {
      res@graph <- dag2cpdag(pdag2dag(as(pdag,"graphNEL"),keepVstruct=FALSE)$graph)
      status <- 2
      ## res@graph <- res2$graph
    }
  }
  list(pcObj=res, evisit=evisit, xtbl=xtbl, xtbl.orig=xtbl.orig,
       amat0=amat0, amat1=amat1, status=status, counter=counter)
}

udag2pdagRelaxed <- function(gInput, verbose=FALSE, unfVect=NULL)
{
  ## Purpose: Transform the Skeleton of a pcAlgo-object to a PDAG using
  ## the rules of Pearl. The output is again a pcAlgo-object. There is
  ## NO CHECK whether the resulting PDAG is really extendable.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gInput: pcAlgo object
  ## - verbose: 0 - no output, 1 - detailed output
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Sep 2006, 15:03

  if (numEdges(gInput@graph) == 0) return(gInput)
  ## else

  ## collider
  g <- as(gInput@graph, "matrix")
  p <- dim(g)[1]
  pdag <- g
  ind <- which(g==1, arr.ind=TRUE)

  ## Create minimal pattern
  for (i in seq_len(nrow(ind))) {
    x <- ind[i,1]
    y <- ind[i,2]
    allZ <- setdiff(which(g[y,]==1),x) ## x-y-z

    for (z in allZ) {
      if (g[x,z] == 0 &&
          !((y %in% gInput@sepset[[x]][[z]]) ||
            (y %in% gInput@sepset[[z]][[x]]))) {
        ## normal version
        if (length(unfVect)==0) {
          if (verbose)
            cat("\n", x, "->", y, "<-", z, "\n",
                "Sxz=", gInput@sepset[[z]][[x]],
                "Szx=", gInput@sepset[[x]][[z]])
          pdag[x, y] <- pdag[z, y] <- 1
          pdag[y, x] <- pdag[y, z] <- 0
        }
        ## conservative version
        else {
          ## check if x-y-z is faithful
          if (!any(unfVect == triple2numb(p,x,y,z)) &&
              !any(unfVect == triple2numb(p,z,y,x))) {
            ## faithful -> as usual; otherwise do nothing
            if (verbose) {
              cat("\n", x, "->", y, "<-", z, "\n")
              cat("Sxz=", gInput@sepset[[z]][[x]], "Szx=",
                  gInput@sepset[[x]][[z]])
            }
            pdag[x, y] <- pdag[z, y] <- 1
            pdag[y, x] <- pdag[y, z] <- 0
          }
        }
      }
    }
  } ## for ( i )

  repeat {
    old_pdag <- pdag

    ## Rule 1
    ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- ((pdag[b, ] == 1 & pdag[, b] == 1) &
              (pdag[a, ] == 0 & pdag[, a] == 0))
      if (any(isC)) {
        indC <- which(isC)
        ## normal version
        if (length(unfVect)==0) {
          pdag[b, indC] <- 1
          pdag[indC, b] <- 0
          if (verbose)
            cat("\nRule 1:", a, "->", b, " and ", b,
                "-", indC, " where ", a, " and ", indC,
                " not connected: ", b, "->", indC,
                "\n")
        }
        ## conservative version
        else {
          for (c in indC) {
            ## check that a-b-c not unfaithful
            if (!any(unfVect == triple2numb(p, a,b,c)) &&
                !any(unfVect == triple2numb(p, c,b,a))) {
              pdag[b, c] <- 1
              pdag[c, b] <- 0
              if (verbose)
                cat("\nRule 1':", a, "->", b, " and ", b, "-", c,
                    " where ", a, " and ", c, " not connected and ",
                    a, b, c," faithful triple: ", b, "->", c, "\n")
            }
          }
        }
      }
    }##end for (i ..)

    ## Rule 2
    ## normal version = conservative version
    ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- ((pdag[a, ] == 1 & pdag[, a] == 0) &
              (pdag[, b] == 1 & pdag[b, ] == 0))
      if (any(isC)) {
        pdag[a, b] <- 1
        pdag[b, a] <- 0
        if (verbose)
          cat("\nRule 2: Chain ", a, "->", which(isC),
              "->", b, ":", a, "->", b, "\n")
      }
    }

    ## Rule 3
    ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      indC <- which((pdag[a, ] == 1 & pdag[, a] == 1) &
                    (pdag[, b] == 1 & pdag[b, ] == 0))
      if (length(indC) >= 2) {
        ## normal version
        if (length(unfVect)==0) {
          g2 <- pdag[indC, indC]
          if (length(g2) <= 1) {
            g2 <- 0
          }
          else {
            diag(g2) <- rep(1, length(indC))
          }
          if (any(g2 == 0)) {
            pdag[a, b] <- 1
            pdag[b, a] <- 0
            if (verbose)
              cat("\nRule 3:", a, "->", b, "\n")
          }
        }
        ## conservative version
        else {
          cmb.C <- combn(indC, 2)
          cC1 <- cmb.C[1,]
          cC2 <- cmb.C[2,]
          for (j in seq_along(cC1)) {
            c1 <- cC1[j]
            c2 <- cC2[j]
            if (c1 != c2 && pdag[c1,c2] == 0 && pdag[c2,c1] == 0 &&
                !any(unfVect == triple2numb(p, c1,a,c2)) &&
                !any(unfVect == triple2numb(p, c2,a,c1))) { ## faithful triple found
              pdag[a, b] <- 1
              pdag[b, a] <- 0
              if (verbose)
                cat("\nRule 3':", a, c1, c2, "faithful triple: ", a, "->", b, "\n")
              break
              ##=== out of loop
            }
          }## for ( j ..)
        }
      }## if( length(indC) >= 2 )
    }## for ( i ..)

    if(all(pdag == old_pdag)) break
  }## repeat

  gInput@graph <- as(pdag, "graphNEL")
  gInput
}


## DEPRECATED! -- use  ida() --
beta.special <- function(dat=NA,x.pos,y.pos,verbose=0,a=0.01,myDAG=NA,myplot=FALSE,perfect=FALSE,
                         method="local",collTest=TRUE,pcObj=NA,all.dags=NA,u2pd="rand")
{
  ## Purpose: Estimate the causal effect of x on y; the pcObj and all DAGs
  ## can be precomputed
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dat: data
  ## - x.pos, y.pos: Column of x and y in d.mat
  ## - verbose: 0=no comments, 1=progress in BB, 2=detail on estimates
  ## - a: significance level of tests for finding CPDAG
  ## - myDAG: needed if bootstrp==FALSE
  ## - myplot: plot estimated graph
  ## - perfect: True cor matrix is calculated from myDAG
  ## - method: "local" - local (all combinations of parents in regr.)
  ##           "global" - all DAGs
  ## - collTest: True - Exclude orientations of undirected edges that
  ##   introduce a new collider
  ## - pcObj: Fit of PC Algorithm (semidirected); if this is available, no
  ##   new fit is done
  ## - all.dags: All DAGs in the format of function allDags; if this is
  ##   available, no new function call allDags is done
  ## - u2pd: Function for converting udag to pdag
  ##   "rand": udag2pdag
  ##   "relaxed": udag2pdagRelaxed
  ##   "retry": udag2pdagSpecial
  ## ----------------------------------------------------------------------
  ## Value: causal values
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 21 Nov 2007, 11:18

  ##  dat=d.mat;x.pos;y.pos;verbose=0;a=0.01;myDAG=NA;myplot=FALSE;perfect=FALSE;method="global";collTest=TRUE;pcObj=NA;all.dags=NA;trueLocal=TRUE;u2pd="rand"

  cat("This function is deprecated and is only kept for backward compatibility.
Please use ida or idaFast instead\n")

  tmpColl <- FALSE

  ## Covariance matrix: Perfect case / standard case
  if (perfect) {
    if(!is(myDAG, "graphNEL")) stop("For perfect-option the true DAG is needed!")
    mcov <- trueCov(myDAG)
    mcor <- cov2cor(mcov)
  } else {
    mcov <- cov(dat)
  }

  ## estimate skeleton and CPDAG of given data
  res <-
    if (is(pcObj, "pcAlgo"))
      pcObj
    else if(perfect)
      pcAlgo.Perfect(mcor, corMethod = "standard",directed=TRUE,u2pd=u2pd)
    else
      pcAlgo(dat, alpha = a, corMethod = "standard",directed=TRUE,u2pd=u2pd)

  ## prepare adjMatrix and skeleton {MM FIXME : can be improved}
  amat <- wgtMatrix(res@graph)
  amat[which(amat!=0)] <- 1 ## i->j if amat[j,i]==1
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel!=0] <- 1

  if (method=="local") {
##############################
    ## local method
    ## Main Input: mcov
##############################
    ## find unique parents of x
    wgt.est <- wgtMatrix(res@graph)
    tmp <- wgt.est-t(wgt.est)
    tmp[which(tmp<0)] <- 0
    wgt.unique <- tmp
    pa1 <- which(wgt.unique[x.pos,]!=0)
    if (y.pos %in% pa1) {
      ## x is parent of y -> zero effect
      beta.hat <- 0
    } else { ## y.pos not in pa1
      ## find ambiguous parents of x
      wgt.ambig <- wgt.est-wgt.unique
      pa2 <- which(wgt.ambig[x.pos,]!=0)
      if (verbose==2) {
        cat("\n\nx=",x.pos,"y=",y.pos,"\n")
        cat("pa1=",pa1,"\n")
        cat("pa2=",pa2,"\n")
      }

      ## estimate beta
      if (length(pa2)==0) {
        beta.hat <- lm.cov(mcov,y.pos,c(x.pos,pa1))
        if (verbose==2) {cat("Fit - y:",y.pos,"x:",c(x.pos,pa1),
              "|b.hat=",beta.hat)}
      } else {
        beta.hat <- NA
        ii <- 1
        ## no member of pa2
        pa2.f <- pa2
        pa2.t <- NA
        ## check for new collider
        if (collTest) {
          tmpColl <- check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)
        }
        if (!tmpColl | !collTest) {
          beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
          if (verbose==2) {cat("\ny:",y.pos,"x:",c(x.pos,pa1),"|b.hat=",
                beta.hat[ii])}
        } else {
          ## cat("\nx:",x.pos," pa1:",pa1," pa2.t:",pa2.t," pa2.f:",pa2.f)
        }
        ## exactly one member of pa2
        for (i2 in seq_along(pa2)) {
          ## check for new collider
          pa2.f <- pa2[-i2]
          pa2.t <- pa2[i2]
          if (collTest) {
            tmpColl <- check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)
          }
          if (!tmpColl | !collTest) {
            ii <-  ii+1
            if (y.pos %in% pa2.t) {
              ## cat("Y in Parents: ",y.pos," in ",pa2.t,"\n")
              beta.hat[ii] <- 0
            } else {
              beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1,pa2[i2]))
            }
            if (verbose==2) {cat("\ny:",y.pos,"x:",c(x.pos,pa1,pa2[i2]),
                  "|b.hat=",beta.hat[ii])}
          } else {
            ## cat("\nx:",x.pos," pa1:",pa1," pa2.t:",pa2.t," pa2.f:",pa2.f)
          }
        }
        ## higher order subsets
        if (length(pa2)>1) {
          for (i in 2:length(pa2)) {
            pa.tmp <- combn(pa2,i,simplify=TRUE)
            n.comb <- ncol(pa.tmp)
            for (j in 1:n.comb) {
              pa2.f <- setdiff(pa2,pa.tmp[,j])
              pa2.t <- pa.tmp[,j]
              ## teste auf neuen collider
              if (collTest) {
                tmpColl <- check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)
              }
              if (!tmpColl | !collTest) {
                ii <- ii+1
                if (y.pos %in% pa2.t) {
                  cat("Y in Parents: ",y.pos," in ",pa2.t,"\n")
                  beta.hat[ii] <- 0
                } else {
                  beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1,pa.tmp[,j]))
                }
                if (verbose==2) {cat("\ny:",y.pos,"x:",c(x.pos,pa1,pa.tmp[,j]),
                      "|b.hat=",beta.hat[ii])}
              } else {
                ## cat("\nx:",x.pos," pa1:",pa1," pa2.t:",pa2.t," pa2.f:",pa2.f)
              }
            }
          }
        }
      } ## if pa2
    } ## if y in pa1
  } else {
##############################
    ## global method
    ## Main Input: mcov
##############################
    p <- numNodes(res@graph)
    am.pdag <- wgtMatrix(res@graph)
    am.pdag[am.pdag!=0] <- 1
    ## find all DAGs if not provided externally
    if (is.na(all.dags)) {
      ad <- allDags(am.pdag,am.pdag,NULL)
    } else {
      ad <- all.dags
    }
    n.dags <- nrow(ad)
    beta.hat <- rep(NA,n.dags)
    if (n.dags>0) {
      if (myplot) {
        ## x11()
        par(mfrow=c(ceiling(sqrt(n.dags)), round(sqrt(n.dags)) ))
      }
      for (i in 1:n.dags) {
        ## compute effect for every DAG
        gDag <- as(matrix(ad[i,],p,p),"graphNEL")
        if (myplot) plot(gDag)
        ## path from y to x
        rev.pth <- RBGL::sp.between(gDag,as.character(y.pos),
                                    as.character(x.pos))[[1]]$path
        if (length(rev.pth)>1) {
          ## if reverse path exists, beta=0
          beta.hat[i] <- 0
        } else {
          ## path from x to y
          pth <- RBGL::sp.between(gDag,as.character(x.pos),
                                  as.character(y.pos))[[1]]$path
          if (length(pth)<2) {
            ## sic! There is NO path from x to y
            beta.hat[i] <- 0
          } else {
            ## There is a path from x to y
            wgt.unique <- t(matrix(ad[i,],p,p)) ## wgt.est is wgtMatrix of DAG
            pa1 <- which(wgt.unique[x.pos,]!=0)
            if (y.pos %in% pa1) {
              cat("Y in Parents: ",y.pos," in ",pa1,"\n")
              beta.hat[i] <- 0
            } else {
              beta.hat[i] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
            }
            if (verbose==2)
              cat("Fit - y:",y.pos,"x:",c(x.pos,pa1), "|b.hat=",beta.hat,"\n")
          } ## if length(pth)
        } ## if rev.pth
      } ## for n.dags
    } ## if n.dags
  } ## if method
  beta.hat
}


## DEPRECATED! -- use  ida() / idafast() --
beta.special.pcObj <- function(x.pos,y.pos,pcObj,mcov=NA,amat=NA,amatSkel=NA,
                               t.amat=NA)
{
  ## Purpose: Estimate the causal effect of x on y; the pcObj has to be
  ## precomputed. This method is intended to be a fast version of
  ##
  ## beta.special(dat=NA,x.pos,y.pos,verbose=0,a=NA,myDAG=NA,myplot=FALSE,
  ## perfect=FALSE,method="local",collTest=TRUE,pcObj=pcObj,all.dags=NA,u2pd="relaxed")
  ##
  ## Thus, this is a faster version for the local method given a
  ## precomputed PC-Algo Object (relaxed udag2pdag, so CPDAG might not
  ## be a real CPDAG; this does not matter, since we try not to extend).
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x.pos, y.pos: Column of x and y in d.mat
  ## - pcObj: Fit of pc Algorithm (semidirected); if this is available, no
  ##   new fit is done
  ## - mcov: covariance matrix of pcObj fit
  ## - amat,amatSkel,g2,t.amat are variants of the adjacency matrix that
  ##   are used internally but can be precomputed; the relevant code
  ##   is commented out
  ## ----------------------------------------------------------------------
  ## Value: List with two elements
  ## - beta.res: beta.causal values
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 21 Nov 2007, 11:18

  cat("This function is deprecated and is only kept for backward compatibility.
Please use ida or idaFast instead\n")

  if (is.na(amat) | is.na(amatSkel) | is.na(t.amat)) {
    ## Code for computing precomputable variables
    ## prepare adjMatrix and skeleton {MM FIXME : can be improved}
    amat <- wgtMatrix(pcObj@graph)
    amat[which(amat!=0)] <- 1 ## i->j if amat[j,i]==1
    t.amat <- t(amat)
    amatSkel <- amat + t.amat
    amatSkel[amatSkel!=0] <- 1
  }

  ## find unique parents of x
  tmp <- amat-t.amat
  tmp[which(tmp<0)] <- 0
  wgt.unique <- tmp
  pa1 <- which(wgt.unique[x.pos,]!=0)
  if (y.pos %in% pa1) {
    cat("Y in Parents: ",y.pos," in ",pa1,"\n")
    beta.hat <- 0
  } else { ## y.pos not in pa1
    ## find ambiguous parents of x
    wgt.ambig <- amat-wgt.unique
    pa2 <- which(wgt.ambig[x.pos,]!=0)
    pa2 <- setdiff(pa2,y.pos)
    ## estimate beta
    if (length(pa2)==0) {
      beta.hat <- lm.cov(mcov,y.pos,c(x.pos,pa1))
    } else {
      beta.hat <- NA
      ii <- 1
      ## no member of pa2
      ## check for new collider
      pa2.f <- pa2
      pa2.t <- NA
      tmpColl <- check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)
      if (!tmpColl) {
        beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
      }
      ## exactly one member of pa2
      for (i2 in seq_along(pa2)) {
        ## check for new collider
        pa2.f <- pa2[-i2]
        pa2.t <- pa2[i2]
        tmpColl <- check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)
        if (!tmpColl) {
          ii <-  ii+1
          if (y.pos %in% pa2.t) {
            cat("Y in Parents: ",y.pos," in ",pa2.t,"\n")
            beta.hat[ii] <- 0
          } else {
            beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1,pa2[i2]))
          }
        }
      }
      ## higher order subsets
      if (length(pa2)>1) {
        for (i in 2:length(pa2)) {
          pa.tmp <- combn(pa2,i,simplify=TRUE)
          n.comb <- ncol(pa.tmp)
          for (j in 1:n.comb) {
            ## teste auf neuen collider
            pa2.f <- setdiff(pa2,pa.tmp[,j])
            pa2.t <- pa.tmp[,j]
            tmpColl <- check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)
            if (!tmpColl) {
              ii <- ii+1
              if (y.pos %in% pa2.t) {
                cat("Y in Parents: ",y.pos," in ",pa2.t,"\n")
                beta.hat[ii] <- 0
              } else {
                beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1,pa.tmp[,j]))
              }
            }
          }
        }
      } ## if pa2
    } ## length(pa2)
  } ## y.pos %in% pa2
  beta.hat
}

##- lm.cov <- function(C,y,x) {
##-   sig <- C[x,x]
##-   beta <- solve(sig)%*%C[x,y,drop=FALSE]
##-   beta[1]
##- }

##' @title
##' @param C covariance matrix
##' @param y column of response
##' @param x columns of expl. vars
##' @return
lm.cov <- function (C, y, x) {
  solve(C[x, x], C[x, y, drop = FALSE])[1, ]
}

causalEffect <- function(g,y,x) {
  ## Compute true causal effect of x on y in g
  wmat <- wgtMatrix(g)
  p <- ncol(wmat)
  vec <- matrix(0,p,1)
  vec[x] <- 1
  ## compute and return  beta_{true} :
  if(y-x > 1) {
    for (i in (x+1):y) vec[i] <- wmat[i,]%*%vec
    vec[y]
  } else {
    wmat[y,x]
  }
}

check.new.coll <- function(amat,amatSkel,x,pa1,pa2.t,pa2.f) {
  ## Check if undirected edges that are pointed to x create a new v-structure
  ## Additionally check, if edges that are pointed away from x create
  ## new v-structure; i.e. x -> pa <- papa would be problematic
  ## pa1 are definit parents of x
  ## pa2 are undirected "parents" of x
  ## pa2.t are the nodes in pa2 that are directed towards pa2
  ## pa2.f are the nodes in pa2 that are directed away from pa2
  ## Value is TRUE, if new collider is introduced
  res <- FALSE
  if ((length(pa2.t)>0) && any(!is.na(pa2.t))) {
    ## check whether all pa1 and all pa2.t are connected;
    ## if no, there is a new collider
    if ((length(pa1)>0) && any(!is.na(pa1))) {
      res <- (min(amatSkel[pa1,pa2.t])==0) ## TRUE if new collider
    }
    ## in addition, all pa2.t have to be connected
    if ((length(pa2.t)>1) && (!res)) {
      tmp <- amatSkel[pa2.t,pa2.t]
      diag(tmp) <- 1
      res2 <- (min(tmp)==0) ## TRUE if new collider
      res <- (res|res2)
    }
  }
  if (!res && ((length(pa2.f)>0) && any(!is.na(pa2.f)))) {
    ## consider here only the DIRECTED Parents of pa2.f
    ## remove undirected edges
    amatTmp <- amat
    amatTmp <- amatTmp-t(amatTmp)
    amatTmp[amatTmp<0] <- 0
    tmp <- amatTmp[pa2.f,,drop=FALSE]
    ## find parents of pa2.f
    papa <- setdiff(which(apply(tmp,2,sum)!=0),x)
    ## if any node in papa is not directly connected to x, there is a new
    ## collider
    if (length(papa) > 0)
      res <- res | (min(amatSkel[x,papa])==0) ## TRUE if new collider
  }
  res
}

allDags <- function(gm,a,tmp, verbose=FALSE)
{
  ## Purpose: Find all DAGs for a given PDAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: Adjacency matrix of initial PDAG; only 0-1 entries
  ##   i -> j iff gm(j,i)=1
  ## - a: copy of gm
  ## - tmp: NULL
  ## ----------------------------------------------------------------------
  ## Value:
  ## - one 0/1 adj.matrix per row
  ## Reversion to graph: as(matrix(res[i,],p,p),"graphNEL")
  ## Reversion to wgtMatrix (i->j iff a[j,i]=1): t(matrix(res[i,],p,p))
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  7 Apr 2008, 14:08
  if (sum(a) == 0) {
    if (verbose) {
      cat("Last Call - Final Graph: \n")
      print(gm)
      cat("#################### \n")
    }
    tmp2 <- rbind(tmp,c(t(gm)))
    if (all(!duplicated(tmp2))) tmp <- tmp2
  } else {
    sinks <- find.sink(a)
    if (verbose) {
      cat("Main Call: ################## \n")
        print(gm)
      print(a)
      cat("Sinks: ",sinks,"\n")
    }
    for(x in sinks) {
      if (verbose) cat("Try removing", x," in a.\n")
      gm2 <- gm
      a2 <- a
      if (adj.check(a,x)) {
        inc.to.x <- a[, x] == 1 & a[x, ] == 1
        if (any(inc.to.x)) {
          real.inc.to.x <- as.numeric(row.names(a)[inc.to.x])
          real.x <- as.numeric(row.names(a)[x])
          gm2[real.x, real.inc.to.x] <- 1
          gm2[real.inc.to.x, real.x] <- 0
        }
        a2 <- a[-x,-x]
        if (verbose) {
          cat("Removed sink",as.numeric(row.names(a)[x]),
              "in g (", x,"in a).\n")
          cat("New graphs: \n")
          print(gm2)
          print(a)
        }
        tmp <- allDags(gm2,a2,tmp, verbose)
      }
    }
  }
  tmp
}

pcAlgo.Perfect <- function(C, cutoff= 1e-8, corMethod = "standard", verbose = 0,
                           directed=FALSE, u2pd="rand", psepset=FALSE) {
  ## Purpose: Perform PC-Algorithm, i.e., estimate skeleton of DAG given data
  ## Output is an unoriented graph object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - C: True Correlation matrix
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - verbose: 0-no output, 1-small output, 2-details
  ## - u2pd: Function for converting udag to pdag
  ##   "rand": udag2pdag
  ##   "relaxed": udag2pdagRelaxed
  ##   "retry": udag2pdagSpecial
  ## - psepset: Also check possible sep sets.
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler
  ## Modification: Diego Colombo, Sept 2009
  ## backward compatibility
  stopifnot((p <- nrow(C)) >= 2)
  if (verbose==FALSE) verbose <- 0
  if (verbose==TRUE) verbose <- 1
  cl <- match.call()
  seq_p <- 1:p
  pcMin <- matrix(Inf, p,p)
  diag(pcMin) <- 0
  sepset <- pl <- vector("list",p)
  for (i in seq_p) sepset[[i]] <- pl
  n.edgetests <- numeric(1)# final length = max { ord}
  ## G := complete graph :
  G <- matrix(TRUE, p,p)
  diag(G) <- FALSE

  done <- FALSE
  ord <- 0
  while (!done && any(G)) {
    n.edgetests[ord+1] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)
    ## For comparison with C++ sort according to first row
    ind <- ind[order(ind[,1]) ,]
    remainingEdgeTests <- nrow(ind)
    if(verbose>=1)
      cat("Order=",ord,"; remaining edges:",remainingEdgeTests,"\n", sep='')
    for (i in 1:remainingEdgeTests) {
      if(verbose>=1) { if(i%%100==0) cat("|i=",i,"|iMax=",nrow(ind),"\n") }
      x <- ind[i,1]
      y <- ind[i,2]
      ##      done <- !G[y,x] # i.e. (x,y) was not already deleted in its (y,x) "version"
      ##      if(!done) {
      if (G[y,x]) {
        nbrsBool <- G[,x]
        nbrsBool[y] <- FALSE
        nbrs <- seq_p[nbrsBool]
        length_nbrs <- length(nbrs)
        ## if(verbose)

        ##        done <- length_nbrs < ord
        ##        if (!done) {
        if (length_nbrs >= ord) {
          if (length_nbrs > ord) done <- FALSE
          S <- seq_len(ord)

          ## now includes special cases (ord == 0) or (length_nbrs == 1):
          repeat { ## condition w.r.to all  nbrs[S] of size 'ord' :
            ##  if (condIndFisherZ(x,y, nbrs[S], C,n, cutoff,verbose)) {
            ## MM: want to use this: --- but it changes the result in some cases!!
            ##            cat("X=",x,"|Y=",y,"|ord=",ord,"|nbrs=",nbrs[S],"|iMax=",nrow(ind),"\n")
            n.edgetests[ord+1] <- n.edgetests[ord+1]+1
            pc.val <- pcorOrder(x,y,nbrs[S],C)
            if (abs(pc.val)<pcMin[x,y]) pcMin[x,y] <- abs(pc.val)
            if (verbose==2) cat(paste("x:",x,"y:",y,"S:"),nbrs[S],paste("pc:",pc.val,"\n"))
            if (abs(pc.val) <= cutoff) {
              ##              ##  pnorm(abs(z), lower.tail = FALSE) is the P-value
              G[x,y] <- G[y,x] <- FALSE
              sepset[[x]][[y]] <- nbrs[S]
              break
            }
            else {
              nextSet <- getNextSet(length_nbrs, ord, S)
              if(nextSet$wasLast)
                break
              S <- nextSet$nextSet
            }
          }

        } } ## end if(!done)

    } ## end for(i ..)
    ord <- ord+1
    ##    n.edgetests[ord] <- remainingEdgeTests
  }

  if (psepset) {
    amat <- G
    amat[amat==TRUE] <- 1
    amat[amat==FALSE] <- 0
    ind <- which(amat==1, arr.ind=TRUE)
    ## Orient colliders
    for (i in 1:dim(ind)[1]) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(amat[y,]==1),x) ## x-y-z

      if (length(allZ)>0) {
        for (j in seq_along(allZ)) {
          z <- allZ[j]
          if ((amat[x,z]==0) && !((y %in% sepset[[x]][[z]]) |(y %in% sepset[[z]][[x]]))) {
            if (verbose == 2) {
              cat("\n",x,"*->",y,"<-*",z,"\n")
              cat("Sxz=",sepset[[z]][[x]],"and","Szx=",sepset[[x]][[z]],"\n")
            }

            ## x o-> y <-o z
            amat[x,y] <- amat[z,y] <- 2

          } ## if
        } ## for
      } ## if
    } ## for

    ## Compute poss. sepsets
    for (x in 1:p) {
      attr(x,'class') <- 'possibledsep'
      if (any(amat[x,]!=0)){
        tf1 <- setdiff(reach(x,-1,-1,amat),x)
        for (y in seq_p[amat[x,]!=0]) {
          ## tf = possible_d_sep(amat,x,y)
          tf <- setdiff(tf1,y)
          ## test
          if (length(tf)>0) {
            pc.val <- pcorOrder(x,y,tf,C)
            if (abs(pc.val)<pcMin[x,y]){
              pcMin[x,y] <- abs(pc.val)
            }
            if (abs(pc.val) <= cutoff) {
              ## delete x-y
              amat[x,y] <- amat[y,x] <- 0
              ## save pos d-sepset in sepset
              sepset[[x]][[y]] <- tf
              if (verbose==2){
                cat("Delete edge",x,"-",y,"\n")
              }
            }
            if (verbose == 2) {
              cat("Possible-D-Sep of", x, "and", y, "is", tf, " - pc = ",pc.val,"\n")
            }
          }
        }
      }
    }
    G[amat==0] <- FALSE
    G[amat==1] <- TRUE
  } ## end if(psepset)

  if(verbose>=1) { cat("Final graph adjacency matrix:\n"); print(symnum(G)) }

  ## transform matrix to graph object :
  if (sum(G) == 0) {
    Gobject <- new("graphNEL", nodes = as.character(seq_p))
  }
  else {
    colnames(G) <- rownames(G) <- as.character(seq_p)
    Gobject <- as(G,"graphNEL")
  }

  for (i in 1:(p-1)) {
    for (j in 2:p) {
      pcMin[i,j] <- pcMin[j,i] <- min(pcMin[i,j],pcMin[j,i])
    }
  }

  res <- new("pcAlgo",
             graph = Gobject,
             call = cl, n = as.integer(1), max.ord = as.integer(ord-1),
             n.edgetests = n.edgetests, sepset = sepset,
             zMin = pcMin)

  if (directed)
    switch(u2pd,
           "rand" = udag2pdag(res),
           "retry" = udag2pdagSpecial(res)$pcObj,
           "relaxed" = udag2pdagRelaxed(res))
  else
    res
}## pcAlgo.Perfect

### reach(): currently only called from  pcAlgo() and pcAlgo.Perfect()
### -------  and only in "possibledsep" version
## Function that computes the Possible d-sepset, done by Spirtes
reach <- function(a,b,c,adjacency)
{
  ## reachable      set of vertices;
  ## edgeslist      array[1..maxvertex] of list of edges
  ## numvertex      integer
  ## labeled        array (by depth) of list of edges that have been labeled

  makeedge <- function(x,y) list(list(x,y))

  legal.pdsep <- function(r,s) {
    ## Modifying global 'edgeslist'
    if ((adjacency[r[[1]],r[[2]]] == 2 &&
         adjacency[s,     r[[2]]] == 2 && r[[1]] != s) ||
        (adjacency[r[[1]],s] != 0 && r[[1]] != s)) {
      edgeslist[[r[[2]]]] <<- setdiff(edgeslist[[r[[2]]]],s)
      makeedge(r[[2]],s)
    }
  }

  initialize.pdsep <- function(x,y) mapply(makeedge, x=x, y=y)

  labeled <- list()
  numvertex <- dim(adjacency)[1]
  edgeslist <- list()
  for (i in 1:numvertex)
    edgeslist <- c(edgeslist,list(which(adjacency[,i] != 0)))
  labeled[[1]] <- initialize.pdsep(a, edgeslist[[a]])
  edgeslist[[a]] <- list()
  depth <- 2
  repeat {
    labeled[[depth]] <- list()
    for (i in seq_along(labeled[[depth-1]])) {
      lab.i <- labeled[[depth-1]][[i]]
      edgestemp <- edgeslist[[lab.i[[2]]]]
      if (length(edgestemp) == 0) break
      for (j in seq_along(edgestemp))
        labeled[[depth]] <- union(legal.pdsep(lab.i, edgestemp[[j]]),
                                  labeled[[depth]])
    }
    if (length(labeled[[depth]]) == 0)
      break
    ## else :
    depth <- depth  + 1
  }
  unique(unlist(labeled))
}

udag2pag <- function(gInput, rules=rep(TRUE,10), verbose=TRUE, unfVect=NULL)
{
  ## Purpose:Transform the Skeleton of a pcAlgo-object to a PAG using
  ## the rules of Zhang. The output is an adjacency matrix.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gInput: pcAlgo object
  ## - rules: array of length 10 wich contains TRUE or FALSE corrsponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - verbose: 0 - no output, 1 - detailed output
  ## - unfVect: Vector with unfaithful triples (coded as number using triple2numb)
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 6 Mar 2009; cleanup: Martin Maechler, 2010

  stopifnot(is.logical(rules), length(rules) == 10)

  if (numEdges(gInput@graph) == 0)
    return(matrix(NA_integer_, 0,0))

  ## Notation:
  ## ----------------------------------------------------------------------
  ## 0: no edge
  ## 1: -o
  ## 2: -> (arrowhead)
  ## 3: - (tail)
  ## a=alpha
  ## b=beta
  ## c=gamma
  ## d=theta
  g <- as(gInput@graph, "matrix")
  p <- dim(g)[1]
  evisit <- amat0 <- amat1 <- matrix(0, p, p)
  pag <- g
  ind <- which(g == 1, arr.ind = TRUE)
  for (i in seq_len(nrow(ind))) {
    x <- ind[i, 1]
    y <- ind[i, 2]
    allZ <- setdiff(which(g[y, ] == 1), x)
    for (z in allZ) {
      if (g[x, z] == 0 &&
          !((y %in% gInput@sepset[[x]][[z]]) ||
            (y %in% gInput@sepset[[z]][[x]]))) {
        ## normal version
        if (length(unfVect)==0) {
          if (verbose) {
            cat("\n", x, "*->", y, "<-*", z, "\n")
            cat("Sxz=", gInput@sepset[[z]][[x]], "and",
                "Szx=", gInput@sepset[[x]][[z]], "\n")
          }
          pag[x, y] <- pag[z, y] <- 2
        }
        ## conservative version
        else {
          ## check if x-y-z is faithful
          if (!any(unfVect == triple2numb(p,x,y,z)) &&
              !any(unfVect == triple2numb(p,z,y,x))) {
            if (verbose) {
              cat("\n", x, "*->", y, "<-*", z, "\n")
              cat("Sxz=", gInput@sepset[[z]][[x]], "and",
                  "Szx=", gInput@sepset[[x]][[z]], "\n")
            }
            pag[x, y] <- pag[z, y] <- 2
          }
        }
      }
    }
  }
  old_pag1 <- matrix(0, p, p)

  while (any(old_pag1 != pag)) {
    old_pag1 <- pag
    if (rules[1]) {
      ind <- which((pag == 2 & t(pag) != 0), arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        a <- ind[i, 1]
        b <- ind[i, 2]
        indC <- which((pag[b, ] != 0 & pag[, b] == 1) &
                      (pag[a, ] == 0 & pag[, a] == 0))
        indC <- setdiff(indC, a)
        if (length(indC) > 0) {
          if (length(unfVect)==0) {
            pag[b, indC] <- 2
            pag[indC, b] <- 3
            if (verbose)
              cat("\nRule 1:", a, "*->", b, "o-*", indC,
                  " where ", a, " and ", indC, " not connected: ", b, "->", indC, "\n")
          }
          else {
            for (c in indC) {
              ## check that a-b-c faithful
              if (!any(unfVect == triple2numb(p,a,b,c)) &&
                  !any(unfVect == triple2numb(p,c,b,a))) {
                pag[b, c] <- 2
                pag[c, b] <- 3
                if(verbose)
                  cat("\nRule 1':", a, "*->", b, "o-*", c,
                      " where ", a, " and ", c, " not connected: ",a, b, c,
                      "faithful triple", b, "->", c, "\n")
              }
            }
          }
        }
      }
    } ## r..[1]

    if (rules[2]) {
      ind <- which((pag == 1 & t(pag) != 0), arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        a <- ind[i, 1]
        c <- ind[i, 2]
        indB <- which((pag[a, ] == 2 & pag[, a] == 3 &
                       pag[c, ] != 0 & pag[, c] == 2)
                      |
                      (pag[a, ] == 2 & pag[, a] != 0 &
                       pag[c, ] == 3 & pag[, c] == 2))
        if (length(indB) > 0) {
          pag[a, c] <- 2
          if (verbose) {
            cat("\nRule 2:", a, "->", indB, "*->",
                c, "or", a, "*->", indB, "->", c, "and",
                a, "*-o", c, ":", a, "*->", c, "\n")
          }
        }
      }
    } ## r..[2]

    if (rules[3]) {
      ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        b <- ind[i, 1]
        d <- ind[i, 2]
        indAC <- which((pag[b, ] != 0 & pag[, b] == 2) &
                       (pag[, d] == 1 & pag[d, ] != 0))
        if (length(indAC) >= 2) {
          ## normal version
          if (length(unfVect)==0) {
            counter <- 0
            while ((counter < (length(indAC) - 1)) &
                   (pag[d, b] != 2)) {
              counter <- counter + 1
              ii <- counter
              while (ii < length(indAC) && pag[d, b] != 2) {
                ii <- ii + 1
                if (pag[indAC[counter], indAC[ii]] == 0 &&
                    pag[indAC[ii], indAC[counter]] == 0) {
                  if(verbose) cat("\nRule 3:", d, "*->", b, "\n")
                  pag[d, b] <- 2
                }
              }
            }
          }
          ## conservative version
          else {
            comb.indAC <- combn(indAC, 2)
            for (j in 1:dim(comb.indAC)[2]) {
              a <- comb.indAC[1,j]
              c <- comb.indAC[2,j]
              if (pag[a,c]==0 && pag[c,a]==0 && c!=a) {
                ## check fatihfulness a-d-c
                if (!any(unfVect == triple2numb(p,a,d,c)) &&
                    !any(unfVect == triple2numb(p,c,d,a))) {
                  pag[d, b] <- 2
                  if(verbose)
                    cat("\nRule 3':", a, d, c, "faithful",  d, "*->", b, "\n")
                }
              }
            } ## for( j )
          }
        }
      }
    } ## r..[3]

    if (rules[4]) {
      ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        b <- ind[i, 1]
        c <- ind[i, 2]
        indA <- which((pag[b, ] == 2 & pag[, b] != 0) &
                      (pag[c, ] == 3 & pag[, c] == 2))
        for(a in indA)
          pag <- discr.path(path = c(c, b, a), pag = pag,
                            gInput = gInput, verbose = verbose)
      }
    } ## r..[4]

    if (rules[5]) {
      ind <- which((pag == 1 & t(pag) == 1), arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        a <- ind[i, 1]
        b <- ind[i, 2]
        indC <- which((pag[a, ] == 1 & pag[, a] == 1) &
                      (pag[b, ] == 0 & pag[, b] == 0))
        indC <- setdiff(indC, b)
        indD <- which((pag[b, ] == 1 & pag[, b] == 1) &
                      (pag[a, ] == 0 & pag[, a] == 0))
        indD <- setdiff(indD, a)
        if (length(indC) > 0 && length(indD) > 0) {
          for (c in indC) {
            for (d in indD) {
              if (pag[c, d] == 1 && pag[d, c] == 1) {
                if (length(unfVect)==0) {
                  pag[a, b] <- pag[b, a] <- 3
                  pag[a, c] <- pag[c, a] <- 3
                  pag[c, d] <- pag[d, c] <- 3
                  pag[d, b] <- pag[b, d] <- 3
                  if(verbose)
                    cat("\nRule 5: There exists an uncovered circle path between",
                        a, "and", b, ":", a, "-", b,
                        "and", a, "-", c, "-", d, "-", b, "\n")
                }
                else {
                  ## check that every triple on the circle is faithful
                  path2check <- c(a,c,d,b)
                  if (faith.check(path2check, unfVect, p)) {
                    pag[a, b] <- pag[b, a] <- 3
                    pag[a, c] <- pag[c, a] <- 3
                    pag[c, d] <- pag[d, c] <- 3
                    pag[d, b] <- pag[b, d] <- 3
                    if(verbose)
                      cat("\nRule 5': There exists a faithful uncovered circle path between",
                          a, "and", b, ":", a, "-", b,
                          "and", a, "-", c, "-", d, "-", b, "\n")
                  }
                }
              }
              else {
                path <- c(a, c, d, b)
                pag <- ucp(path = path, pag = pag,
                           verbose = verbose, unfVect=unfVect, p=p)
              }
            } ## for d
          }## for c
        }
      }
    } ## r..[5]

    if (rules[6]) {
      ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        b <- ind[i, 1]
        c <- ind[i, 2]
        indA <- which(pag[b, ] == 3 & pag[, b] == 3)
        if (length(indA) > 0) {
          pag[c, b] <- 3
          if (verbose)
            cat("\nRule 6:", a, "-", b, "o-*", c, ":", b, "-*", c, "\n")
        }
      }
    } ## r..[6]

    if (rules[7]) {
      ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        b <- ind[i, 1]
        c <- ind[i, 2]
        indA <- which((pag[b, ] == 3 & pag[, b] ==
                       1) & (pag[c, ] == 0 & pag[, c] == 0))
        indA <- setdiff(indA, c)
        if (length(indA) > 0) {
          if (length(unfVect)==0) {
            pag[c, b] <- 3
            if(verbose)
              cat("\nRule 7:", indA, "-o", b, "o-*", c, "and", indA, " and", c,
                  " are not adjacent:", b, "-*", c, "\n")
          }
          else {
            for (a in indA) {
              ## check fatihfulness of a-b-c
              if (!any(unfVect == triple2numb(p,a,b,c)) &&
                  !any(unfVect == triple2numb(p,c,b,a))) {
                pag[c, b] <- 3
                if(verbose)
                  cat("\nRule 7':", a, "-o", b, "o-*", c, "and", a, " and", c,
                      " are not adjacent and", a,b,c, "faithful:", b, "-*", c, "\n")
              }
            }
          }
        }
      }
    } ## r..[7]

    if (rules[8]) {
      ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        a <- ind[i, 1]
        c <- ind[i, 2]
        indB <- which(((pag[a, ] == 2 & pag[, a] == 3) |
                       (pag[a, ] == 1 & pag[, a] == 3)) &
                      (pag[c, ] == 3 & pag[, c] == 2))
        if (length(indB) > 0) {
          pag[c, a] <- 3
          if(verbose)
            cat("\nRule 8:", a, "->", indB, "->", c,
                "or", a, "-o", indB, "->", c, "and",
                a, "o->", c, ":", a, "->", c, "\n")
        }
      }
    } ## r..[8]

    if (rules[9]) {
      ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        a <- ind[i, 1]
        c <- ind[i, 2]
        indB <- which((pag[a, ] == 2 | pag[a, ] == 1) &
                      (pag[, a] == 1 | pag[, a] == 3) &
                      (pag[c, ] == 0 & pag[, c] == 0))
        indB <- setdiff(indB, c)
        for (b in indB) {
          indD <- which((pag[b, ] == 2 | pag[b, ] == 1) &
                        (pag[, b] == 1 | pag[, b] == 3) &
                        (pag[a, ] == 0 & pag[, a] == 0))
          indD <- setdiff(indD, a)
          for (d in indD) {
            if (pag[c, a] != 3) {
              if ((pag[c, d] == 1 | pag[c, d] == 3) &
                  (pag[d, c] == 1 | pag[d, c] == 2)) {
                ## normal version
                if (length(unfVect)==0) {
                  pag[c, a] <- 3
                  if(verbose)
                    cat("\nRule 9: There exists an upd between",
                        a, "and", c, ":", a, " ->", c, "\n")
                }
                ## conservative version
                else {
                  path2check <- c(a,b,d,c)
                  if (faith.check(path2check, unfVect, p)) {
                    pag[c, a] <- 3
                    if(verbose)
                      cat("\nRule 9': There exists a faithful upd between",
                          a, "and", c, ":", a, " ->", c, "\n")
                  }
                }
              }
              else {
                path <- c(c, a, b, d)
                pag <- upd(path=path, pag=pag, verbose=verbose, unfVect=unfVect, p=p)
              }
            }
          } ## for( d )
        }   ## for( b )
      }
    } ## r..[9]

    if (rules[10]) {
      ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        a <- ind[i, 1]
        c <- ind[i, 2]
        indB <- which(pag[c, ] == 3 & pag[, c] == 2)
        for (b in indB) {
          indD <- setdiff(indB, b)
          if (length(indD) > 0 && pag[c, a] != 3) {
            for (d in indD) {
              if ((pag[a, b] == 1 || pag[a, b] == 2) &&
                  (pag[b, a] == 1 || pag[b, a] == 3) &&
                  (pag[a, d] == 1 || pag[a, d] == 2) &&
                  (pag[d, a] == 1 || pag[d, a] == 3) &&
                  pag[d, b] == 0 && pag[b, d] == 0) {
                ## normal version
                if (length(unfVect)==0) {
                  pag[c, a] <- 3
                  if(verbose)
                    cat("\nRule 10 with mu = beta = ", b,
                        "and omega = theta =", d, ":", a, "->", c, "\n")
                }
                ## conservative version
                else {
                  ## check faithfulness of b-a-d
                  if (!any(unfVect == triple2numb(p,b,a,d)) &
                      !any(unfVect == triple2numb(p,d,a,b))) {
                    pag[c, a] <- 3
                    if(verbose)
                      cat("\nRule 10' with mu = beta = ",
                          b, "and omega = theta =", d,
                          "and",b,a,d,"faithful:", a, "->", c, "\n")
                  }
                }
              }
              else {
                indA <- which((pag[a, ] == 1 | pag[a, ] == 2) &
                              (pag[, a] == 1 | pag[, a] == 3), arr.ind = TRUE)
                indA <- setdiff(indA, c)
                for (first.pos in indA) {
                  indAA <- setdiff(indA, first.pos)
                  if ((length(indAA) > 0) && (pag[c, a] != 3)) {
                    for (s in seq_along(indAA)) {
                      sec.pos <- indAA[s]
                      p1 <- find.upd(path = c(first.pos, b), a = a, pag = pag,
                                     verbose = verbose, unfVect=unfVect, p=p)
                      p2 <- find.upd(path = c(sec.pos, d), a = a, pag = pag,
                                     verbose = verbose, unfVect=unfVect, p=p)
                      if (p1$res && p2$res) {
                        mu <- p1$uncov.path[1]
                        omega <- p2$uncov.path[1]
                        if (mu != omega && pag[mu, omega] == 0 && pag[omega, mu] == 0) {
                          ## normal version
                          if (length(unfVect)==0) {
                            pag[c, a] <- 3
                            if(verbose) cat("\nRule 10:", a, "->", c, "\n")
                          }
                          ## conservative version
                          else {
                            if (!any(unfVect == triple2numb(p,mu,a,omega)) &&
                                !any(unfVect == triple2numb(p,omega,a,mu))) {
                              pag[c, a] <- 3
                              if (verbose)
                                cat("\nRule 10':",mu,a,omega,"faithful:", a,
                                    "->", c, "\n")
                            }
                          }
                        }
                      }
                    } ## for s
                  }
                } ## for first.pos
              }
            } ## for d
          }
        } ## for b
      }
    } ## r..[10]

  } ## {while}

  return(pag)

} ## udag2



plotAG <- function(amat)
{
  ## Purpose: Plot ancestral graph
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: Adjacency matrix
  ##   amat[i,j]=3 & amat[j,i]=1 iff i 1-3 j
  ##   "0": no edge; "1": circle; "2": arrow; "3": tail
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 16 Feb 2009, 18:01
  check.Rgraphviz()

  g <- as(amat,"graphNEL")
  nn <- nodes(g)
  p <- numNodes(g)
  n.edges <- numEdges(g)
  ah.list <- at.list <- rep("none",n.edges)
  counter <- 0
  list.names <- NULL
  amat[amat==1] <- "odot"
  amat[amat==2] <- "normal"
  amat[amat==3] <- "none"
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      x <- nn[i]
      y <- nn[j]
      if (amat[x,y]!=0) {
        counter <- counter + 1
        ah.list[[counter]] <- amat[x,y]
        at.list[[counter]] <- amat[y,x]
        list.names <- c(list.names,paste(x,"~",y,sep=""))
      }
    }
  }
  names(ah.list) <- names(at.list) <- list.names

  edgeRenderInfo(g) <- list(arrowhead=ah.list, arrowtail=at.list)
  Rgraphviz::renderGraph(Rgraphviz::layoutGraph(g))
}

skeleton <- function(suffStat, indepTest, p, alpha, verbose = FALSE,
                     fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                     m.max = Inf) {

  ## Purpose: Perform undirected part of PC-Algorithm, i.e.,
  ## estimate skeleton of DAG given data
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - C: correlation matrix (only for continuous)
  ## - n: sample size
  ## - p: number of variables !! NEU
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - fixedGaps: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical); gaps fixed here are not changed
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## - datatype: distinguish between discrete and continuous data
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - gTrue: Graph object of true DAG
  ## ----------------------------------------------------------------------
  ## Value:
  ## - G, sepset, pMax, ord, n.edgetests
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 09.12.2009

  ## x,y,S konstruieren
  ##-   tst <- try(indepTest(x,y,S, obj))
  ##-   if(inherits(tst, "try-error"))
  ##-     stop("the 'indepTest' function does not work correctly with 'obj'")

  stopifnot((p <- as.integer(p)) >= 2)
  cl <- match.call()
  ## start skeleton

  ## fixed gaps
  if (is.null(fixedGaps)) {
    ## G := complete graph :
    G <- matrix(TRUE, p,p)
    diag(G) <- FALSE
  } else {
    if (!identical(dim(fixedGaps),c(p,p))) {
      stop("Dimensions of the dataset and fixedGaps do not agree.")
    } else {
      if (!all(fixedGaps == t(fixedGaps)))
        stop("fixedGaps must be symmetric")
      G <- !fixedGaps
    }
  } ## if(is.null(G))

  ## fixed edges
  if (is.null(fixedEdges)) {
    fixedEdges <- matrix(FALSE, p,p)
  } else {
    if (!(identical(dim(fixedEdges),c(p,p))))
      stop("Dimensions of the dataset and fixedEdges do not agree.")
    if (fixedEdges != t(fixedEdges))
      stop("fixedEdges must be symmetric")
  }

  seq_p <- seq_len(p)
  sepset <- pl <- vector("list",p)
  for (i in seq_p) sepset[[i]] <- pl
  ## save maximal p value
  pMax <- matrix(-Inf, p,p)
  diag(pMax) <- 1

  done <- FALSE
  ord <- 0L
  n.edgetests <- numeric(1)# final length = max { ord}

  while (!done && any(G) && ord <= m.max) {
    n.edgetests[ord1 <- ord+1L] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)
    ## For comparison with C++ sort according to first row
    ind <- ind[order(ind[,1]) ,]
    remainingEdgeTests <- nrow(ind)
    if(verbose)
      cat("Order=",ord,"; remaining edges:",remainingEdgeTests,"\n", sep='')
    for (i in 1:remainingEdgeTests) {
      if(verbose) { if(i%%100==0) cat("|i=",i,"|iMax=",nrow(ind),"\n") }
      x <- ind[i,1]
      y <- ind[i,2]
      if (G[y,x] && !fixedEdges[y,x]) {
        nbrsBool <- G[,x]
        nbrsBool[y] <- FALSE
        nbrs <- seq_p[nbrsBool]
        length_nbrs <- length(nbrs)
        if (length_nbrs >= ord) {
          if (length_nbrs > ord) done <- FALSE
          S <- seq_len(ord)
          repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
            n.edgetests[ord1] <- n.edgetests[ord1]+1
            pval <- indepTest(x,y, nbrs[S], suffStat)
            ## pval <- dsepTest(x,y,nbrs[S],gTrue,jp = jp)
            if (verbose) cat("x=",x," y=",y," S=",nbrs[S],": pval =",pval,"\n")
            if (is.na(pval)) pval <- if(NAdelete) 1 else 0
            if (pval > pMax[x,y]) pMax[x,y] <- pval
            if(pval >= alpha) { # independent
              G[x,y] <- G[y,x] <- FALSE
              sepset[[x]][[y]] <- nbrs[S]
              break
            } else {
              nextSet <- getNextSet(length_nbrs, ord, S)
              if(nextSet$wasLast)
                break
              S <- nextSet$nextSet
            }
          } ## repeat
        } ## if (length_nbrs >= ord)
      } ## if(!done)

    } ## for(i in 1:remainingEdgeTests)
    ord <- ord1
  } ## while

  for (i in 1:(p-1)) {
    for (j in 2:p) {
      pMax[i,j] <- pMax[j,i] <- max(pMax[i,j],pMax[j,i])
    } ## for (j in 2:p)
  } ## for (i in 1:(p-1))

  ## transform matrix to graph object :
  nnms <- as.character(seq_p)
  Gobject <-
    if (sum(G) == 0) {
      new("graphNEL", nodes = nnms)
    } else {
      colnames(G) <- rownames(G) <- nnms
      as(G,"graphNEL")
    }

  ## final object
  new("pcAlgo",
      graph = Gobject,
      call = cl, n = integer(0), max.ord = as.integer(ord-1),
      n.edgetests = n.edgetests, sepset = sepset,
      pMax = pMax, zMin = matrix(NA,1,1))

}## end{ skeleton }



pc <- function(suffStat, indepTest, p, alpha, verbose = FALSE, fixedGaps = NULL,
               fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,
               u2pd = "rand", conservative=FALSE) {

  ## Purpose: Perform PC-Algorithm, i.e., estimate skeleton of DAG given data
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - C: correlation matrix (only for continuous)
  ## - n: sample size
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - G: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical)
  ## - datatype: distinguish between discrete and continuous data
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - u2pd: Function for converting udag to pdag
  ##   "rand": udag2pdag
  ##   "relaxed": udag2pdagRelaxed
  ##   "retry": udag2pdagSpecial
  ## - gTrue: Graph suffStatect of true DAG
  ## - conservative: If TRUE, conservative PC is done
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler
  ## Modifications: Sarah Gerster, Date: July 2007
  ## Modifications: Diego Colombo, Date: Sept 2009
  ## Modifications: Markus Kalisch, Date: Dec 2009

  ## Initial Checks
  cl <- match.call()

  ## Skeleton
  skel <- skeleton(suffStat, indepTest, p, alpha, verbose = verbose,
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                   NAdelete = NAdelete, m.max = m.max)

  ## Orient edges
  if (!conservative) {
    switch (u2pd,
            "rand" = udag2pdag(skel),
            "retry" = udag2pdagSpecial(skel)$pcObj,
            "relaxed" = udag2pdagRelaxed(skel))
  }
  else {
    if (u2pd != "relaxed") stop("Conservative PC can only be run with 'u2pd = relaxed'")
        ## version.unf defined per default
    ## Tetrad CPC works with version.unf=c(2,1)
    ## see comment on pc.cons.intern for description of version.unf
    tmp <- pc.cons.intern(skel, suffStat, indepTest, alpha, verbose=verbose, version.unf=c(2,2))
    tripleList <- tmp$unfTripl
    ## vers <- tmp$vers
    udag2pdagRelaxed(gInput=skel, verbose=verbose, unfVect=tripleList)
  }
}


fci <- function(suffStat, indepTest, p, alpha, verbose = FALSE,
                fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                m.max = Inf, rules = rep(TRUE, 10), doPdsep = TRUE,
                conservative=c(FALSE,FALSE), biCC=FALSE, cons.rules=FALSE,
                labels = NA)
{
  ## Purpose: Perform PC-Algorithm, i.e., estimate skeleton of DAG given data
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - C: correlation matrix (only for continuous)
  ## - n: sample size
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - G: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical)
  ## - datatype: distinguish between discrete and continuous data
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - doPdsep: compute possible dsep
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Dec 2009

##################################################
  ## Initial Checks
##################################################
  if (all(!is.na(labels))) {
    stopifnot(length(labels) == p)
  } else {
    labels <- as.character(1:p)
  }
  cl <- match.call()
  if (verbose)
    cat("Compute Skeleton\n================\n")
  skel <- skeleton(suffStat, indepTest, p, alpha, verbose = verbose,
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges, NAdelete = NAdelete,
                   m.max = m.max)
  G <- (as(skel@graph, "matrix") != 0)
  sepset <- skel@sepset
  pMax <- skel@pMax
  n.edgetestsSKEL <- skel@n.edgetests
  max.ordSKEL <- skel@max.ord

##################################################
  ## Possible D-Sep
##################################################
  allPdsep <- NA
  tripleList <- NULL
  ## conservative version
  if (conservative[1]) {
    ## version.unf defined per default
    ## Tetrad CFCI works with version.unf=c(1,2)
    sk <- pc.cons.intern(skel, suffStat, indepTest, alpha, verbose=verbose, version.unf=c(1,2))
    tripleList <- sk$unfTripl
    ## cat("\nfirst conservative list=", tripleList,"\n")
    ## vers <- sk$vers
  }
  if (doPdsep) {
    if (verbose) {
      cat("\nCompute PDSEP\n=============\ncompute collider...done\n")
    }
    ## it will be "pdsep"
    pdsepRes <- pdsep(skel@graph, suffStat, indepTest, p, sepset,
                          pMax, NAdelete, verbose=verbose, alpha, unfVect=tripleList, biCC=biCC)
    G <- pdsepRes$G
    sepset <- pdsepRes$sepset
    pMax <- pdsepRes$pMax
    allPdsep <- pdsepRes$allPdsep
    n.edgetestsPD <- pdsepRes$n.edgetests
    max.ordPD <- pdsepRes$max.ord

    ## conservative version
    if (conservative[2]) {
      ## compute conservative again, because skelet may have changed
      colnames(G) <- rownames(G) <- labels
      Gobject <- as(G, "graphNEL")
      tmp.pdsep <- new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
                       max.ord = as.integer(max.ordSKEL), n.edgetests = n.edgetestsSKEL,
                       sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
      sk.pdsep <- pc.cons.intern(tmp.pdsep, suffStat, indepTest, alpha, verbose=verbose, version.unf=c(1,2))
      tripleList <- sk.pdsep$unfTripl
      ## cat("\nsecond conservative list=", tripleList,"\n")
      ## vers <- sk.pdsep$vers
    }
  }
  else {
    n.edgetestsPD <- 0
    max.ordPD <- 0
    allPdsep <- vector("list", p)
  }

  ## if(verbose) { cat("Final graph adjacency matrix:\n"); print(symnum(G)) }

##################################################
  ## Make pcAlgo object
##################################################
  ## transform matrix to graph object
  if (sum(G) == 0) {
    Gobject <- new("graphNEL", nodes = labels)
  } else {
    colnames(G) <- rownames(G) <- labels
    Gobject <- as(G, "graphNEL")
  }
  tmp <- new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
             max.ord = as.integer(max.ordSKEL), n.edgetests = n.edgetestsSKEL,
             sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
  if (verbose)
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules),
        "\nCompute collider:\n")
  res <-
    if (numEdges(tmp@graph) > 0)
      ## it will be "udag2pag"
      udag2pag(gInput = tmp, rules = rules, verbose = verbose,
               unfVect= if(cons.rules) tripleList)
    else G

##################################################
  ## Return FCI object
##################################################
  new("fciAlgo",
      amat = res, call = cl, n = integer(0),
      max.ord = as.integer(max.ordSKEL),
      max.ordPDSEP = as.integer(max.ordPD),
      n.edgetests = n.edgetestsSKEL,
      n.edgetestsPDSEP = n.edgetestsPD,
      sepset = sepset, pMax = pMax, allPdsep = allPdsep)
}


gSquareBin <- function(x, y, S, dm, verbose=FALSE, adaptDF=FALSE)
{
  ## Purpose: G^2 statistic to test for (conditional) independence
  ##          of X and Y given S
  ## -------------------------------------------------------------------------
  ## Arguments:
  ## - x,y,S: Are x,y conditionally independent given S (S can be NULL)?
  ## - dm: data matrix (rows: samples, columns: variables) with binary entries
  ## - verbose: if TRUE, some additional info is outputted during the
  ##            computations
  ## - adaptDF: lower the degrees of freedom by one for each zero count.
  ##            The value for the DF cannot go below 1.
  ## -------------------------------------------------------------------------

  if(verbose) cat('\nEdge ',x,' -- ',y,' with subset: ',S,'\n')

  n <- dim(dm)[1] # nr of samples
  lenS <- length(S)

  ## degrees of freedom assuming no structural zeros
  df <- 2^lenS

  if (n < 10*df) { ## not enough samples to perform the test, assume
    ## independence
    if (verbose) warning("Not enough samples...\n")
    return( 1 )   ## gerster prob=0
  }
  ## else --  enough data to perform the test
  if(lenS < 6) {

    if (lenS == 0) {
      nij <- array(0,c(2,2))
      for (i in 1:2) {
        for (j in 1:2) {
          nij[i,j] <- sum(dm[,x] == i-1 &
                          dm[,y] == j-1)
        }
      }
      ## marginal counts
      t.X <- rowSums(nij)
      dim(t.X) <- c(length(t.X),1)
      t.Y <- colSums(nij)
      dim(t.Y) <- c(1,length(t.Y))

      ## compute G^2
      dij <- t.X %*% t.Y                # s_ia * s_jb
      t.G2 <- 2*nij*log(nij*n/dij)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    }                                   #end lenS=0
    else if (lenS==1) {
      ## a.pos <- sort(c(x,y,S))

      nijk <- array(0,c(2,2,2))
      for(i in 1:2) for(j in 1:2) for(k in 1:2) {
        nijk[i,j,k] <- sum((dm[,x]==i-1)&
                           (dm[,y]==j-1)&
                           (dm[,S]==k-1))
      }

      alt <- c(x,y,S)
      c <- which(alt==S)
      nik <- apply(nijk,c,rowSums)
      njk <- apply(nijk,c,colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0,c(2,2,2))
      if(c==3){
        for (k in 1:2) {
          t.X <- array(nik[,k],dim=c(dim(nik)[1],1))
          t.Y <- array(njk[,k],dim=c(1,dim(njk)[1]))
          t.dijk <- t.X %*% t.Y
          t.log[,,k] <- nijk[,,k]*nk[k]/t.dijk
        }
      } else if(c==1) {
        for (k in 1:2) {
          t.X <- array(nik[,k],dim=c(dim(nik)[1],1))
          t.Y <- array(njk[,k],dim=c(1,dim(njk)[1]))
          t.dijk <- t.X %*% t.Y
          t.log[k,,] <- nijk[k,,]*nk[k]/t.dijk
        }
      } else { ## c == 2 (?)
        for (k in 1:2) {
          t.X <- array(nik[,k],dim=c(dim(nik)[1],1))
          t.Y <- array(njk[,k],dim=c(1,dim(njk)[1]))
          t.dijk <- t.X %*% t.Y
          t.log[,k,] <- nijk[,k,]*nk[k]/t.dijk
        }
      }

      t.G2 <- 2 * nijk * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)

    }                                   # end lenS=1
    else if(lenS==2) {
      ## a.pos <- sort(c(x,y,S))

      nijk2 <- array(NA,c(2,2,2,2))
      for(i in 1:2) for(j in 1:2) for(k in 1:2) for(l in 1:2){
        nijk2[i,j,k,l] <- sum((dm[,x]==i-1)&(dm[,y]==j-1)&(dm[,S[1]]==k-1)&(dm[,S[2]]==l-1))
      }

      ## alt <- c(x,y,S)

      nijk <- array(NA,c(2,2,4))
      for(i in 1:2) for(j in 1:2){
        nijk[,,2*(i-1)+j] <- nijk2[,,i,j]
      }

      nik <- apply(nijk,3,rowSums)
      njk <- apply(nijk,3,colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0,c(2,2,4))
      for (k in 1:4) {
        t.X <- array(nik[,k],dim=c(dim(nik)[1],1))
        t.Y <- array(njk[,k],dim=c(1,dim(njk)[1]))
        t.dijk <- t.X %*% t.Y
        t.log[,,k] <- nijk[,,k]*nk[k]/t.dijk
      }

      t.G2 <- 2 * nijk * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    }                                   #end lenS=2
    else if(lenS==3) {
      nijk <- array(NA,c(2,2,8))
      for(i1 in 1:2) for(i2 in 1:2) for(i3 in 1:2) for(i4 in 1:2) for(i5 in 1:2){
        nijk[i1,i2,4*(i3-1)+2*(i4-1)+i5] <-
          sum((dm[,x]==i1-1)&(dm[,y]==i2-1)&(dm[,S[1]]==i3-1)&(dm[,S[2]]==i4-1)&(dm[,S[3]]==i5-1))
      }

      nik <- apply(nijk,3,rowSums)
      njk <- apply(nijk,3,colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0,c(2,2,8))
      for (k in 1:8){
        t.X <- array(nik[,k],dim=c(dim(nik)[1],1))
        t.Y <- array(njk[,k],dim=c(1,dim(njk)[1]))
        t.dijk <- t.X %*% t.Y
        t.log[,,k] <- nijk[,,k]*nk[k]/t.dijk
      }

      t.G2 <- 2 * nijk * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    }                                   #end lenS=3
    else if(lenS==4) {
      nijk <- array(NA,c(2,2,16))
      for(i1 in 1:2) for(i2 in 1:2) for(i3 in 1:2) for(i4 in 1:2) for(i5 in 1:2) for(i6 in 1:2){
        nijk[i1,i2,8*(i3-1)+4*(i4-1)+2*(i5-1)+i6] <-
          sum((dm[,x]==i1-1) &
              (dm[,y]==i2-1) &
              (dm[,S[1]]==i3-1) &
              (dm[,S[2]]==i4-1) &
              (dm[,S[3]]==i5-1) &
              (dm[,S[4]]==i6-1))
      }

      nik <- apply(nijk,3,rowSums)
      njk <- apply(nijk,3,colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0,c(2,2,16))
      for (k in 1:16){
        t.X <- array(nik[,k],dim=c(dim(nik)[1],1))
        t.Y <- array(njk[,k],dim=c(1,dim(njk)[1]))
        t.dijk <- t.X %*% t.Y
        t.log[,,k] <- nijk[,,k]*nk[k]/t.dijk
      }

      t.G2 <- 2 * nijk * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    }                                   #end lens=4
    else if(lenS==5) {
      nijk <- array(NA,c(2,2,32))
      for(i1 in 1:2) for(i2 in 1:2) for(i3 in 1:2) for(i4 in 1:2) for(i5 in 1:2) for(i6 in 1:2) for(i7 in 1:2){
        nijk[i1,i2,16*(i3-1)+8*(i4-1)+4*(i5-1)+2*(i6-1)+i7] <-
          sum((dm[,x]==i1-1)&
              (dm[,y]==i2-1)&
              (dm[,S[1]]==i3-1)&
              (dm[,S[2]]==i4-1)&
              (dm[,S[3]]==i5-1)&
              (dm[,S[4]]==i6-1)&
              (dm[,S[5]]==i7-1))
      }

      nik <- apply(nijk,3,rowSums)
      njk <- apply(nijk,3,colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0,c(2,2,32))
      for (k in 1:32){
        t.X <- array(nik[,k],dim=c(dim(nik)[1],1))
        t.Y <- array(njk[,k],dim=c(1,dim(njk)[1]))
        t.dijk <- t.X %*% t.Y
        t.log[,,k] <- nijk[,,k]*nk[k]/t.dijk
      }

      t.G2 <- 2 * nijk * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    }                                   # end lenS=5

  } else { # --- lenS >= 6 ---------
    nijk <- tijk <- array(0,c(2,2,1))
    ## first sample 'by hand' to avoid if/else in the for-loop
    i <- dm[1,x]+1
    j <- dm[1,y]+1
    ## create directly a list of all k's
    k <- NULL
    sapply(as.list(S),function(x){k <<- cbind(k,dm[,x]+1);return(TRUE)})
    ## first set of subset values
    parents.count <- 1 ## counter variable corresponding to the number
    ## of value combinations for the subset varibales
    ## observed in the data
    parents.val <- t(k[1,])
    nijk[i,j,parents.count] <- 1        # cell counts

    ## Do the same for all other samples. If there is already a table
    ## for the subset values of the sample, increase the corresponding
    ## cell count. If not, create a new table and set the corresponding
    ## cell count to 1.
    for (it.sample in 2:n) {
      flag <- 0
      i <- dm[it.sample,x]+1
      j <- dm[it.sample,y]+1
      ## comparing the current values of the subset variables to all
      ## already existing combinations of subset variables values
      t.comp <- t(parents.val[1:parents.count,])==k[it.sample,]
      ## Have to be careful here. When giving dimension to a list,
      ## R fills column after column, and NOT row after row.
      dim(t.comp) <- c(lenS,parents.count)
      for (it.parents in 1:parents.count) {
        ## check if the present combination of value alreay exists
        if(all(t.comp[,it.parents])) {
          ## if yes, increase the corresponding cell count
          nijk[i,j,it.parents] <- nijk[i,j,it.parents] + 1
          flag <- 1
          break
        }
      }
      ## if the combination of subset values is new...
      if (flag==0) {
        if (verbose)
          cat('\n Adding a new combination of parents at sample ',
              it.sample,'\n')
        ## ...increase the number of subset 'types'
        parents.count <- parents.count + 1
        ## ...add the new subset to the others
        parents.val <- rbind(parents.val,k[it.sample,])
        ## ...update the cell counts (add new array)
        nijk <- abind(nijk,array(0,c(2,2,1)))
        nijk[i,j,parents.count] <- 1
      }## end if(flag==0)
    }## end for(it.sample ..)

    nik <- apply(nijk,3,rowSums)
    njk <- apply(nijk,3,colSums)
    nk <- colSums(njk)
    ## compute G^2
    t.log <- array(0,c(2,2,parents.count))
    for (k in 1:parents.count) {
      t.X <- array(nik[,k],dim=c(dim(nik)[1],1))
      t.Y <- array(njk[,k],dim=c(1,dim(njk)[1]))
      t.dijk <- t.X %*% t.Y
      t.log[,,k] <- nijk[,,k]*nk[k]/t.dijk
    }
    t.G2 <- 2 * nijk * log(t.log)
    t.G2[is.nan(t.G2)] <- 0
    G2 <- sum(t.G2)
  }

  if (adaptDF && lenS > 0) {
    ## lower the degrees of freedom according to the amount of zero
    ## counts; add zero counts corresponding to the number of parents
    ## combinations that are missing
    zero.counts <- sum(nijk == 0) + 4*(2^lenS-dim(nijk)[3])
    df <- max(1, df-zero.counts)
  }

  pchisq(G2, df, lower.tail=FALSE)# i.e. == 1 - P(..)

}## gSquareBin()

gSquareDis <- function(x, y, S, dm, nlev, verbose=FALSE,adaptDF=FALSE){

  ## Purpose: G^2 statistic to test for (conditional) independence
  ##          of X and Y given S
  ## -------------------------------------------------------------------
  ## Arguments:
  ## - x,y,S: Are x,y conditionally independent given S (S can be NULL)?
  ## - dm: data matrix (rows: samples, columns: variables) with
  ##       discrete entries
  ## - nlev: vector with numbers of levels for each variable
  ## - verbose: if TRUE, some additional info is outputted during the
  ##            computations
  ## - adaptDF: lower the degrees of freedom by one for each zero count.
  ##            The value for the DF cannot go below 1.
  ## -------------------------------------------------------------------

  if(verbose) cat('\nEdge ',x,' -- ',y,' with subset: ',S,'\n')

  n <- dim(dm)[1] # nr of samples
  lenS <- length(S)

  ## degrees of freedom assuming no structural zeros
  df <- (nlev[x]-1)*(nlev[y]-1)*prod(nlev[S])

  if (n < 10*df) { ## not enough samples to perform the test, assume
    ## independence
    if (verbose) warning("Not enough samples...\n")
    return( 1 )   ## gerster prob=0
  }
  ## else --  enough data to perform the test

  if(lenS < 5) { #bei gSquareBin lenS<6

    if (lenS == 0) {
      nij <- array(0,c(nlev[x],nlev[y]))
      for (i in 1:nlev[x]) {
        for (j in 1:nlev[y]) {
          nij[i,j] <- sum((dm[,x]==i-1)&(dm[,y]==j-1))
        }
      }
      ## marginal counts
      t.X <- rowSums(nij)
      dim(t.X) <- c(length(t.X),1)
      t.Y <- colSums(nij)
      dim(t.Y) <- c(1,length(t.Y))

      ## compute G^2
      dij <- t.X %*% t.Y                # s_ia * s_jb
      t.log <- nij*n/dij
      t.G2 <- 2*nij*log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    }                                   # end lenS=0
    else if(lenS==1) {
      nijk <- array(0,c(nlev[x],nlev[y],nlev[S]))
      for(i in 1:nlev[x]) for(j in 1:nlev[y]) for(k in 1:nlev[S]) {
        nijk[i,j,k] <- sum((dm[,x]==i-1)&(dm[,y]==j-1)&(dm[,S]==k-1))
      }

      nik <- apply(nijk,3,rowSums)
      njk <- apply(nijk,3,colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0,c(nlev[x],nlev[y],prod(nlev[S])))
      for (k in 1:prod(nlev[S])) {
        t.X <- array(nik[,k],dim=c(dim(nik)[1],1))
        t.Y <- array(njk[,k],dim=c(1,dim(njk)[1]))
        t.dijk <- t.X %*% t.Y
        t.log[,,k] <- nijk[,,k]*nk[k]/t.dijk
      }

      t.G2 <- 2 * nijk * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    }                                   # end lenS=1
    else if(lenS == 2) {
      nijk <- array(NA,c(nlev[x],nlev[y],nlev[S[1]]*nlev[S[2]]))
      for(i in 1:nlev[x]) for(j in 1:nlev[y]) for(k in 1:nlev[S[1]]) for(l in 1:nlev[S[2]]){
        nijk[i,j,nlev[S[2]]*(k-1)+l] <- sum((dm[,x]==i-1)&(dm[,y]==j-1)&(dm[,S[1]]==k-1)&(dm[,S[2]]==l-1))
      }

      nik <- apply(nijk,3,rowSums)
      njk <- apply(nijk,3,colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0,c(nlev[x],nlev[y],prod(nlev[S])))
      for (k in 1:prod(nlev[S])) {
        t.X <- array(nik[,k],dim=c(dim(nik)[1],1))
        t.Y <- array(njk[,k],dim=c(1,dim(njk)[1]))
        t.dijk <- t.X %*% t.Y
        t.log[,,k] <- nijk[,,k]*nk[k]/t.dijk
      }

      t.G2 <- 2 * nijk * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    }                                   #end lenS=2
    else if(lenS == 3) {
        nijk <- array(NA,c(nlev[x],nlev[y],prod(nlev[S])))
        for(i1 in 1:nlev[x]) for(i2 in 1:nlev[y]) for(i3 in 1:nlev[S[1]]) for(i4 in 1:nlev[S[2]]) for(i5 in 1:nlev[S[3]]){
          nijk[i1,i2,nlev[S[3]]*nlev[S[2]]*(i3-1)+nlev[S[3]]*(i4-1)+i5] <- sum((dm[,x]==i1-1)&(dm[,y]==i2-1)&(dm[,S[1]]==i3-1)&(dm[,S[2]]==i4-1)&(dm[,S[3]]==i5-1))
        }

        nik <- apply(nijk,3,rowSums)
        njk <- apply(nijk,3,colSums)
        nk <- colSums(njk)

        ## compute G^2
        t.log <- array(0,c(nlev[x],nlev[y],prod(nlev[S])))
        for (k in 1:prod(nlev[S])){
          t.X <- array(nik[,k],dim=c(dim(nik)[1],1))
          t.Y <- array(njk[,k],dim=c(1,dim(njk)[1]))
          t.dijk <- t.X %*% t.Y
          t.log[,,k] <- nijk[,,k]*nk[k]/t.dijk
        }

        t.G2 <- 2 * nijk * log(t.log)
        t.G2[which(is.nan(t.G2),arr.ind=TRUE)] <- 0
        G2 <- sum(t.G2)
      } #end lenS=3
      else if(lenS == 4) {
        nijk <- array(NA,c(nlev[x],nlev[y],prod(nlev[S])))
        for(i1 in 1:nlev[x]) for(i2 in 1:nlev[y]) for(i3 in 1:nlev[S[1]]) for(i4 in 1:nlev[S[2]]) for(i5 in 1:nlev[S[3]]) for(i6 in 1:nlev[S[4]]){
          nijk[i1,i2,nlev[S[4]]*nlev[S[3]]*nlev[S[2]]*(i3-1)+nlev[S[4]]*nlev[S[3]]*(i4-1)+nlev[S[4]]*(i5-1)+i6] <- sum((dm[,x]==i1-1)&(dm[,y]==i2-1)&(dm[,S[1]]==i3-1)&(dm[,S[2]]==i4-1)&(dm[,S[3]]==i5-1)&(dm[,S[4]]==i6-1))
        }

        nik <- apply(nijk,3,rowSums)
        njk <- apply(nijk,3,colSums)
        nk <- colSums(njk)

        ## compute G^2
        t.log <- array(0,c(nlev[x],nlev[y],prod(nlev[S])))
        for (k in 1:prod(nlev[S])){
          t.X <- array(nik[,k],dim=c(dim(nik)[1],1))
          t.Y <- array(njk[,k],dim=c(1,dim(njk)[1]))
          t.dijk <- t.X %*% t.Y
          t.log[,,k] <- nijk[,,k]*nk[k]/t.dijk
        }

        t.G2 <- 2 * nijk * log(t.log)
        t.G2[is.nan(t.G2)] <- 0
        G2 <- sum(t.G2)
      } #end lens=4
    }
    else{ #  lenS > 4 (gSquareBin: lenS>5)
      nijk <- tijk <- array(0,c(nlev[x],nlev[y],1))
      ## first sample 'by hand' to avoid if/else in the for-loop
      i <- dm[1,x]+1
      j <- dm[1,y]+1
      ## create directly a list of all k's
      k <- NULL
      sapply(as.list(S),(function(x){k <<- cbind(k,dm[,x]+1);return(TRUE)}))
      ## first set of subset values
      parents.count <- 1 ## counter variable corresponding to the number
      ## of value combinations for the subset varibales
      ## observed in the data
      parents.val <- t(k[1,])
      nijk[i,j,parents.count] <- 1 # cell counts

      ## Do the same for all other samples. If there is already a table
      ## for the subset values of the sample, increase the corresponding
      ## cell count. If not, create a new table and set the corresponding
      ## cell count to 1.
      for (it.sample in 2:n) {
        flag <- 0
        i <- dm[it.sample,x]+1
        j <- dm[it.sample,y]+1
        ## comparing the current values of the subset variables to all
        ## already existing combinations of subset variables values
        t.comp <- t(parents.val[1:parents.count,])==k[it.sample,]
                                        # Have to be careful here. When giving dimension to a list,
                                        # R fills column after column, and NOT row after row.
        dim(t.comp) <- c(lenS,parents.count)
        for (it.parents in 1:parents.count) {
          ## check if the present combination of value alreay exists
          if(all(t.comp[,it.parents])) {
            ## if yes, increase the corresponding cell count
            nijk[i,j,it.parents] <- nijk[i,j,it.parents] + 1
            flag <- 1
            break
          }
        }# end for(it.parents...)
        ## if the combination of subset values is new...
        if (flag==0) {
          if (verbose) {
            cat('\n Adding a new combination of parents at sample ',
                it.sample,'\n')
          }
          ## ...increase the number of subset 'types'
          parents.count <- parents.count + 1
          ## ...add the new subset to the others
          parents.val <- rbind(parents.val,k[it.sample,])
          ## ...update the cell counts (add new array)
          nijk <- abind(nijk,array(0,c(nlev[x],nlev[y],1)))
          nijk[i,j,parents.count] <- 1
        } # end if(flag==0)
      }

      nik <- apply(nijk,3,rowSums)
      njk <- apply(nijk,3,colSums)
      nk <- colSums(njk)
      ## compute G^2
      t.log <- array(0,c(nlev[x],nlev[y],parents.count))
      for (k in 1:parents.count) {
        t.X <- array(nik[,k],dim=c(dim(nik)[1],1))
        t.Y <- array(njk[,k],dim=c(1,dim(njk)[1]))
        t.dijk <- t.X %*% t.Y
        t.log[,,k] <- nijk[,,k]*nk[k]/t.dijk
      }
      t.G2 <- 2 * nijk * log(t.log)
      t.G2[which(is.nan(t.G2),arr.ind=TRUE)] <- 0
      G2 <- sum(t.G2)
    }

    if (adaptDF && lenS>0) {
      ## lower the degrees of freedom according to the amount of
      ## zero counts
      if(lenS==0){
        zero.counts <- length(which(nij==0))
      }
      else{
        zero.counts <- length(which(nijk == 0))
        zero.counts <- zero.counts + 4*(2^lenS-dim(nijk)[3])
      }
      ## add zero counts corresponding to the number of parents
      ## combinations that are missing
      df <- max((df-zero.counts),1)
    } # end adaptDF

    pchisq(G2, df, lower.tail=FALSE)# i.e. == 1 - P(..)

}## gSquareDis()

gaussCItest <- function(x,y,S,suffStat) {
  ## suffStat$C: correlation matrix
  ## suffStat$n: sample size
  z <- zStat(x,y,S, suffStat$C, suffStat$n)
  2*pnorm(abs(z), lower.tail=FALSE)
}


dsep <- function(a,b,S,g,john.pairs=NA)
{
  ## Purpose: Are the set a and the set b d-separeted given the set S?
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - a,b,S: vectors of node names
  ## - g: graphNEL object
  ## - john.pairs: matrix from johnson.all.pairs.sp
  ## ----------------------------------------------------------------------
  ## Value:
  ## Boolean decision
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch

  ## Check that g is a DAG
  amatTmp <- wgtMatrix(g)
  amatTmp[amatTmp!=0] <- 1
  if (max(amatTmp+t(amatTmp))>1) stop("dsep: Undirected edge in input graph!")
  p <- numNodes(g)
  ## build node union of a,b,S
  if (any(is.na(john.pairs))) john.pairs <- johnson.all.pairs.sp(g)
  if (length(S) > 0) {
    nodeUnion <- c(a,b,S)
  } else {
    nodeUnion <- c(a,b)
  } ## if (length(S) > 0)
  my.nodes <- nodes(g)

  ## find ancestor graph of nodeUnion
  anc.set <- NULL
  for (i in 1:p) {
    desc.nodes <- my.nodes[which(john.pairs[i,]<Inf)]
    if (any(desc.nodes %in% nodeUnion)) anc.set <- c(anc.set,my.nodes[i])
  } ## for (i in 1:p)
  gS <- subGraph(anc.set,g)

  ## Moralize in amatM
  amat <- wgtMatrix(gS, transp=FALSE)
  amat[amat!=0] <- 1
  amatM <- amat
  ind <- which(amat==1,arr.ind=TRUE)
  if (length(ind)>0) {
    for (i in 1:dim(ind)[1]) {
      ## input is guaranteed to be directed
      x <- ind[i,1]
      y <- ind[i,2] ## x->y
      allZ <- setdiff(which((amat[y,]==0)&(amat[,y]==1)),x) ## x->y<-z
      for (z in allZ) {
        if ((amat[x,z]==0) && (amat[z,x]==0)) amatM[x,z] <- 1 ## moralize
      } ## for (z)
    } ## for (i in 1:dim(ind)[1])

    ## make undirected graph
    ## (up to now, there is NO undirected edge -> just add t(amat))
    gSM <- as(amatM+t(amatM),"graphNEL")

    ## check separation
    if (length(S)>0) {
      res <- separates(a,b,S,gSM)
    } else {
      bfs.res <- bfs(gSM,a)
      if (!is.list(bfs.res)) {
        res <- !(b %in% bfs.res)
      } else {
        res <- !(b %in% bfs.res[[1]])
      }
    } ## if (length(S)>0)
  } else { ## if (length(ind)>0)
    ## if no edge in graph, nodes are d-separated
    res <- TRUE
  }
  res
}


## Orakel
dsepTest <- function(x,y,S,suffStat) {
  ## suffStat$g: True graph (graphNEL suffStatect)
  ## suffStat$jp: johnson all pairs
  ## Return "P-value"
  ## p == 0: keep edge / d-connected
  ## p == 1; drop edge / d-separated
  g <- suffStat$g
  jp <- suffStat$jp
  stopifnot(is(g, "graph"))

  if (x==y || (x %in% S) || (y %in% S)) {
    0
  } else {
    dSepTrue <- dsep(a = as.character(x), b = as.character(y),
                     S = if(length(S)) as.character(S),# else NULL
                     g = g, john.pairs = jp)
    if (dSepTrue) ## delete edge
      1 else 0
  }
}

## disCItest
disCItest <- function(x,y,S,suffStat) {
  dm <- suffStat$dm
  nlev <- suffStat$nlev
  adaptDF <- suffStat$adaptDF
  ## P-value:
  gSquareDis(x = x, y = y, S = S, dm = dm, nlev = nlev,
             verbose = FALSE, adaptDF = adaptDF)
}

## binCItest
binCItest <- function(x,y,S,suffStat) {
  dm <- suffStat$dm
  adaptDF <- suffStat$adaptDF
  ## P-value:
  gSquareBin(x = x, y = y, S = S, dm = dm, verbose = FALSE,
             adaptDF = adaptDF)
}

pdsep <- function (skel, suffStat, indepTest, p, sepset, pMax, NAdelete = TRUE,
                       verbose = FALSE, alpha, unfVect=NULL, biCC=FALSE)
{
  ## Purpose: Compute possible D-sep for each node and adapt graph accordingly
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - skel: Graph object returned by function skeleton
  ## - sepset: Sepset that was used for finding the skeleton
  ## - pMax: Maximal p-values during estimation of skeleton
  ## - testType: Type of cond. independence test used
  ## - gTrue: True graph [graph object]
  ## ----------------------------------------------------------------------
  ## Value:
  ## - G: Updated boolean adjacency matrix
  ## - sepset: Updated sepsets
  ## - pMax: Updated pMax
  ## - allPdsep: Possible d-sep for each node [list]
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  9 Dec 2009, 16:01
  G <- (as(skel, "matrix") != 0)
  n.edgetests <- rep(0, 1000)
  max.ord <- 0
  allPdsep <- vector("list", p)
  if (biCC) {
    conn.comp <- lapply(biConnComp(skel), as.numeric)
  }
  if (any(G)) {
    amat <- G
    amat[amat==TRUE] <- 1
    amat[amat==FALSE] <- 0
    ind <- which(amat==1, arr.ind=TRUE)
    ## Orient colliders
    for (i in 1:dim(ind)[1]) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(amat[y, ] != 0), x)
      if (length(allZ) > 0) {
        for (j in seq_along(allZ)) {
          z <- allZ[j]
          if (amat[x, z] == 0 &&
              !((y %in% sepset[[x]][[z]]) |
                (y %in% sepset[[z]][[x]]))) {
            if (length(unfVect)==0) { ## normal version:
              amat[x, y] <- amat[z, y] <- 2
              if (verbose) cat(x,"o->", y, "<-o", z, "\n")
            }
            else { ## conservative version
              ## check if x-y-z is faithful
              if (!any(unfVect == triple2numb(p,x,y,z)) &&
                  !any(unfVect == triple2numb(p,z,y,x))) {
                amat[x, y] <- amat[z, y] <- 2
                if (verbose) cat(x,"o->", y, "<-o", z, "\n")
              }
            }
          }
        } ## for( j )
      }
    } ## for( i )

    for (x in 1:p) {
      if (any(x.has <- amat[x, ] != 0)) {
        allPdsep[[x]] <- qreach(x, amat)
        tf1 <- setdiff(allPdsep[[x]], x)
        if (verbose) {
          cat("Possible D-Sep of", x, " is: ", allPdsep[[x]], "\n")
        }
        adj.x <- (1:p)[amat[x, ] != 0]
        for (y in adj.x) {
          tf <- setdiff(tf1, y)
          diff.set <- setdiff(tf, adj.x)
          ## bi-connected components
          if (biCC) {
            index.conncomp <- 0
            found.conncomp <- FALSE
            while ((!found.conncomp) & (index.conncomp < length(conn.comp))) {
              index.conncomp <- index.conncomp + 1
              if (x %in% conn.comp[[index.conncomp]] && y %in% conn.comp[[index.conncomp]]) {
                found.conncomp <- TRUE
              }
            }
            bi.conn.comp <- setdiff(conn.comp[[index.conncomp]],c(x,y))
            tmp.tf <- intersect(tf,bi.conn.comp)
            tf <- tmp.tf
            if (verbose) {
              cat("Possible D-Sep of", x,"and", y,"intersected with the biconnected component is", tf, "\n")
            }
            ## if (length(tf)>15) {
              ## cat("Size of pds bigger than 15: break the survey between",x,"and",y,"\n")
              ## break
            ##}
          }
          if (length(diff.set) > 0) {
            done <- FALSE
            ord <- 0
            while (!done && ord < length(tf)) {
              ord <- ord + 1
              if (ord > max.ord)
                max.ord <- ord
              if (ord == 1) {
                for (j in seq_along(diff.set)) {
                  pval <- indepTest(x, y, diff.set[j], suffStat)
                  n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                  if (is.na(pval))
                    pval <- if(NAdelete) 1 else 0
                  if (pval > pMax[x, y])
                    pMax[x, y] <- pval
                  if (pval >= alpha) {
                    amat[x, y] <- amat[y, x] <- 0
                    sepset[[x]][[y]] <- sepset[[y]][[x]] <- diff.set[j]
                    done <- TRUE
                    if (verbose)
                      cat("x=", x, " y=", y, " S=", diff.set[j],
                          ": pval =", pval, "\n")
                    break
                  }
                }
              }
              else if (ord <= length(adj.x)) {
                  tmp.combn <- combn(tf, ord)
                  for (k in 1:dim(tmp.combn)[2]) {
                    tmp.ii <- tmp.combn[, k] %in% adj.x
                    if (any(!tmp.ii)) {
                      pval <- indepTest(x, y, tmp.combn[, k], suffStat)
                      n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                      if (is.na(pval))
                        pval <- if(NAdelete) 1 else 0
                      if (pval > pMax[x, y])
                        pMax[x, y] <- pval
                      if (pval >= alpha) {
                        amat[x, y] <- amat[y, x] <- 0
                        sepset[[x]][[y]] <- sepset[[y]][[x]] <- tmp.combn[,k]
                        done <- TRUE
                        if (verbose) {
                          cat("x=", x, " y=", y, " S=",
                              tmp.combn[, k], ": pval =",
                              pval, "\n")
                        }
                        break
                      }
                    }
                  }
              }
              else {
                  ## check all combinations; no combination has been
                  ## tested before, since ord > adj.x
                  tmp.combn <- combn(tf,ord)
                  for (k in 1:dim(tmp.combn)[2]) {
                    pval <- indepTest(x, y, tmp.combn[, k], suffStat)
                    n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                    if (is.na(pval))
                      pval <- if(NAdelete) 1 else 0
                    if (pval > pMax[x, y])
                      pMax[x, y] <- pval
                    if (pval >= alpha) {
                      amat[x, y] <- amat[y, x] <- 0
                      sepset[[x]][[y]] <- sepset[[y]][[x]] <- tmp.combn[,k]
                      done <- TRUE
                      if(verbose)
                        cat("x=", x, " y=", y, " S=", tmp.combn[, k], ": pval =", pval, "\n")
                      break
                    }
                  }
              }
            }
          }
        }
      }
    }
    G[amat == 0] <- FALSE
    G[amat == 1] <- TRUE
    G[amat == 2] <- TRUE
  }
  list(G = G, sepset = sepset, pMax = pMax, allPdsep = allPdsep,
       max.ord = max.ord, n.edgetests = n.edgetests[1:(max.ord + 1)])
}


ida <- function(x.pos,y.pos,mcov,graphEst,method="local",
                y.notparent = FALSE, verbose=FALSE, all.dags=NA)
{
  ## Purpose: Estimate the causal effect of x on y; the graphEst and correlation
  ## matrix have to be precomputed; all DAGs can be precomputed;
  ## Orient undirected edges at x in a way so that no new collider
  ## is introduced
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x.pos, y.pos: Column of x and y in d.mat
  ## - mcov: Covariance matrix that was used to estimate graphEst
  ## - graphEst: Fit of PC Algorithm (semidirected)
  ## - method: "local" - local (all combinations of parents in regr.)
  ##           "global" - all DAGs
  ## - y.notparent: if TRUE, the effect of x <- y is ignored;
  ##                (remove y from all parents set pa1 or pa2)
  ##                if FALSE, the effect of x <- y is set to zero
  ## - verbose: if TRUE, details on regressions that were used
  ## - all.dags: All DAGs in the format of function allDags; if this is
  ##   available, no new function call allDags is done
  ## ----------------------------------------------------------------------
  ## Value: causal values
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 7 Jan 2010, 11:18

  tmpColl <- FALSE

  ## prepare adjMatrix and skeleton
  amat <- wgtMatrix(graphEst)
  amat[which(amat!=0)] <- 1 ## i->j if amat[j,i]==1
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel!=0] <- 1

  if (method=="local") {
##############################
    ## local method
    ## Main Input: mcov, graphEst
##############################
    ## find unique parents of x
    wgt.est <- (wgtMatrix(graphEst)!=0)
    if (y.notparent) {
      ## Direct edge btw. x.pos and y.pos towards y.pos
      wgt.est[x.pos, y.pos] <- FALSE
    }
    tmp <- wgt.est-t(wgt.est)
    tmp[which(tmp<0)] <- 0
    wgt.unique <- tmp
    pa1 <- which(wgt.unique[x.pos,]!=0)
    if (y.pos %in% pa1) {
      ## x is parent of y -> zero effect
      beta.hat <- 0
    } else { ## y.pos not in pa1
      ## find ambiguous parents of x
      wgt.ambig <- wgt.est-wgt.unique
      pa2 <- which(wgt.ambig[x.pos,]!=0)
      if (verbose) {
        cat("\n\nx=",x.pos,"y=",y.pos,"\n")
        cat("pa1=",pa1,"\n")
        cat("pa2=",pa2,"\n")
      }

      ## estimate beta
      if (length(pa2)==0) {
        beta.hat <- lm.cov(mcov,y.pos,c(x.pos,pa1))
        if (verbose) cat("Fit - y:",y.pos,"x:",c(x.pos,pa1),
                         "|b.hat=",beta.hat,"\n")
      } else {
        ## at least one undirected parent
        beta.hat <- NA
        ii <- 1

        ## no member of pa2
        pa2.f <- pa2
        pa2.t <- NA
        tmpColl <- check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)
        if (!tmpColl) {
          beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
          if (verbose) cat("Fit - y:",y.pos,"x:",c(x.pos,pa1),
                           "|b.hat=",beta.hat[ii],"\n")
        }
        ## exactly one member of pa2
        for (i2 in seq_along(pa2)) {
          pa2.f <- pa2[-i2]
          pa2.t <- pa2[i2]
          tmpColl <- check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)
          if (!tmpColl) {
            ii <-  ii+1
            if (y.pos %in% pa2.t) {
              beta.hat[ii] <- 0
            } else {
              beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1,pa2[i2]))
              if (verbose) cat("Fit - y:",y.pos,"x:",c(x.pos,pa1,pa2[i2]),
                               "|b.hat=",beta.hat[ii],"\n")
            } ## if (y.pos %in% pa2.t)
          } ## if (!tmpColl)
        } ## for (i2 in seq_along(pa2))

        ## higher order subsets of pa2
        if (length(pa2)>1) {
          for (i in 2:length(pa2)) {
            pa.tmp <- combn(pa2,i,simplify=TRUE)
            n.comb <- ncol(pa.tmp)
            for (j in 1:n.comb) {
              pa2.f <- setdiff(pa2,pa.tmp[,j])
              pa2.t <- pa.tmp[,j]
              tmpColl <- check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)
              if (!tmpColl) {
                ii <- ii+1
                if (y.pos %in% pa2.t) {
                  beta.hat[ii] <- 0
                } else {
                  beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1,pa.tmp[,j]))
                  if (verbose)
                    cat("Fit - y:",y.pos,"x:",c(x.pos,pa1,pa.tmp[,j]),
                        "|b.hat=",beta.hat[ii],"\n")
                }
              } ## if (!tmpColl)
            } ## for (j in 1:n.comb)
          } ## for (i in 2:length(pa2))
        } ## if (length(pa2)>1)
      } ## if (length(pa2) == 0)
    } ## if (y.pos %in% pa1)

  } else {
##############################
    ## global method
    ## Main Input: mcov, graphEst
##############################
    p <- numNodes(graphEst)
    am.pdag <- wgtMatrix(graphEst)
    am.pdag[am.pdag!=0] <- 1
    if (y.notparent) {
      ## Direct edge btw. x.pos and y.pos towards y.pos
      am.pdag[x.pos, y.pos] <- 0
    }

    ## find all DAGs if not provided externally
    ad <- if(is.na(all.dags)) allDags(am.pdag,am.pdag,NULL) else all.dags
    n.dags <- nrow(ad)
    beta.hat <- rep(NA,n.dags)
    for (i in 1:n.dags) {
      ## compute effect for every DAG
      gDag <- as(matrix(ad[i,],p,p),"graphNEL")
      ## path from y to x
      ## rev.pth <- RBGL::sp.between(gDag,as.character(y.pos),
      ##                    as.character(x.pos))[[1]]$path
      ## if (length(rev.pth)>1) {
      ## if reverse path exists, beta=0
      ##  beta.hat[i] <- 0
      ## } else {
      ## path from x to y
      ##       pth <- RBGL::sp.between(gDag,as.character(x.pos),
      ##                       as.character(y.pos))[[1]]$path
      ##   if (length(pth)<2) {
      ## sic! There is NO path from x to y
      ##   beta.hat[i] <- 0
      ## } else {
      ## There is a path from x to y
      wgt.unique <- t(matrix(ad[i,],p,p)) ## wgt.est is wgtMatrix of DAG
      pa1 <- which(wgt.unique[x.pos,]!=0)
      if (y.pos %in% pa1) {
        beta.hat[i] <- 0
      } else {
        beta.hat[i] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
        if (verbose) cat("Fit - y:",y.pos,"x:",c(x.pos,pa1),
                         "|b.hat=",beta.hat[i],"\n")
      }
    } ## for ( i  n.dags)
  } ## else : method = "global"
  beta.hat
}

idaFast <- function(x.pos,y.pos.set,mcov,graphEst)
{
  ## Purpose: Estimate the causal effect of x on each element in the
  ## set y using the local method; graphEst and correlation matrix
  ## have to be precomputed; orient
  ## undirected edges at x in a way so that no new collider is
  ## introduced; if there is an undirected edge between x and y, both directions are considered;
  ## i.e., y might be partent of x in which case the effect is 0.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x.pos, y.pos: Column of x and y in d.mat
  ## - mcov: Covariance matrix that was used to estimate graphEst
  ## - graphEst: Fit of PC Algorithm (semidirected)
  ## ----------------------------------------------------------------------
  ## Value: list of causal values; one list element for each element of
  ## y.pos.set
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 7 Jan 2010, 11:18

  ## prepare adjMatrix and skeleton
  amat <- wgtMatrix(graphEst)
  amat[which(amat!=0)] <- 1 ## i->j if amat[j,i]==1
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel!=0] <- 1

  ## find unique and ambiguous parents of x
  wgt.est <- (wgtMatrix(graphEst)!=0)
  tmp <- wgt.est-t(wgt.est)
  tmp[which(tmp<0)] <- 0
  wgt.unique <- tmp
  wgt.ambig <- wgt.est-wgt.unique
  pa1 <- which(wgt.unique[x.pos,]!=0)
  pa2 <- which(wgt.ambig[x.pos,]!=0)

  ## estimate beta
  if (length(pa2)==0) {
    beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1)) ####
    beta.tmp[y.pos.set %in% pa1] <- 0
    beta.hat <- cbind(beta.tmp)
  } else {    ## at least one undirected parent
    ## no member of pa2
    pa2.f <- pa2
    pa2.t <- NA
    beta.hat <-
      if (!check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
	beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1)) ####
	beta.tmp[y.pos.set %in% pa1] <- 0
	cbind(beta.tmp)
      } else NULL

    ## exactly one member of pa2
    for (i2 in seq_along(pa2)) {
      pa2.f <- pa2[-i2]
      pa2.t <- pa2[i2]
      if (!check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
        beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1,pa2.t)) ####
        beta.tmp[y.pos.set %in% c(pa1,pa2.t)] <- 0
        beta.hat <- cbind(beta.hat, beta.tmp)
      }
    } ## for (i2 in seq_along(pa2))

    ## higher order subsets of pa2
    if (length(pa2)>1) {
      for (i in 2:length(pa2)) {
        pa.tmp <- combn(pa2,i,simplify=TRUE)
        n.comb <- ncol(pa.tmp)
        for (j in seq_len(n.comb)) {
          pa2.t <- pa.tmp[,j]
          pa2.f <- setdiff(pa2, pa2.t)
          if (!check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
            beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1,pa2.t)) ####
            beta.tmp[y.pos.set %in% c(pa1,pa2.t)] <- 0
            beta.hat <- cbind(beta.hat, beta.tmp)
          }
        } ## for (j )
      } ## for (i )
    } ## if (length(pa2) > 1)
  } ## if .. else length(pa2) > 0)

  ## MM: for now, maybe in the future get sensible column names:
  colnames(beta.hat) <- NULL
  if (nrow(beta.hat) > 0) rownames(beta.hat) <- as.character(y.pos.set)
  beta.hat
}

legal.psep <- function(a,b,c,amat)
{
  ## Purpose: Is path a-b-c legal (either collider in b or a,b,c is triangle)
  ## !! a-b-c must be in a path !! this is not checked !!
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - a, b, c: nodes
  ## - amat: adj matrix (coding 0,1,2 for no edge, circle, arrowhead)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 29 Oct 2009, 16:57
  if(a == c || (a.b <- amat[a,b]) == 0 || amat[b,c] == 0)
    return(FALSE)
  ## else  a != c  and  amat[a,b] != 0  and   amat[b,c] != 0
  ## return  TRUE iff
  (amat[a,c] != 0 || ## triangle
   ## need not check [c,a], since there must be SOME edgemark !=0 at [a,c], if
   ## edge is present
   (a.b == 2 && amat[c,b] == 2)) ## a collider
}


qreach <- function(x,amat,verbose=FALSE)
{
  ## Purpose: Compute possible-d-sep(x) ("psep")
  ## !! The non-zero entries in amat must be symmetric !!
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x: node of which psep berechnet werden soll
  ## - amat: adjacency matrix
  ##         amat[i,j] = 0 iff no edge btw i,j
  ##         amat[i,j] = 1 iff i *-o j
  ##         amat[i,j] = 2 iff i *-> j
  ## - verbose: Show checked node sequence
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 29 Oct 2009, 11:54
  ## Stopping:
  ## =========
  ## At every iteration, Q get's reduced by one. It is only increased by
  ## at least one, if there are edges in amat.tmp and then, at least one
  ## edge in amat.tmp is removed. Edges are never inserted into amat.tmp.
  ## Thus, either Q or amat.tmp becomes empty and the loop stops.
  ## Runtime:
  ## ========
  ## At least O(|V|), since look up in adjacency matrix is made. Assume O(|E|)>O
  ## Every edge can be visited at most twice. At each visit, there are no
  ## more than max(deg(V_i)) neighboring edges to look at. Thus, the runtime is
  ## O(2*|E| * max(deg(V_i))) = O(|E|^2) [worst case]; O(|E|) if sparse in the
  ## sense that max(deg(V_i)) is constant.
  ## Correctness:
  ## ============
  ## (A) All nodes in PSEP have a path of legal triples from x.
  ## (B) If there is a legal path from x to y in amat, at least one of them
  ## is found and y is recorded in PSEP:
  ## Suppose there is a node y != x that has a legal path from x, but is not in
  ## PSEP. y cannot be in nbrs(x), because they are added and nothing is
  ## deleted from PSEP. Hence, there must be a node z which
  ## is in PSEP but has a neighbor w in amat that has a legal path from
  ## x but is not in PSEP.
  ## Assuming that the function legal is correct, and noting that at (*) all
  ## neighbors of z in amat.tmp (tmp!) are analyzed, it follows that w is not
  ## in adj(z) in amat.tmp (but in amat). Thus, w must have been removed
  ## before from adj(z). Because of (+), it could only have been removed if
  ## u-z-w was found legal at some point. But then, w would have been added
  ## to PSEP. This is a contradiction to the assumption that w is not in PSEP.

  ## check: quadratic; x in V; edgemarks ok; non-zeroes symmetric
  stopifnot((ncol(amat)==nrow(amat)),x<=ncol(amat),all(amat %in% c(0,1,2)),
            all((amat!=0)==(t(amat!=0))))
  amat.tmp <- amat

  nb <- which(amat[x,]!=0)
  Q <- nb
  P <- rep(x,length(Q))
  PSEP <- nb

  amat.tmp[x,nb] <- 0 ## delete edge to nbrs

  while(length(Q) > 0) {
    ## Invariants:
    ## ===========
    ## (A1) length(Q) == length(P) > 0
    ## (A2) non-zero in amat.tmp -> non-zero in amat [no non-zero added]
    ## (A3) Q[i] and P[i] are adjacent in amat [using (A2)]
    if (verbose) {
      cat("\n-------------","\n")
      cat("Queue Q:",Q,"\n")
      cat("Queue P:",P,"\n")
    }
    a <- Q[1]
    Q <- Q[-1]
    pred <- P[1] ## not empty because of (A1)
    P <- P[-1]
    if (verbose) cat("Select",pred,"towards",a,"\n")
    nb <- which(amat.tmp[a,]!=0) ## (*)
    if (verbose) cat("Check nbrs",nb,"\nLegal:")

    if (length(nb)>0) {
      for (i in seq_along(nb)) {
        b <- nb[i]
        ## Guaranteed: pred-a-b are a path because of (A3)
        if (lres <- legal.psep(pred,a,b,amat)) {
          amat.tmp[a,b] <- 0 ## remove b out of adj(a) in amat.tmp (+)
          Q <- c(Q,b)
          P <- c(P,a)
          PSEP <- c(PSEP,b)
        }
        if (verbose) cat(lres," ")
      }
    }
  }
  sort(unique(setdiff(PSEP,x)))
}

### Classes -->  AllClasses.R

##################################################
## Methods
##################################################
setMethod("summary", "pcAlgo",
          function(object) {
 	    cat("\nObject of class 'pcAlgo', from Call: \n",
                deparse(object@call),
 		"\n\nNmb. edgetests during skeleton estimation:\n")
            cat("===========================================\n")
            cat("Max. order of algorithm: ",object@max.ord,
                "\nNumber of edgetests from m = 0 up to m =",object@max.ord,
                ": ",object@n.edgetests)
            tmp <- object@graph@edgeL
            nbrs <- rep(0,length(tmp))
            for (i in seq_along(tmp)) {
              nbrs[i] <- length(tmp[[i]]$edges)
            }
            cat("\n\nGraphical properties of skeleton:\n")
            cat("=================================\n")
            cat("Max. number of neighbours: ",max(nbrs),
                "at node(s)",which(nbrs==max(nbrs)),
                "\nAvg. number of neighbours: ",mean(nbrs),"\n")
          })

setMethod("summary", "fciAlgo",
          function(object) {
 	    cat("Object of class 'fciAlgo'\n\n")
            cat("Call: \n=====\n", deparse(object@call))
            cat("\n\nNmb. edgetests during skeleton estimation:\n==========================================")
            cat("\nMax. order of algorithm: ",object@max.ord,
                "\nNumber of edgetests from m = 0 up to m =",object@max.ord,
                ": ",object@n.edgetests)
            cat("\n\nAdd. nmb. edgetests when using PDSEP:\n=====================================")
            cat("\nMax. order of algorithm: ",object@max.ordPDSEP,
                "\nNumber of edgetests from m = 0 up to m =",object@max.ordPDSEP,
                ": ",object@n.edgetestsPDSEP)

            myfun <- function(x) if(is.null(x)) NA else length(x)
            cat("\n\nSize distribution of SEPSET:")
            myTab <- table(sapply(object@sepset,function(x) sapply(x,myfun)),
                           useNA="always")
            print(myTab)

            cat("\nSize distribution of PDSEP:")
            print(table(sapply(object@allPdsep, length)))


          })


setMethod("show", "pcAlgo",
	  function(object) {
	    cat("Object of class 'pcAlgo', from Call: \n", deparse(object@call),"\n")
	    print(object@graph)
	    invisible(object)
	  })

setMethod("show", "fciAlgo",
	  function(object) {
	    cat("Object of class 'fciAlgo', from Call:", deparse(object@call),
                "\nAdjacency Matrix G:",
                "G[i,j] = 1/2/3 if edge mark of edge i-j at j is circle/head/tail.",
                "", sep="\n")
	    print(object@amat)
	    invisible(object)
	  })

setMethod("plot", signature(x = "pcAlgo"),
          function(x, y, main = NULL, zvalue.lwd = FALSE,
                   lwd.max = 7, labels = NULL, ...)
	{
          check.Rgraphviz()

          if(is.null(main))
              main <- deparse(x@call)
          attrs <- list()
          nodeAttrs <- list()
          if (!is.null(labels)) {
              attrs$node <- list(shape = "ellipse", fixedsize = FALSE)
              names(labels) <- nodes(x@graph)
              nodeAttrs$label <- labels
          }

          if (zvalue.lwd && numEdges(x@graph)!=0) {
              lwd.Matrix <- x@zMin
              lwd.Matrix <- ceiling(lwd.max*lwd.Matrix/max(lwd.Matrix))
              z <- agopen(x@graph,
                          name="lwdGraph",
                          nodeAttrs = nodeAttrs,
                          attrs = attrs)
              eLength <- length(z@AgEdge)
              for (i in 1:eLength) {
                  x <- as.numeric(z@AgEdge[[i]]@head)
                  y <- as.numeric(z@AgEdge[[i]]@tail)
                  z@AgEdge[[i]]@lwd <- lwd.Matrix[x,y]
              }
              Rgraphviz::plot(z, main = main, ...)
          } else {
              Rgraphviz::plot(x@graph, nodeAttrs = nodeAttrs, main = main,
                              attrs = attrs, ...)
          }
      })

setMethod("plot", signature(x = "fciAlgo"),
          function(x, y, main = NULL, ...)
      {
          check.Rgraphviz()

          if(is.null(main))
	      main <- deparse(x@call)
	  else ## see also below
	      warning("main title cannot *not* be set yet [Rgraphviz::plot() deficiency]")
          amat <- x@amat
          g <- as(amat,"graphNEL")
          nn <- nodes(g)
          p <- numNodes(g)
          n.edges <- numEdges(g)
          ah.list <- at.list <- rep("none",n.edges)
          counter <- 0
          list.names <- NULL
          amat[amat==1] <- "odot"
          amat[amat==2] <- "normal"
          amat[amat==3] <- "none"
          for (i in seq_len(p-1)) {
              for (j in (i+1):p) {
                  x <- nn[i]
                  y <- nn[j]
                  if (amat[x,y]!=0) {
                      counter <- counter + 1
                      ah.list[[counter]] <- amat[x,y]
                      at.list[[counter]] <- amat[y,x]
                      list.names <- c(list.names,paste(x,"~",y,sep=""))
                  }
              }
          }
          names(ah.list) <- names(at.list) <- list.names
	  edgeRenderInfo(g) <- list(arrowhead= ah.list,
				    arrowtail= at.list)
	  ## Rgraphviz::plot(g, main = main, ...)
          ## XXX undid change by MM, since edge marks didn't work anymore
          ## "known bug in Rgraphviz, but not something they may fix soon"
          Rgraphviz::renderGraph(Rgraphviz::layoutGraph(g))
      })


plotSG <- function(graphObj, y, dist, amat = NA, directed = TRUE,
                   main = paste("Subgraph of ", deparse(substitute(graphObj)),
                   "\nfrom ", y, " with distance ", dist ))
{
  ## Title:
  ## Plot (directed) subgraph
  ## ----------------------------------------------------------------------
  ## Description:
  ## This function plots a subgraph for a specified starting node and a
  ## given graph.
  ## ----------------------------------------------------------------------
  ## Usage:
  ## plotSG( graphObj, y, dist, amat = NA, directed = TRUE )
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## graphObj: Graph object
  ##        y: Starting node
  ##     dist: Distance of nodes included in subgraph from starting node y
  ##     amat: Adjacency matrix of skeleton graph (optional)
  ## directed: Boolean; should the plotted subgraph be directed?
  ## ----------------------------------------------------------------------
  ## Details:
  ## Commencing at the starting point 'y' the function looks for the
  ## neighbouring nodes. Beginning with direct parents and children it
  ## will continue hierarchically through the distances to 'y'. If
  ## 'directed' is set to TRUE the orientation of the edges is taken from
  ## the initial graph.
  ## ----------------------------------------------------------------------
  ## Author: Daniel Stekhoven, Date: 26 Jan 2010, 08:56
  ## ----------------------------------------------------------------------

  check.Rgraphviz()
  stopifnot(dist >= 1)

  ## Extract adjacency matrix (if necessary)
  if (any( is.na(amat)))
    amat <- wgtMatrix(graphObj)

  ## Diagonalise (has no effect if already diagonal)
  amat[amat != 0] <- 1
  amat <- amat + t(amat)
  diag(amat) <- 0 # can that happen anyway??
  amat[amat == 2] <- 1

  ## Find connected nodes hierarchically
  nh <- which( amat[y,] == 1 )
  rel.nodes <- c( y, nh )
  if (dist > 1) {
    for (i in seq_len(dist-1L)) {
      nh <- { if ( length(nh) == 1 )
		  which( amat[nh,] == 1 )
	      else if (length(nh) != 0)
		  which( amat[nh,] == 1, arr.ind = TRUE )[,"col"]
	      ## else NULL
       }
      rel.nodes <- unique( c( rel.nodes, nh ) )
    }
  }

  ## Name nodes
  if(is(graphObj, "graphNEL"))
    names(rel.nodes) <- graphObj@nodes[rel.nodes]

  ## subgraph - distinguish between directed edges or not
  sg <- { if (directed)
	      subGraph(as.character(rel.nodes), graphObj) else
	  as( amat[rel.nodes, rel.nodes], "graphNEL") }
  ## Plot subgraph
  Rgraphviz::plot( sg )
  if(!is.null(main) && !is.na(main))
      title(main=main)
  invisible(sg)
}

triple2numb <- function(p,i,j,k)
{
  ## Purpose: transform a triple i-j-k into a unique number
  ## ----------------------------------------------------------------------
  ## Arguments:-p:number of nodes;-triple: i-j-k
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 11 May 2010, 15:35
  k + p*(j + p*i)
}

find.upd <- function (path, a=NULL, pag, verbose=FALSE, unfVect=NULL, p)
{
  stopifnot((n <- length(path)) >= 2)
  if (n > nrow(pag)) return(list(res = FALSE, uncov.path = NA))
  ## else  --- "real work"

  uncov.path <- NA
  res <- FALSE
  if (n <= nrow(pag)) {
    final <- path[n]
    mittle <- path[n - 1]
    if (n == 2) {
      if (path[1] == path[n]) {
        uncov.path <- path
        res <- TRUE
      }
      else {
        if ((pag[path[1], path[n]] == 1 || pag[path[1], path[n]] == 2) &&
            (pag[path[n], path[1]] == 1 || pag[path[n], path[1]] == 3)) {
          ## check uncovered
          if (pag[a,final]==0 && pag[final,a]==0 && final!=a) {
            ## normal version
            if (length(unfVect)==0) {
              uncov.path <- path
              res <- TRUE
            }
            ## conservative version
            else {
              ## check faithfulness of a-mittle-final
              if (!any(unfVect == triple2numb(p,a,mittle,final)) &&
                  !any(unfVect == triple2numb(p,final,mittle,a))) {
                uncov.path <- path
                res <- TRUE
              }
            }
          }
        }
      }
    }
    else { ## n != 2
      indB <- which((pag[mittle, ] == 1 | pag[mittle, ] == 2) &
                    (pag[, mittle] == 1 | pag[, mittle] == 3))
      indB <- setdiff(indB, a)
      if (length(indB) > 0) {
        counter <- 0
        while (!res && counter < length(indB)) {
          counter <- counter + 1
          b <- indB[counter]
          if (!(b %in% path)) {
            new.path <- c(path[1:(n - 1)], b, path[n])
            check.uncov <- 0
            for (l in 1:(n - 1)) {
              if (pag[new.path[l], new.path[l + 2]] != 0 ||
                  pag[new.path[l + 2], new.path[l]] != 0)

                check.uncov <- check.uncov + 1
            }
            if (check.uncov == 0) {
              if ((pag[b, final] == 1 | pag[b, final] == 2) &&
                  (pag[final, b] == 1 | pag[final, b] == 3)) {
                ## normal version
                if (length(unfVect)==0) {
                  uncov.path <- new.path
                  res <- TRUE
                }
                ## conservative version
                else {
                  ## check faithfulness
                  check.faith <- 0
                  for (l in 1:(n - 1)) {
                    if (any(unfVect == triple2numb(p,new.path[l  ],new.path[l+1],new.path[l+2])) ||
                        any(unfVect == triple2numb(p,new.path[l+2],new.path[l+1],new.path[l])))

                      check.faith <- check.faith + 1
                  }
                }
              }
              else {
		tmp <- find.upd(path = new.path,
				pag=pag, verbose=verbose, unfVect=unfVect, p=p)
                uncov.path <- tmp[[2]]
                res <- tmp[[1]]
              }
            }
            else {
	      tmp <- find.upd(path = new.path,
			      pag=pag, verbose=verbose, unfVect=unfVect, p=p)
              uncov.path <- tmp[[2]]
              res <- tmp[[1]]
            }
          } ## if (b not in path)
        }
      } ## if length(.) > 0
    }
  }
  list(res = res, uncov.path = uncov.path)
}## {find.upd}

pc.cons.intern <- function(sk, suffStat, indepTest, alpha, verbose=FALSE, version.unf=c(NA,NA))
{
  ## Purpose:  For any unshielded triple A-B-C, consider all subsets D of
  ## the neighbors of A and of the neighbors of C, and record the sets
  ## D for which A and C are conditionally independent given D. If B
  ## is in none of these sets, do nothing (it is a
  ## v-structure). If B is in all sets, do nothing (it is not a
  ## v-structure). If B is in some but not all sets, mark the triple as
  ## "unfaithful".
  ## ----------------------------------------------------------------------
  ## Arguments: sk: output returned by function "skeleton"
  ## suffStat: Sufficient statistics for independent tests
  ## indepTest: Function for independence test
  ## alpha: Significance level of test
  ## version.unf[1]: Consider the case, where a indep c given S in the skeleton; furthermore, suppose that a
  ## and c are dependent given every subset of the neigbors in the conservative step; then, some error must
  ## have occured in the pc, but not necessarily in the fci (since here sepsets can contain nodes outside
  ## the neighborhood). If this option is 1 the corresponding triple is marked 'faithful'; if this option
  ## is 2, the corresponding triple is marked unfaithful.
  ## version.unf[2]: 1 do not consider the initial sepset (same as in Tetrad),
  ##         2 also consider the initial sepset
  ## ----------------------------------------------------------------------
  ## Value: unfTripl: Triple that were marked as unfaithful
  ## vers: vector containing the version (1 or 2) of the
  ## corresponding triple saved in unfTripl (1=normal unfaithful triple
  ## that is, B is in some sepsets; 2=triple coming from version.unf[1]==2
  ## that is a and c are indep given the initial sepset but there doesn't
  ## exist a subset of the neighbours that d-separates them)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 12 Feb 2010, 10:43
  ## Modification: Diego Colombo, Date: 6 May 2010, 9:13

  skelObj <- list(G = as(sk@graph,"matrix"), sepset = sk@sepset)
  g <- skelObj$G
  stopifnot(all(g==t(g))) ## g is guaranteed to be symmetric
  p <- dim(g)[1]
  unfTripl <- vers <- rep(NA,min(p*p,100000))
  if (sum(skelObj$G)>0) {
    ind <- which(g==1,arr.ind=TRUE)
    tripleMatrix <- NULL
    ## Go through all edges
    for (i in 1:dim(ind)[1]) {
      a <- ind[i,1]
      b <- ind[i,2]
      allC <- setdiff(which(g[b,]==1),a) ## a-b-c
      newC <- allC[g[a,allC]==0]
      tmpMatrix <- cbind(rep(a,length(newC)),rep(b,length(newC)),newC)
      tripleMatrix <- rbind(tripleMatrix,tmpMatrix)
      colnames(tripleMatrix) <- c("","","")
    }
    if (dim(tripleMatrix)[1]>=1) {
      ## make sure that a < c
      deleteDupl <- rep(0,dim(tripleMatrix)[1])
      for (i in 1:dim(tripleMatrix)[1]) {
        if (tripleMatrix[i,1]>tripleMatrix[i,3]) {
          deleteDupl[i] <- 1
        }
      }
      deletedupl <- which(deleteDupl==1)
      if (length(deletedupl)>0) {
        finalMatrix <- tripleMatrix[-deletedupl,, drop=FALSE]
      }
      else {
        finalMatrix <- tripleMatrix
      }

      counter <- 0
      if (dim(finalMatrix)[1]>=1) {
        for (i in 1:dim(finalMatrix)[1]) {
          ## pay attention to the size of counter
          if (counter==(length(unfTripl)-1)) {
            tmp.vec <- rep(NA,min(p*p,100000))
            unfTripl <- c(unfTripl,tmp.vec)
            vers <- c(vers,tmp.vec)
          }
          a <- finalMatrix[i,1]
          b <- finalMatrix[i,2]
          c <- finalMatrix[i,3]
          nbrsA <- which(g[,a]!=0) ## G symm; c no nbr of a
          nbrsC <- which(g[,c]!=0)
          if (verbose) {
            cat("Sepset by skelet:", skelObj$sepset[[a]][[c]],"and",skelObj$sepset[[c]][[a]],"\n")
          }
	  resTriple <- checkTriple(a,b,c, nbrsA,nbrsC,
				   skelObj$sepset[[a]][[c]],
				   skelObj$sepset[[c]][[a]],
				   suffStat,indepTest,alpha,verbose=verbose,version.unf=version.unf)
          ## 1: in NO set; 2: in ALL sets; 3: in SOME but not all
          ## Take action only if case "3"
          if (resTriple$decision == 3) {
            ## record unfaithful triple
            counter <- counter + 1
            unfTripl[counter] <- triple2numb(p,a,b,c)
            vers[counter] <- resTriple$version
            if (verbose) {
              cat("new unfTriple:", a,b,c, "\n")
            }
          }
          ## can happen the case in Tetrad, so we must save the triple as unfaithful
          ## a and c independent given S but not given subsets of the adj(a) or adj(c)
          if ((version.unf[1]==2) && (resTriple$version==2) && (resTriple$decision!=3)){
            counter <- counter + 1
            unfTripl[counter] <- triple2numb(p,a,b,c)
            vers[counter] <- resTriple$version
            if (verbose) {
              cat("new unfTriple:", a,b,c, "\n")
            }
          }
        }
      }
    }
  }
  length(unfTripl) <- length(vers) <- counter
  list(unfTripl = unfTripl, vers=vers)
} ## function

checkTriple <- function(a,b,c,nbrsA,nbrsC,sepsetA,sepsetC,suffStat,indepTest,alpha,
                        verbose=FALSE,version.unf=c(NA,NA))
{
  ## Purpose: For each subset of nbrsA and nbrsC where a and c are cond. independent,
  ## it is checked if b is in the conditioning set.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## a,b,c: Nodes (positions in adjacency matrix)
  ## nbrsA: Neighbors of a
  ## nbrsC: Neighbors of c
  ## sepsetA, sepsetsC: sepsets of a and c
  ## suffStat: Sufficient statistics for independent tests
  ## indepTest: Function for independence test
  ## alpha: Significance level of test
  ## version.unf[1]: 1 it checks if b is in some sepsets, 2 it also checks if
  ## there exists a sepset which is a subset of the neighbours.
  ## version.unf[2]: 1 same as Tetrad (do not consider the initial sepset),
  ##         2 consider if b is in sepsetA or sepsetC
  ## ----------------------------------------------------------------------
  ## Value: res, vers
  ## res = 1: b is in NO sepset (-> collider)
  ## res = 2: b is in ALL sepsets (-> non-collider)
  ## res = 3: b is in SOME but not all sepsets (-> unfaithful triplet)
  ## vers: version (1 or 2) of the unfaithful triple (1=normal unfaithful triple
  ## that is b is in some sepsets; 2=triple coming from version.unf[1]==2
  ## that is a and c are indep given the initial sepset but there doesn't
  ## exist a subset of the neighbours that d-separates them)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 12 Feb 2010, 12:13
  ## Modification: Diego Colombo, Date: 23 Apr 2010, 11:48

  if ((length(nbrsA) == 0) && (length(nbrsC)==0)) {
    res <- 1
    vers <- 0
  }
  else {
    ## make matrix that encodes all subsets of parents
    tmpLA <- vector("list",length(nbrsA))
    tmpLC <- vector("list",length(nbrsC))
    for (i in seq_along(nbrsA)) {
      tmpLA[[i]] <- c(0,1)
    } ## for (i in seq_along(nbrsA))
    allCombA <- expand.grid(tmpLA) ## all subsets of nbrsA
    for (i in seq_along(nbrsC)) {
      tmpLC[[i]] <- c(0,1)
    } ## for (i in seq_along(nbrsC))
    allCombC <- expand.grid(tmpLC) ## all subsets of nbrsC
    ## loop through all subsets of parents
    indep.true <- NULL
    ## Tetrad
    if (version.unf[2]==1) {
      tmp <- NULL
    }
    ## our version
    else {
      ## check whether b was in original sepset
      tmp <- (b %in% sepsetA | b %in% sepsetC)
    }
    vers <- 0
    ## loop through all subsets of parents of A
    for (i in 1:nrow(allCombA)) {
      S <- nbrsA[which(allCombA[i,]!=0)]
      pval <- indepTest(a, c, S, suffStat)
      if (verbose) {
        cat("(a,b,c)=", a, b, c," S =",S," - pval =",pval,"\n")
      }
      if (pval >= alpha) {
        indep.true <- c(indep.true,TRUE)
        ## is b in set?
        tmp <- c(tmp,(b %in% S))
        vers <- 1
      } ## if (pval >= alpha)
      if (length(tmp)>=2) {
        if ((!all(tmp)) && (sum(tmp)>0)) {
          ## result should be 3
          break
        } ## if ((!all(tmp)) && (sum(tmp)>0))
      } ## if (length(tmp)>1)
    } ## for (i in 1:nrow(allComb))

    ## loop through all subsets of parents of C
    for (i in 1:nrow(allCombC)) {
      S <- nbrsC[which(allCombC[i,]!=0)]
      pval <- indepTest(a, c, S, suffStat)
      if (verbose) {
        cat("(a,b,c)=", a, b, c," S =",S," - pval =",pval,"\n")
      }
      if (pval >= alpha) {
        indep.true <- c(indep.true,TRUE)
        ## is b in set?
        tmp <- c(tmp,(b %in% S))
        vers <- 1
      } ## if (pval >= alpha)
      if (length(tmp)>=2) {
        if ((!all(tmp)) && (sum(tmp)>0)) {
          ## result should be 3
          break
        } ## if ((!all(tmp)) && (sum(tmp)>0))
      } ## if (length(tmp)>1)
    } ## for (i in 1:nrow(allCombC))

    if (version.unf[1]==2 && length(indep.true)==0) {
      vers <- 2
    }
    if (is.null(tmp)) tmp <- FALSE
    res <-
      if (all(tmp))
        2 ## in ALL sets
      else if (all(!tmp))
        1 ## in NO set
      else
        3 ## in SOME sets
  }
  if (verbose) cat("res=", res, "vers=", vers, "\n")
  list(decision=res, version=vers)
}

faith.check <- function(cpath, unfVect, p)
{
  ## Purpose: check if every triple on the 'cpath' is faithful
  ## ----------------------------------------------------------------------
  ## Arguments: cpath: circle path to check for faithfulness
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 25 May 2010, 13:57
  n <- length(cpath)
  for (j in seq_len(n)) {
    if (any(unfVect == triple2numb(p,cpath[j%%n],    cpath[(j+1)%%n],cpath[(j+2)%%n])) ||
	any(unfVect == triple2numb(p,cpath[(j+2)%%n],cpath[(j+1)%%n],cpath[j%%n])))

      return(FALSE)
  }
  TRUE
}

discr.path <- function (path, pag, gInput, verbose = FALSE)
{
  stopifnot((n <- length(path)) >= 3)
  if (n > nrow(pag)) return(pag)
  ## else  --- "update" pag and return it in the end

  pag <- pag
  c <- path[1]
  b <- path[2]
  a <- path[3]
  first.pos <- path[n]
  del.pos <- path[n - 1]
  indD <- which(pag[first.pos, ] != 0 &
                pag[, first.pos] == 2)
  indD <- setdiff(indD, del.pos)
  for (d in indD)  if(pag[c, b] == 1 && all(d != path)) {
    if (pag[c, d] == 0 && pag[d, c] == 0) {
      if ((b %in% gInput@sepset[[c]][[d]]) ||
          (b %in% gInput@sepset[[d]][[c]])) {
        if (verbose)
          cat("\nRule 4: There is a discriminating path between:",
              d, "and", c, "for", b, "and", b, "is in Sepset of",
              c, "and", d, ":", b, "->", c, "\n")
        pag[b, c] <- 2
        pag[c, b] <- 3
      }
      else {
        if (verbose) {
          cat("\nRule 4: There is a discriminating path between:",
              d, "and", c, "for", b, "and", b, "is not in Sepset of",
              c, "and", d, ":", a, "<->", b, "<->", c, "\n")
        }
        pag[a, b] <- pag[b, c] <- pag[c, b] <- 2
      }
    }
    else {
      if (pag[first.pos, d] == 2 &&
          pag[d, c] == 2 && pag[c, d] == 3) {
        pag <- discr.path(path = c(path, d),
                          pag = pag, gInput = gInput, verbose = verbose)
      } ## else : keep 'pag'
    }
  } ## for( d )

  pag
}## {discr.path}

##################################################
## New
##################################################
combn2 <- function(n,k)
{
  ## Purpose: generate all possible combinations of the set n
  ##          of length k without duplicates
  ## ----------------------------------------------------------------------
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 17 Aug 2010, 14:51
  if (length(n)==1) {
    if (k==1) {
      res <- matrix(n,1,1)
      return(res)
    }
  }
  else {
    return(combn(n,k))
  }
}


skeleton.rfci <- function (suffStat, indepTest, p, alpha, verbose = FALSE, fixedGaps = NULL, 
    fixedEdges = NULL, NAdelete = TRUE, m.max = Inf) 
{
  ## Purpose: find an initial skeleton conditioning on the adjacency sets
  ##          It is the same idea as the pc but we compute all the adjacency
  ##          sets at the same time to avoid ordering issues
  ## ----------------------------------------------------------------------
  ## Arguments: same as skeleton
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 23 Nov 2010, 10:14
  
  cl <- match.call()
  pval <- NULL
  sepset <- vector("list", p)
  n.edgetests <- numeric(1)
  if (is.null(fixedGaps)) {
    G <- matrix(rep(TRUE, p * p), nrow = p, ncol = p)
    diag(G) <- FALSE
  }
  else {
    if (!(identical(dim(fixedGaps), c(p, p)))) {
      stop("Dimensions of the dataset and fixedGaps do not agree.")
    }
    else {
      if (fixedGaps != t(fixedGaps)) {
        stop("fixedGaps must be symmetric")
      }
      G <- !fixedGaps
    }
  }
  if (is.null(fixedEdges)) {
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else {
    if (!(identical(dim(fixedEdges), c(p, p)))) {
      stop("Dimensions of the dataset and fixedEdges do not agree.")
    }
    if (fixedEdges != t(fixedEdges)) {
      stop("fixedEdges must be symmetric")
    }
  }
  seq_p <- 1:p
  for (iList in 1:p) sepset[[iList]] <- vector("list", p)
  pMax <- matrix(rep(-Inf, p * p), nrow = p, ncol = p)
  diag(pMax) <- 1
  done <- FALSE
  ord <- 0
  while (!done && any(G) && ord <= m.max) {
    n.edgetests[ord + 1] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)
    ind <- ind[order(ind[, 1]), ]
    remainingEdgeTests <- nrow(ind)
    if (verbose) 
      cat("Order=", ord, "; remaining edges:", remainingEdgeTests, 
          "\n", sep = "")
    ## compute the adjacency sets for any vertex
    nbrsBool.tmp <- vector("list",p) 
    for (i in 1:p) {
      nbrsBool.tmp[[i]] <- G[, i]
    }
    for (i in 1:remainingEdgeTests) {
      if (verbose) {
        if (i%%100 == 0) 
          cat("|i=", i, "|iMax=", nrow(ind), "\n")
      }
      x <- ind[i, 1]
      y <- ind[i, 2]
      if (G[y, x] && (!fixedEdges[y, x])) {
        nbrsBool <- nbrsBool.tmp[[x]]
        nbrsBool[y] <- FALSE
        nbrs <- seq_p[nbrsBool]
        length_nbrs <- length(nbrs)
        if (length_nbrs >= ord) {
          if (length_nbrs > ord) 
            done <- FALSE
          S <- seq(length = ord)
          repeat {
            n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
            pval <- indepTest(x, y, nbrs[S], suffStat)
            if (verbose) 
              cat("x=", x, " y=", y, " S=", nbrs[S], 
                  ": pval =", pval, "\n")
            if (is.na(pval)) 
              pval <- ifelse(NAdelete, 1, 0)
            if (pval > pMax[x, y]) 
              pMax[x, y] <- pval
            if (pval >= alpha) {
              G[x, y] <- G[y, x] <- FALSE
              sepset[[x]][[y]] <- nbrs[S]
              break
            }
            else {
              nextSet <- getNextSet(length_nbrs, ord, S)
              if (nextSet$wasLast) 
                break
              S <- nextSet$nextSet
            }
          }
        }
      }
    }
    ord <- ord + 1
  }
  for (i in 1:(p - 1)) {
    for (j in 2:p) {
      pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j, 
                                                       i])
    }
  }
  if (sum(G) == 0) {
    Gobject <- new("graphNEL", nodes = as.character(1:p))
  }
  else {
    colnames(G) <- rownames(G) <- as.character(1:p)
    Gobject <- as(G, "graphNEL")
  }
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0), 
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests, 
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}


rfci <- function (suffStat, indepTest, p, alpha, verbose = FALSE, fixedGaps = NULL, 
    fixedEdges = NULL, NAdelete = TRUE, m.max = Inf) 
{
  cl <- match.call()
  skel <- skeleton.rfci(suffStat, indepTest, p, alpha, verbose=verbose, 
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges, NAdelete = NAdelete, 
                   m.max = m.max)
  g <- as(skel@graph, "matrix")
  sepset <- skel@sepset
  ##the list of all ordered unshielded triples (the graph g does not change it is just a search!)
  tmp <- find.unsh.triple(g,p)
  listM <- tmp$unshTripl
  vectM <- tmp$unshVect
 
    tripleList <- NULL

  ##check and orient v-structures recursively
  tmp1 <- rfci.vstructures(suffStat, indepTest, p, alpha, sepset, g, listM, vectM, unfVect=tripleList, verbose=verbose)
  graph <- tmp1$graph
  sepset <- tmp1$sepset
  ##orient as many edge marks as possible
  res <- udag2apag(suffStat, indepTest, alpha, graph, sepset, verbose=verbose, unfVect=tripleList)
 
  max.ordSKEL <- skel@max.ord  
  max.ordPD <- 0
  n.edgetestsSKEL <- skel@n.edgetests
  n.edgetestsPD <- 0
  pMax <- skel@pMax
  allPdsep <- vector("list", p)
  sepset <- res$sepset
  Amat <- res$graph
  rfciRes <- new("fciAlgo", amat = Amat, call = cl, n = integer(0), 
        max.ord = as.integer(max.ordSKEL), max.ordPDSEP = as.integer(max.ordPD), 
        n.edgetests = n.edgetestsSKEL, n.edgetestsPDSEP = n.edgetestsPD, 
        sepset = sepset, pMax = pMax, allPdsep = allPdsep)
  rfciRes
}

find.unsh.triple <- function(g,p)
{
  ## Purpose: find the ordered (<x,y,x> with x<z) list of all the unshielded
  ##          triples in the graph
  ## ----------------------------------------------------------------------
  ## Arguments: g: adjacency matrix
  ## ----------------------------------------------------------------------
  ## Values: - unshTripl: matrix with 3 rows containing in each column
  ##                    an unshielded triple
  ##         - unshVect: containing the unique number for each column in unshTripl
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 21 Oct 2010, 14:05

  if (any(g!=0)) {
    unshtmp <- NULL
    ##find all unshielded triple in g
    indS <- which(g == 1, arr.ind = TRUE) ##x-y
    for (i in 1:dim(indS)[1]) {
      x <- indS[i, 1]
      y <- indS[i, 2]
      allZ <- setdiff(which(g[y, ] == 1), x) ##x-y-z
      if (length(allZ) > 0) {
        for (j in 1:length(allZ)) {
          z <- allZ[j]
          if (g[x,z]==0 & g[z,x]==0) {
            ##save the matrix
            unshtmp <- cbind(unshtmp,c(x,y,z))
          }
        }
      }
    }
    ##delete duplicates in the matrix
    if (dim(unshtmp)[2] > 0) {
      deleteDupl <- rep(0,dim(unshtmp)[2])
      for (i in 1:dim(unshtmp)[2]) {
        if (unshtmp[1,i]>unshtmp[3,i]) {
          deleteDupl[i] <- 1
        }
      }
      deletedupl <- which(deleteDupl==1)
      if (length(deletedupl)>0) {
        unshTripl <- unshtmp[,-deletedupl, drop=FALSE]
      }
      else {
        unshTripl <- unshtmp
      }
    }
    ##define the vector with the unique number for each triple
    if (dim(unshTripl)[2] > 0) {
      unshVect <- rep(NA, dim(unshTripl)[2])
      for (k in 1 : dim(unshTripl)[2]) {
        unshVect[k] <- triple2numb(p, unshTripl[1,k], unshTripl[2,k], unshTripl[3,k])
      }
    }
    else {
      unshVect <- NULL
    }
  }
  list(unshTripl=unshTripl, unshVect=unshVect)
}
  

rfci.vstructures <- function(suffStat, indepTest, p, alpha, sepset, graph, unshTripl, unshVect, unfVect=NULL, verbose=FALSE)
{
  ## Purpose: check the unshielded triples in unshTripl as v-structures
  ##          recursively check if new unshielded triples have been found
  ##          then save and orient the final ones in finalList and finalVect
  ## ----------------------------------------------------------------------
  ## Arguments: - suffStat, indepTest, p ,alpha: arguments from the algorithm
  ##            - skel: output of the function skeleton
  ##            - unshTripl, unshVect: list/numbers of the unshielded triples
  ##                                   in graph
  ##            - unfVect: list of unfaithful triples for the conservative version
  ## ----------------------------------------------------------------------
  ## Values: - updated sepset and graph
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 21 Oct 2010, 14:15

  ##check every column (unshielded triple) of unshTripl
  ##NOTE: what's done is done!!!! I don't look back in the matrix
  if (dim(unshTripl)[2] > 0) {
    if (verbose) {
      cat("\nCheck any unshielded triple for conditional dependence","\n")
    }
    ##per default every triple is defined as a v-structure
    trueVstruct <- rep(TRUE,dim(unshTripl)[2])
    ##check all the columns in unshTripl as v-structure or not
    checktmp <- dep.triple(suffStat, indepTest, p, alpha, sepset, graph, unshTripl, unshVect, trueVstruct, verbose=verbose)
    ##save the updated objects
    graph <- checktmp$graph
    sepset <- checktmp$sepset
    ##note that no column is deleted from the matrix, we can add new triples instead
    unshTripl <- checktmp$triple
    unshVect <- checktmp$vect
    trueVstruct <- checktmp$trueVstruct
    ##orient the v-structures
    ##if there is at least one triple with the desired properties
    if (any(trueVstruct)){
      if (verbose) {
        cat("\nOrient the v-structures","\n")
      }
      for (i in 1:dim(unshTripl)[2]) {
        if (trueVstruct[i]) {
          x <- unshTripl[1, i]
          y <- unshTripl[2, i]
          z <- unshTripl[3, i]
          if ((graph[x, z] == 0) & !((y %in% sepset[[x]][[z]]) | (y %in% sepset[[z]][[x]]))) {
            ##normal version
            if (length(unfVect)==0) {
              ##this is to avoid the problem of:
              ##assume that <a,b,c> was saved in finalList because a "true" unshielded triple
              ##but after in the matrix appears <b,c,d> and the edge b-c is deleted
              ##of course the triple <a,b,c> stays in the finalList but we cannot orient the edge
              ##b-c because it doesn'e exist anymore
              if (graph[x,y]!=0 & graph[z,y]!=0) {
                graph[x, y] <- 2
                graph[z, y] <- 2
                if (verbose) {
                  cat("\n", x, "*->", y, "<-*", z, "\n")
                  cat("Sxz=", sepset[[z]][[x]], "and", 
                      "Szx=", sepset[[x]][[z]], "\n")
                }
              }
            }
            ##conservative version
            else {
              ##check if x-y-z is faithful
              if (!any(unfVect==triple2numb(p,x,y,z)) & !any(unfVect==triple2numb(p,z,y,x))) {
                if (graph[x,y]!=0 & graph[z,y]!=0) {
                  graph[x, y] <- 2
                  graph[z, y] <- 2
                  if (verbose) {
                    cat("\n", x, "*->", y, "<-*", z, "\n")
                    cat("Sxz=", sepset[[z]][[x]], "and", 
                        "Szx=", sepset[[x]][[z]], "\n")
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  list(sepset = sepset, graph = graph)
}


dep.triple <- function(suffStat, indepTest, p, alpha, sepset, apag, unshTripl, unshVect, trueVstruct, verbose=verbose)
{
  ## Purpose: test the two edges of any unshielded triple in unshTripl (column)
  ##          for dependence given SepSet. If independent find the minimal sepset,
  ##          delete the edge, define the triple as FALSE in trueVstruct and search
  ##          in the graph for new triple or destroyed ones. Otherwise do nothing. 
  ## ----------------------------------------------------------------------
  ## Arguments: - suffStat, indepTest, p, alpha, sepset, apag: skeleton
  ##              parameters
  ##            - unshTripl: matrix containing the unshielded triples (columns)
  ##            - unshVect: triple2numb of unshTripl
  ##            - trueVstruct: vector containing T/F for the v-structures
  ##            - k: recursive parameter for the matrix unshTripl
  ## ----------------------------------------------------------------------
  ## Values: - updated sepset, apag, unshTripl, unshVect, trueVstruct
  ## Author: Diego Colombo, Date: 16 Aug 2010, 15:07

  k <- 0
  while (k < dim(unshTripl)[2]) {
    k <- k + 1
    if (trueVstruct[k]) {
      ##triple under inspection x-y-z
      x <- unshTripl[1, k]
      y <- unshTripl[2, k]
      z <- unshTripl[3, k]
      SepSet <- unique(c(sepset[[x]][[z]],sepset[[z]][[x]]))
      SepSet <- setdiff(SepSet,y)
      if (verbose) {
        cat("\nTriple=","<",x,",",y,",",z,">", ";SepSet=", SepSet, "\n")
      }
      if (length(SepSet)!=0) {
        del1 <- FALSE
        ##check x-y
        ##check first x and y given the whole sepset(x,z)
        pvalue1.whole <- indepTest(x, y, SepSet, suffStat)
        if (pvalue1.whole >= alpha) {
          ##afterwards we have to define this triple as FALSE in trueVstruct
          del1 <- TRUE
          ## x and y are independent, then find the minimal sepset
          done <- FALSE
          ord <- 0
          while (!done && ord < length(SepSet)) {
            ord <- ord + 1
            ##all combinations of SepSet of size ord
            tmp.combn <- combn2(SepSet, ord)
            for (i in 1:dim(tmp.combn)[2]) {
              pval <- indepTest(x, y, tmp.combn[,i], suffStat)
              if (verbose) {
                cat("x=", x, " y=", y, "S=", tmp.combn[,i], ": pval =", pval, "\n")
              }
              if (pval >= alpha) {
                ##delete edge and save set in sepset
                if (verbose) {
                  cat("Independence found: delete edge between",x,"and",y,"\n")
                }
                apag[x, y] <- apag[y, x] <- 0
                sepset[[x]][[y]] <- sepset[[y]][[x]] <- tmp.combn[,i]
                done <- TRUE
                break
              }
            }
          }
      
          ##case 1: before we had a triangle and now it is an unshielded triple x-m-y
          indM <- which((apag[x, ] == 1 & apag[, x] == 1) & (apag[y, ] == 1 & apag[, y] == 1))
          indM <- setdiff(indM,c(x,y,z))##just to be sure
          if (length(indM) > 0) {
            for (j in 1:length(indM)) {
              m <- indM[j]
              ##in the matrix the first column is always smaller than the third
              if (x < y) {
                unshTripl <- cbind(unshTripl, c(x,m,y))
                unshVect <- c(unshVect, triple2numb(p, x, m, y))
              }
              else {
                unshTripl <- cbind(unshTripl, c(y,m,x))
                unshVect <- c(unshVect, triple2numb(p, y, m, x))
              }
              ##per default this new triple is set as TRUE in trueVstruct
              trueVstruct <- c(trueVstruct,TRUE)
            }
          }
          ##case 2: an existent unshielded triple has been destroyed
          ##case 2.a): we had q-x-y or y-x-q in the matrix but not anymore
          indQ <- which((apag[x, ] == 1 & apag[, x] == 1) & (apag[y, ] == 0 & apag[, y] == 0))
          indQ <- setdiff(indQ,c(x,y,z))##just to be sure
          if (length(indQ) > 0) {
            for (j in 1:length(indQ)) {
              q <- indQ[j]
              ##define the triple as FALSE in trueVstruct
              if (q < y) {
                ##first check if it exists
                delTripl <- unshVect == triple2numb(p, q, x, y)  
                if (any(delTripl)) {
                  ##if we haven't checked the triple yet, define it as FALSE
                  ##if (k < which.max(delTripl)) {
                    trueVstruct[which.max(delTripl)] <- FALSE
                  ##}
                }
              }
              else {
                delTripl <- unshVect == triple2numb(p, y, x, q)
                if (any(delTripl)) {
                  ##if (k < which.max(delTripl)) {
                    trueVstruct[which.max(delTripl)] <- FALSE
                  ##}
                }
              }
            }
          }
          ##case 2.b): we had r-y-x or x-y-r in the matrix but not anymore
          indR <- which((apag[x, ] == 0 & apag[, x] == 0) & (apag[y, ] == 1 & apag[, y] == 1))
          indR <- setdiff(indQ,c(x,y,z))##just to be sure
          if (length(indR) > 0) {
            for (j in 1:length(indR)) {
              r <- indR[j]
              ##define the triple as FALSE in trueVstruct
              if (r < x) {
                ##first check if it exists
                delTripl <- unshVect == triple2numb(p, r, y, x)
                ##if we haven't checked the triple yet, save it as FALSE
                if (any(delTripl)) {
                  ##if (k < which.max(delTripl)) {
                    trueVstruct[which.max(delTripl)] <- FALSE
                  ##}
                }
              }
              else {
                delTripl <- unshVect == triple2numb(p, x, y, r)
                if (any(delTripl)) {
                  ##if (k < which.max(delTripl)) {
                    trueVstruct[which.max(delTripl)] <- FALSE
                  ##}
                }
              }
            }
          }
        }
    
        ##check z-y
        del2 <- FALSE
        ##check first z and y given the whole sepset(x,z)
        pvalue2.whole <- indepTest(z, y, SepSet, suffStat)
        if (pvalue2.whole >= alpha) {
          ##afterwards we have to set this triple as FALSE in trueVstruct
          del2 <- TRUE
          ## z and y are independent, then find the minimal sepset
          Done <- FALSE
          Ord <- 0
          while (!Done && Ord < length(SepSet)) {
            Ord <- Ord + 1
            ##all combinations of SepSet of size ord
            tmp.combn <- combn2(SepSet, Ord)
            for (i in 1:dim(tmp.combn)[2]) {
              pval <- indepTest(z, y, tmp.combn[,i], suffStat)
              if (verbose) {
                cat("x=", z, " y=", y, " S=", tmp.combn[,i], ": pval =", pval, "\n")
              }
              if (pval >= alpha) {
                ##delete edge and save set in sepset
                if (verbose) {
                  cat("Independence found: delete edge between",y,"and", z,"\n")
                }
                apag[z, y] <- apag[y, z] <- 0
                sepset[[z]][[y]] <- sepset[[y]][[z]] <- tmp.combn[,i]
                Done <- TRUE
                break
              }
            }
          }
      
          ##case 1: before we had a triangle and now it is an unshielded triple z-m-y
          indM <- which((apag[z, ] == 1 & apag[, z] == 1) & (apag[y, ] == 1 & apag[, y] == 1))
          indM <- setdiff(indM,c(x,y,z))##just to be sure
          if (length(indM) > 0) {
            for (j in 1:length(indM)) {
              m <- indM[j]
              ##in the matrix the first column is always smaller than the third
              if (z < y) {
                unshTripl <- cbind(unshTripl, c(z,m,y))
                unshVect <- c(unshVect, triple2numb(p, z, m, y))
              }
              else {
                unshTripl <- cbind(unshTripl, c(y,m,z))
                unshVect <- c(unshVect, triple2numb(p, y, m, z))
              }
              ##per default the triple is defined as TRUE
              trueVstruct <- c(trueVstruct,TRUE)
            }
          }
          ##case 2: an existent unshielded triple has been destroyed
          ##case 2.a): we had q-z-y or y-z-q in the matrix but not anymore
          indQ <- which((apag[z, ] == 1 & apag[, z] == 1) & (apag[y, ] == 0 & apag[, y] == 0))
          indQ <- setdiff(indQ,c(x,y,z))##just to be sure
          if (length(indQ) > 0) {
            for (j in 1:length(indQ)) {
              q <- indQ[j]
              ##define the triple as FALSE in trueVstruct
              if (q < y) {
                ##first check if it exists
                delTripl <- unshVect == triple2numb(p, q, z, y)
                ##if we haven't checked the triple yet, save it as FALSE
                if (any(delTripl)) {
                  ##if (k < which.max(delTripl)) {
                    trueVstruct[which.max(delTripl)] <- FALSE
                  ##}
                }
              }
              else {
                delTripl <- unshVect == triple2numb(p, y, z, q)
                if (any(delTripl)) {
                  ##if (k < which.max(delTripl)) {
                    trueVstruct[which.max(delTripl)] <- FALSE
                  ##}
                }
              }
            }
          }
          ##case 2.b): we had r-y-z or z-y-r in the matrix but not anymore
          indR <- which((apag[z, ] == 0 & apag[, z] == 0) & (apag[y, ] == 1 & apag[, y] == 1))
          indR <- setdiff(indQ,c(x,y,z))##just to be sure
          if (length(indR) > 0) {
            for (j in 1:length(indR)) {
              r <- indR[j]
              ##define the triple as FALSE in trueVstruct
              if (r < z) {
                ##first check if it exists
                delTripl <- unshVect == triple2numb(p, r, y, z)
                ##if we haven't checked the triple yet, save it as FALSE
                if (any(delTripl)) {
                  ##if (k < which.max(delTripl)) {
                    trueVstruct[which.max(delTripl)] <- FALSE
                  ##}
                }
              }
              else {
                delTripl <- unshVect == triple2numb(p, z, y, r)
                if (any(delTripl)) {
                  ##if (k < which.max(delTripl)) {
                    trueVstruct[which.max(delTripl)] <- FALSE
                  ##}
                }
              }
            }
          }
        }
        ##if at least one edge has been deleted this is not a future v-structure
        if (any(del1,del2)) {
          trueVstruct[k] <- FALSE
        }
      }
      ##the sepset is empty
      ##so surely <x,y,z> is an unshielded triple because they cannot be indep given the empty set
      ##nothing changed in sepset and in graph
      else{
        sepset <- sepset
        apag <- apag
      }
      ##recursion on the next column of unshTripl
      ##rec.res <- dep.triple(suffStat, indepTest, p, alpha, sepset, apag, unshTripl, unshVect, trueVstruct, k, verbose=verbose)
      ##save the modified objects
      ##unshTripl <- rec.res$triple
      ##unshVect <- rec.res$vect
      ##sepset <- rec.res$sepset
      ##apag <- rec.res$graph
      ##trueVstruct <- rec.res$trueVstruct
    }
  }
  list(triple = unshTripl, vect = unshVect, sepset = sepset, graph = apag, trueVstruct = trueVstruct)
}



udag2apag <- function (suffStat, indepTest, alpha, apag, sepset, rules = rep(TRUE, 10), verbose = FALSE, unfVect=NULL) 
{
  ## Purpose: use the 10 orientation rules to orient the skeleton.
  ##          Note that the preliminaries v-structures are already oriented
  ##          in apag
  ## ----------------------------------------------------------------------
  ## Arguments: outputs of the function rfci.vstructures
  ## ----------------------------------------------------------------------
  ## Values: updated apag (oriented) and sepset
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 21 Oct 2010, 15:13

  if (any(apag!=0)) {
    if (verbose) {
      cat("\nOrient the edge marks using the 10 orientation rules","\n")
    }
    p <- dim(apag)[1]
    old_apag1 <- matrix(rep(0, p^2), nrow = p, ncol = p)
    while (sum(!(old_apag1 == apag)) > 0) {
      old_apag1 <- apag
      if (rules[1]) {
        ind <- which((apag == 2 & t(apag) != 0), arr.ind = TRUE)
        if (length(ind) > 0) {
          for (i in 1:dim(ind)[1]) {
            a <- ind[i, 1]
            b <- ind[i, 2]
            indC <- which((apag[b, ] != 0 & apag[, b] == 
                           1) & (apag[a, ] == 0 & apag[, a] == 0))
            indC <- setdiff(indC, a)
            if (length(indC) > 0) {
              ##normal version
              if (length(unfVect)==0) {
                apag[b, indC] <- 2
                apag[indC, b] <- 3
                if (verbose) {
                  cat("\nRule 1:", a, "*->", b, "o-*", indC, 
                      " where ", a, " and ", indC, " not connected: ", b, "->", indC, "\n")
                }
              }
              ##conservative
              else {
                for (j in 1:length(indC)) {
                  c <- indC[j]
                  ##check that a-b-c faithful
                  if (!any(unfVect==triple2numb(p,a,b,c)) & !any(unfVect==triple2numb(p,c,b,a))) {
                    apag[b, c] <- 2
                    apag[c, b] <- 3
                    if (verbose) {
                      cat("\nRule 1':", a, "*->", b, "o-*", c, 
                          " where ", a, " and ", c, " not connected: ",a, b, c, "faithful triple", b, "->", c, "\n")
                    }
                  }
                }
              }
            }
          }
        }
      }
      if (rules[2]) {
        ind <- which((apag == 1 & t(apag) != 0), arr.ind = TRUE)
        if (length(ind) > 0) {
          for (i in 1:dim(ind)[1]) {
            a <- ind[i, 1]
            c <- ind[i, 2]
            indB <- which(((apag[a, ] == 2 & apag[, a] == 
                            3) & (apag[c, ] != 0 & apag[, c] == 2)) | 
                          ((apag[a, ] == 2 & apag[, a] != 0) & (apag[c, 
                                                 ] == 3 & apag[, c] == 2)))
            if (length(indB) > 0) {
              apag[a, c] <- 2
              if (verbose) {
                cat("\nRule 2:", a, "->", indB, "*->", 
                    c, "or", a, "*->", indB, "->", c, "and", 
                    a, "*-o", c, ":", a, "*->", c, "\n")
              }
            }
          }
        }
      }
      if (rules[3]) {
        ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)
        if (length(ind) > 0) {
          for (i in 1:dim(ind)[1]) {
            b <- ind[i, 1]
            d <- ind[i, 2]
            indAC <- which((apag[b, ] != 0 & apag[, b] == 
                            2) & (apag[, d] == 1 & apag[d, ] != 0))
            if (length(indAC) >= 2) {
              ##normal version
              if (length(unfVect)==0) {
                counter <- 0
                while ((counter < (length(indAC) - 1)) & 
                       (apag[d, b] != 2)) {
                  counter <- counter + 1
                  ii <- counter
                  while ((ii < length(indAC)) & (apag[d, 
                                                     b] != 2)) {
                    ii <- ii + 1
                    if (apag[indAC[counter], indAC[ii]] == 
                        0 & apag[indAC[ii], indAC[counter]] == 
                        0) {
                      apag[d, b] <- 2
                      if (verbose) {
                        cat("\nRule 3:", d, "*->", b, "\n")
                      }
                    }
                  }
                }
              }
              ##conservative version
              else {
                comb.indAC <- combn2(indAC,2)
                for (j in 1:dim(comb.indAC)[2]) {
                  a <- comb.indAC[1,j]
                  c <- comb.indAC[2,j]
                  if (apag[a,c]==0 & apag[c,a]==0 & c!=a) {
                    ##check fatihfulness a-d-c
                    if (!any(unfVect==triple2numb(p,a,d,c)) & !any(unfVect==triple2numb(p,c,d,a))) {
                      apag[d, b] <- 2
                      if (verbose) {
                        cat("\nRule 3':", a, d, c, "faithful",  d, "*->", b, "\n")
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      if (rules[4]) {
        ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)## b o-* c
        while (length(ind) > 0) {
          b <- ind[1, 1]
          c <- ind[1, 2]
          ind <- ind[-1,, drop=FALSE]
          ##find all a s.t. a -> c and a <-* b
          indA <- which((apag[b, ] == 2 & apag[, b] != 0) & (apag[c, ] == 3 & apag[, c] == 2))
          ##chose one a s.t. the initial triangle structure exists and the edge hasn't been oriented yet
          while ((length(indA) > 0) & (apag[c,b] == 1)) {
            a <- indA[1]
            indA <- indA[-1]
            ##path is the initial triangle
            path <- c(a, b, c)
            ##Done is TRUE if either we found a minimal path or no path exists for this triangle
            Done <- FALSE
            while ((!Done) & (apag[a,b] != 0) & (apag[a,c] != 0) & (apag[b,c] != 0)) {
              ##find a minimal disciminating path for a,b,c
              tmp.path <- min.discr.path(p, pag = apag, path = path, verbose = verbose)
              length.path <- length(tmp.path)
              ##if the path doesn't exists, we are done with this triangle
              if (length.path == 1) {
                Done <- TRUE
              }
              else {
                ##a path exists and needs to be checked and maybe oriented
                ##first check every single edge for independence
                tmp.checked <- CheckEdges(suffStat, indepTest, p, alpha, apag, sepset, tmp.path)
                ##save updated graph and sepset
                sepset <- tmp.checked$sepset
                apag <- tmp.checked$apag
                ##no edge deleted ----> orient the edges
                if (!tmp.checked$deleted) {
                  ##if b is in sepset
                  if ((b %in% sepset[[tmp.path[1]]][[tmp.path[length.path]]]) | (b %in% sepset[[tmp.path[length.path]]][[tmp.path[1]]])) {
                    if (verbose) {
                      cat("\nNew Rule 4: There is a discriminating path between:", 
                          tmp.path[1], "and", c, "for", b, "and", b, "is in Sepset of", 
                          c, "and", tmp.path[1], ":", b, "->", c, "\n")
                    }
                    apag[b, c] <- 2
                    apag[c, b] <- 3
                  }
                  else {
                    ##if b is not in sepset
                    if (verbose) {
                      cat("\nNew Rule 4: There is a discriminating path between:", 
                          tmp.path[1], "and", c, "for", b, "and", b, "is not in Sepset of", 
                          c, "and", tmp.path[1], ":", a, "<->", b, "<->", 
                          c, "\n")
                    }
                    apag[a, b] <- apag[b, c] <- apag[c, b] <- 2
                  }
                  Done <- TRUE
                }
              }
            }
          }
        }
      }
      if (rules[5]) {
        ind <- which((apag == 1 & t(apag) == 1), arr.ind = TRUE)
        if (length(ind) > 0) {
          for (i in 1:dim(ind)[1]) {
            a <- ind[i, 1]
            b <- ind[i, 2]
            indC <- which((apag[a, ] == 1 & apag[, a] == 
                           1) & (apag[b, ] == 0 & apag[, b] == 0))
            indC <- setdiff(indC, b)
            indD <- which((apag[b, ] == 1 & apag[, b] == 
                           1) & (apag[a, ] == 0 & apag[, a] == 0))
            indD <- setdiff(indD, a)
            if (length(indC) > 0 & length(indD) > 0) {
              for (j in 1:length(indC)) {
                c <- indC[j]
                for (l in 1:length(indD)) {
                  d <- indD[l]
                  if (apag[c, d] == 1 & apag[d, c] == 1) {
                    ##normal version
                    if (length(unfVect)==0) {
                      apag[a, b] <- apag[b, a] <- 3
                      apag[a, c] <- apag[c, a] <- 3
                      apag[c, d] <- apag[d, c] <- 3
                      apag[d, b] <- apag[b, d] <- 3
                      if (verbose) {
                        cat("\nRule 5: There exists an uncovered circle path between", 
                            a, "and", b, ":", a, "-", b, 
                            "and", a, "-", c, "-", d, "-", 
                            b, "\n")
                      }
                    }
                    ##conservative
                    else {
                      ##check that every triple on the circle is faithful
                      path2check <- c(a,c,d,b)
                      faithres <- faith.check(path2check, unfVect, p)
                      if (faithres==0) {
                        apag[a, b] <- apag[b, a] <- 3
                        apag[a, c] <- apag[c, a] <- 3
                        apag[c, d] <- apag[d, c] <- 3
                        apag[d, b] <- apag[b, d] <- 3
                        if (verbose) {
                          cat("\nRule 5': There exists a faithful uncovered circle path between", 
                              a, "and", b, ":", a, "-", b, 
                              "and", a, "-", c, "-", d, "-", 
                              b, "\n")
                        }
                      }
                    }
                  }
                  else {
                    path <- c(a, c, d, b)
                    apag <- ucp(path = path, pag = apag, 
                                verbose = verbose, unfVect=unfVect, p=p)
                  }
                }
              }
            }
          }
        }
      }
      if (rules[6]) {
        ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)
        if (length(ind) > 0) {
          for (i in 1:dim(ind)[1]) {
            b <- ind[i, 1]
            c <- ind[i, 2]
            indA <- which(apag[b, ] == 3 & apag[, b] == 
                          3)
            if (length(indA) > 0) {
              apag[c, b] <- 3
              if (verbose) {
                cat("\nRule 6:", a, "-", b, "o-*", c, ":", 
                    b, "-*", c, "\n")
              }
            }
          }
        }
      }
      if (rules[7]) {
        ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)
        if (length(ind) > 0) {
          for (i in 1:dim(ind)[1]) {
            b <- ind[i, 1]
            c <- ind[i, 2]
            indA <- which((apag[b, ] == 3 & apag[, b] == 
                           1) & (apag[c, ] == 0 & apag[, c] == 0))
            indA <- setdiff(indA, c)
            if (length(indA) > 0) {
              ##normal version
              if (length(unfVect)==0) {
                apag[c, b] <- 3
                if (verbose) {
                  cat("\nRule 7:", indA, "-o", b, "o-*", 
                      c, "and", indA, " and", c, " are not adjacent:", 
                      b, "-*", c, "\n")
                }
              }
              ##conservative
              else {
                for (j in 1:length(indA)) {
                  a <- indA[j]
                  ##check fatihfulness of a-b-c
                  if (!any(unfVect==triple2numb(p,a,b,c)) & !any(unfVect==triple2numb(p,c,b,a))) {
                    apag[c, b] <- 3
                    if (verbose) {
                      cat("\nRule 7':", a, "-o", b, "o-*", 
                          c, "and", a, " and", c, " are not adjacent and", a,b,c, "faithful:",  
                          b, "-*", c, "\n")
                    }
                  }
                }
              }
            }
          }
        }
      }
      if (rules[8]) {
        ind <- which((apag == 2 & t(apag) == 1), arr.ind = TRUE)
        if (length(ind) > 0) {
          for (i in 1:dim(ind)[1]) {
            a <- ind[i, 1]
            c <- ind[i, 2]
            indB <- which(((apag[a, ] == 2 & apag[, a] == 
                            3) | (apag[a, ] == 1 & apag[, a] == 3)) & 
                          (apag[c, ] == 3 & apag[, c] == 2))
            if (length(indB) > 0) {
              apag[c, a] <- 3
              if (verbose) {
                cat("\nRule 8:", a, "->", indB, "->", c, 
                    "or", a, "-o", indB, "->", c, "and", 
                    a, "o->", c, ":", a, "->", c, "\n")
              }
            }
          }
        }
      }
      if (rules[9]) {
        ind <- which((apag == 2 & t(apag) == 1), arr.ind = TRUE)
        if (length(ind) > 0) {
          for (i in 1:dim(ind)[1]) {
            a <- ind[i, 1]
            c <- ind[i, 2]
            indB <- which((apag[a, ] == 2 | apag[a, ] == 
                           1) & (apag[, a] == 1 | apag[, a] == 3) & 
                          (apag[c, ] == 0 & apag[, c] == 0))
            indB <- setdiff(indB, c)
            if (length(indB) > 0) {
              for (k in 1:length(indB)) {
                b <- indB[k]
                indD <- which((apag[b, ] == 2 | apag[b, 
                                     ] == 1) & (apag[, b] == 1 | apag[, b] == 
                                                      3) & (apag[a, ] == 0 & apag[, a] == 0))
                indD <- setdiff(indD, a)
                if (length(indD) > 0) {
                  for (l in 1:length(indD)) {
                    d <- indD[l]
                    if (apag[c, a] != 3) {
                      if ((apag[c, d] == 1 | apag[c, d] == 
                           3) & (apag[d, c] == 1 | apag[d, 
                                       c] == 2)) {
                        ##normal version
                        if (length(unfVect)==0) {
                          apag[c, a] <- 3
                          if (verbose) {
                            cat("\nRule 9: There exists an upd between", 
                                a, "and", c, ":", a, " ->", 
                                c, "\n")
                          }
                        }
                        ##conservative version
                        else {
                          path2check <- c(a,b,d,c)
                          faithres <- faith.check(path2check, unfVect, p)
                          if (faithres==0) {
                            apag[c, a] <- 3
                            if (verbose) {
                              cat("\nRule 9': There exists a faithful upd between", 
                                  a, "and", c, ":", a, " ->", 
                                  c, "\n")
                            }
                          }
                        }
                      }
                      else {
                        path <- c(c, a, b, d)
                        apag <- upd(path = path, pag = apag, 
                                    verbose = verbose, unfVect=unfVect, p=p)
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      if (rules[10]) {
        ind <- which((apag == 2 & t(apag) == 1), arr.ind = TRUE)
        if (length(ind) > 0) {
          for (i in 1:dim(ind)[1]) {
            a <- ind[i, 1]
            c <- ind[i, 2]
            indB <- which((apag[c, ] == 3 & apag[, c] == 
                           2))
            if (length(indB) >= 2) {
              for (j in 1:length(indB)) {
                b <- indB[j]
                indD <- setdiff(indB, b)
                if (length(indD) > 0 & apag[c, a] != 3) {
                  for (k in 1:length(indD)) {
                    d <- indD[k]
                    if ((apag[a, b] == 1 | apag[a, b] == 
                         2) & (apag[b, a] == 1 | apag[b, a] == 
                               3) & (apag[a, d] == 1 | apag[a, d] == 
                                     2) & (apag[d, a] == 1 | apag[d, a] == 
                                           3) & (apag[d, b] == 0 & apag[b, d] == 
                                                 0)) {
                      ##normal version
                      if (length(unfVect)==0) {
                        apag[c, a] <- 3
                        if (verbose) {
                          cat("\nRule 10 with mu = beta = ", 
                              b, "and omega = theta =", d, 
                              ":", a, "->", c, "\n")
                        }
                      }
                      ##conservative version
                      else {
                        ##check faithfulness of b-a-d
                        if (!any(unfVect==triple2numb(p,b,a,d)) & !any(unfVect==triple2numb(p,d,a,b))) {
                          apag[c, a] <- 3
                          if (verbose) {
                            cat("\nRule 10' with mu = beta = ", 
                                b, "and omega = theta =", d, 
                                "and",b,a,d,"faithful:", a, "->", c, "\n")
                          }
                        }
                      }
                    }
                    else {
                      indA <- which((apag[a, ] == 1 | 
                                     apag[a, ] == 2) & (apag[, a] == 
                                           1 | apag[, a] == 3), arr.ind = TRUE)
                      indA <- setdiff(indA, c)
                      if (length(indA >= 2)) {
                        for (l in 1:length(indA)) {
                          first.pos <- indA[l]
                          indAA <- setdiff(indA, first.pos)
                          if ((length(indAA) > 0) & (apag[c, 
                                                          a] != 3)) {
                            for (s in 1:length(indAA)) {
                              sec.pos <- indAA[s]
                              p1 <- find.upd(path = c(first.pos, 
                                               b), a = a, pag = apag, 
                                             verbose = verbose, unfVect=unfVect, p=p)
                              p2 <- find.upd(path = c(sec.pos, 
                                               d), a = a, pag = apag, 
                                             verbose = verbose, unfVect=unfVect, p=p)
                              if (p1$res == TRUE & p2$res == TRUE) {
                                mu <- p1$uncov.path[1]
                                omega <- p2$uncov.path[1]
                                if ((mu != omega) & (apag[mu, 
                                       omega] == 0) & (apag[omega, 
                                                           mu] == 0)) {
                                  ##normal version
                                  if (length(unfVect)==0) {
                                    apag[c, a] <- 3
                                    if (verbose) {
                                      cat("\nRule 10:", a, 
                                          "->", c, "\n")
                                    }
                                  }
                                  ##conservative version
                                  else {
                                    if (!any(unfVect==triple2numb(p,mu,a,omega)) & !any(unfVect==triple2numb(p,omega,a,mu))) {
                                      apag[c, a] <- 3
                                      if (verbose) {
                                        cat("\nRule 10':",mu,a,omega,"faithful:", a, 
                                            "->", c, "\n")
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }## for k
                }
              }
            }
          }
        }
      }
    }
  }

  list(graph = apag, sepset = sepset)
}

CheckEdges <- function(suffStat, indepTest, p, alpha, apag, sepset, path, unfVect=NULL, verbose=FALSE)
{
  ## Purpose: check if every edge on the path should exist in R4 
  ## ----------------------------------------------------------------------
  ## Values: - updated sepset and apag
  ##         - checked==FALSE no edge has been deleted on the path
  ##                  ==TRUE the discriminating path doesn't exist anymore
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 17 Aug 2010, 16:21
  
  n.path <- length(path)
  ##did we delete an edge?
  done <- FALSE
  ##define the sepset
  SepSet.tot <- unique(c(sepset[[path[1]]][[path[n.path]]],sepset[[path[n.path]]][[path[1]]]))
  if (length(SepSet.tot)!=0) {
    if (verbose) {
      cat("\nCheck discriminating path:", path, "for dependence of any edge given Sepset=", SepSet.tot,"\n")
    }
    ##check every edge on tha path for independence given every possible subset of SepSet
    for (i in 1 : (n.path-1)) {
      x <- path[i]
      y <- path[i+1]
      SepSet <- setdiff(SepSet.tot,c(x,y))
      if (verbose) {
        cat("\nEdge:",x,"*-*",y,";Sepset=",SepSet,"\n")
      }
      if (length(SepSet)!=0) {
        for (j in 1 : length(SepSet)) {
          ##all combinations of SepSet of size ord
          tmp.combn <- combn2(SepSet, j)
          for (ii in 1:dim(tmp.combn)[2]) {
            pval <- indepTest(x, y, tmp.combn[,ii], suffStat)
            if (verbose) {
              cat("x=", x, " y=", y, " S=", tmp.combn[,ii], ": pval =", pval,"\n")
            }
            if (pval >= alpha) {
              ##delete edge and save set in sepset
              if (verbose) {
                cat("Independence found: delete edge between",x,"and",y,"\n")
              }
              apag[x, y] <- apag[y, x] <- 0
              sepset[[x]][[y]] <- sepset[[y]][[x]] <- tmp.combn[,ii]
              ##before we had a triangle and now it is an unshielded triple x-m-y
              newList <- NULL
              newVect <- NULL
              indM <- which((apag[x, ] != 0 & apag[, x] != 0) & (apag[y, ] != 0 & apag[, y] != 0))
              indM <- setdiff(indM,c(x,y))##just to be sure
              ##create the list with all the new unshielded triples to be tested
              if (length(indM) > 0) {
                for (jj in 1:length(indM)) {
                  m <- indM[jj]
                  if (x < y) {
                    newList <- cbind(newList, c(x,m,y))
                    newVect <- c(newVect, triple2numb(p,x,m,y))
                  }
                  else {
                    newList <- cbind(newList, c(y,m,x))
                    newVect <- c(newVect, triple2numb(p,y,m,x))
                  }
                }
                ##new unshielded triple to be tested
                tmp <- rfci.vstructures(suffStat, indepTest, p, alpha, sepset, apag, newList, newVect, unfVect=unfVect, verbose=verbose)
                ##save the modified graph g in apag
                apag <- tmp$graph
                ##save the new sepset
                sepset <- tmp$sepset
              }              
              done <- TRUE    
            }
          }
        }
      }
    }
  }
  ##if SepSet is the empty set do nothing because surely the vertices are dependent
  list(deleted = done, apag = apag, sepset=sepset)
}


min.discr.path <- function(p, pag = NA, path = NA, verbose = FALSE)
{
  ## Purpose: find a minimal discrimating path for a,b,c saved in path.
  ## If a path exists this is the output, otherwise NA
  ## ----------------------------------------------------------------------
  ## Arguments: p: number of nodes in the graph
  ##            pag: adjacency matrix
  ##            path: a,b,c under interest
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 25 Jan 2011, 10:49

  visited <- rep(FALSE, p)
  a <- path[1]
  b <- path[2]
  c <- path[3]
  visited[a] <- visited[b] <- visited[c] <- TRUE
  ##find all neighbours of a not visited yet
  indD <- which(pag[a,] != 0 & pag[,a] == 2 & !visited) ## d *-> a
  ##queue tmp.Q=next element in the queue
  ##tmp.P is a vector of length p where the node j in position i means that j is the predecessor of i on the path
  tmp.Q <- indD
  tmp.P <- rep(FALSE, p)
  tmp.P[tmp.Q] <- a
  tmp.P[a] <- b
  tmp.P[b] <- c
  done <- FALSE
  min.path <- NA
  while ((length(tmp.Q) > 0) & (!done)) {
    ##next element in the queue
    d <- tmp.Q[1]
    tmp.Q <- tmp.Q[-1]
    visited[d] <- TRUE
    pred <- tmp.P[d]
    if (pag[c,d] == 0 & pag[d,c] == 0) {
      ##discriminating path found
      ##find the whole path using tmp.P starting from d, since a,b,c are the initial structure
      min.path <- d
      counter <- pred
      while (counter != a) {
        min.path <- c(min.path,counter)
        counter <- tmp.P[counter]
      }
      min.path <- c(min.path,path)
      done <- TRUE
    }
    else {
      ##d is connected to c -----> search iteratively
      if (pag[d,c] == 2 & pag[c,d] == 3 & pag[pred,d] == 2) {
        ##find all neighbourd of d not visited yet
        indR <- which(pag[d,] != 0 & pag[,d] == 2 & !visited) ## r *-> d
        ##update the queues
        tmp.Q <- c(tmp.Q, indR)
        tmp.P[indR] <- d
      }
    }
  }
  return(min.path)
}

ucp <- function (path, pag, verbose = FALSE, unfVect=NULL, p)
{
  stopifnot((n <- length(path)) >= 3)
  if (n > nrow(pag)) return(pag)
  ## else  --- "update" pag and return it in the end

  n1 <- n - 1L
  indX <- which((pag[path[n - 2], ] == 1 &
                 pag[, path[n - 2]] == 1))
  for (x in indX) if (!any(x == path) && pag[path[1], path[n]] != 3) {
    new.path <- c(path[c(1:(n - 2))], x, path[c(n1,n)])
    if (pag[x, path[n1]] == 1 &&
        pag[path[n1], x] == 1) {
      check.uncov <- FALSE
      for (j in 1:n1) {
        if (pag[new.path[j], new.path[j + 2]] != 0 ||
            pag[new.path[j + 2], new.path[j]] != 0)
          {
            check.uncov <- TRUE
            break
          }
      }
      if (!check.uncov) {
        ## normal version
        if (length(unfVect)==0) {
          pag[new.path[1], new.path[n + 1]] <- pag[new.path[n + 1], new.path[1]] <- 3
          for (j in 1:n) {
            pag[new.path[j], new.path[j + 1]] <- pag[new.path[j + 1], new.path[j]] <- 3
          }
          if (verbose)
            cat("\nRule 5: There exists an uncovered circle path between",
                new.path[1], "and", new.path[n + 1], ":", new.path[1], "-", new.path[n + 1],
                "and for each edge on the path",  new.path, "\n")
        }
        ## conservative version
        else if (faith.check(new.path, unfVect, p)) {
          pag[new.path[1], new.path[n + 1]] <- pag[new.path[n + 1], new.path[1]] <- 3
          for (j in 1:n) {
            pag[new.path[j], new.path[j + 1]] <- pag[new.path[j + 1], new.path[j]] <- 3
          }
          if (verbose) {
            cat("\nRule 5': There exists a faithful uncovered circle path between",
                new.path[1], "and", new.path[n + 1], ":", new.path[1], "-", new.path[n + 1],
                "and for each edge on the path",  new.path, "\n")
          }
        }
      }
      else if (pag[new.path[1], new.path[3]] == 0 &&
               pag[new.path[3], new.path[1]] == 0) {
        pag <- ucp(path = new.path, pag = pag,
                   verbose=verbose, unfVect=unfVect, p=p)
      }
    }
    else {
      pag <- ucp(path = new.path, pag = pag,
                 verbose=verbose, unfVect=unfVect, p=p)
    }
  } ## for(x ..) if(...)
  pag
} ## {ucp}

## Recursive (!) updating :
upd <- function (path, pag, verbose = FALSE, unfVect=NULL, p)
{
  stopifnot((n <- length(path)) >= 2)
  if (n > nrow(pag)) return(pag)
  ## else  --- "update" pag and return it in the end

  n1 <- n - 1L
  c <- path[1]
  a <- path[2]
  b <- path[n1]
  d <- path[n]
  indX <- which((pag[d, ] == 2 | pag[d, ] == 1) &
                (pag[, d] == 1 | pag[, d] == 3) &
                (pag[b, ] == 0 & pag[, b] == 0))
  for (x in indX) if (!any(x == path) && pag[c, a] != 3) {
    new.path <- c(path[2:n], x, path[1])

    if ((pag[x, c] == 1 || pag[x, c] == 2) &&
        (pag[c, x] == 1 || pag[c, x] == 3)) {

      check.uncov <- FALSE
      for (j in 1:n1) {
        if (pag[new.path[j], new.path[j + 2]] != 0 ||
            pag[new.path[j + 2], new.path[j]] != 0)
          {
            check.uncov <- TRUE
            break
          }
      }
      if (!check.uncov) {
        if (length(unfVect)==0) { ## normal version
          pag[c, a] <- 3
          if (verbose)
            cat("\nRule 9: There exists an upd between",
                new.path, ":", a, "->", c, "\n")
        }
        else if (faith.check(new.path, unfVect, p)) { ## conservative version
          pag[c, a] <- 3
          if (verbose)
            cat("\nRule 9': There exists a faithful upd between",
                new.path, ":", a, "->", c, "\n")
        }
      }
      else {
        pag <- upd(path = c(path, x), pag = pag,
                   verbose = verbose, unfVect=unfVect, p=p)
      }
    }
    else {
      pag <- upd(path = c(path, x), pag = pag,
                 verbose = verbose, unfVect=unfVect, p=p)
    }

  } ## for(x ..) if(...)
  pag
}## {upd}
