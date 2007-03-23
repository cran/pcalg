
randomDAG <- function(n, prob, lB = 0.1, uB = 1) {
  ## Purpose: Randomly generate a DAG (graph object as in graph-package).
  ## The resulting graph is topologically ordered from low to high node
  ## numbers.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - n: Number of nodes
  ## - prob: Probability of connecting a node to another node with higher
  ##         topological ordering
  ## - lB, uB: Lower and upper limit of weights between connected nodes
  ##   (chosen uniformly at random)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006, 15:59

  stopifnot(n >= 2,
            is.numeric(prob), length(prob) == 1, 0 <= prob, prob <= 1,
            is.numeric(lB), is.numeric(uB), lB <= uB)
  ## create DAG
  V <- as.character(1:n)
  edL <- as.list(V)
  names(edL) <- V
  for (i in seq(length = n-2)) {
    listSize <- rbinom(1, n-i, prob)
    edgeList <- sample(seq(i+1,n), size = listSize)
    ## print(list(i=i,eseq=seq(i+1,n),listSize=listSize,elist=edgeList))
    weightList <- runif(length(edgeList), min = lB, max = uB)
    edL[[i]] <- list(edges = edgeList, weights = weightList)
  }
  if (rbinom(1,1, prob) == 1) ## then generate a last edge
    edL[[n-1]] <- list(edges = n,
                       weights = runif(1, min = lB, max = uB))

  new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
}

## A version of this is also in	/u/maechler/R/MM/Pkg-ex/graph/weightmatrix.R

wgtMatrix <- function(g)
{
    ## Purpose: work around "graph" package's  as(g, "matrix") bug
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

    if(!"weight" %in% (edNms <- names(edgeDataDefaults(g))))
        edgeDataDefaults(g, "weight") <- 1:1
    ew <- edgeData(g, attr = "weight")
    w <- unlist(ew)
    ## we need the *transposed* matrix anyway:
    tm <- t(as(g, "matrix"))
    ## now is a 0/1 - matrix (instead of 0/wgts) with the 'graph' bug
    if(any(w != 1)) ## fix it
        tm[tm != 0] <- w
    tm
}

rmvDAG <- function(n, dag, errDist = c("normal", "cauchy", "mix", "mixt3", "mixN100","t4"), mix = 0.1, errMat = NULL)
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

## Result of pcAlgo(), i.e. PC-Algorithm :
setClass("pcAlgo", representation =
	 list(graph = "graph",
	      call = "call",
	      n	   = "integer",
	      max.ord = "integer",
	      n.edgetests= "numeric",
              sepset= "list",
              zMin= "matrix"))
              

pcAlgo <- function(dm, alpha, corMethod = "standard", verbose = FALSE, directed=FALSE) {
  ## Purpose: Perform PC-Algorithm, i.e., estimate skeleton of DAG given data
  ## Output is an unoriented graph object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler

  stopifnot((n <- nrow(dm)) >= 1,
            (p <- ncol(dm)) >= 1)
  cl <- match.call()
  sepset <- vector("list",p)
  zMin <- matrix(rep(Inf,p*p),nrow=p,ncol=p)
  diag(zMin) <- rep(0,p)
  for (iList in 1:p) sepset[[iList]] <- vector("list",p)
  C <- mcor(dm, method = corMethod)
  cutoff <- qnorm(1 - alpha/2)
  n.edgetests <- numeric(1)# final length = max { ord}
  ## G := complete graph :
  G <- matrix(rep(TRUE,p*p), nrow = p, ncol = p)
  diag(G) <- FALSE
  seq_p <- 1:p

  done <- FALSE
  ord <- 0
  while (!done && any(G)) {
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
          S <- seq(length = ord)

          ## now includes special cases (ord == 0) or (length_nbrs == 1):
          repeat { ## condition w.r.to all  nbrs[S] of size 'ord' :
          ##  if (condIndFisherZ(x,y, nbrs[S], C,n, cutoff,verbose)) {
## MM: want to use this: --- but it changes the result in some cases!!
##            cat("X=",x,"|Y=",y,"|ord=",ord,"|nbrs=",nbrs[S],"|iMax=",nrow(ind),"\n")
            n.edgetests[ord+1] <- n.edgetests[ord+1]+1
            z <- zStat(x,y, nbrs[S], C,n)
            if(abs(z)<zMin[x,y]) zMin[x,y] <- abs(z)
##            cat(paste("x:",x,"y:",y,"S:"),nbrs[S],paste("z:",z,"\n"))
            if (abs(z) <= cutoff) {
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

  if(verbose) { cat("Final graph adjacency matrix:\n"); print(symnum(G)) }

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
      zMin[i,j] <- zMin[j,i] <- min(zMin[i,j],zMin[j,i])
    }
  }
### TODO:
### - for each edge[i,j], save the largest P-value ``seen''
### - for each [i,j],     save (the size of) the neighborhood << as option only!
  res <- new("pcAlgo",
      graph = Gobject,
      call = cl, n = n, max.ord = as.integer(ord-1),
      n.edgetests = n.edgetests, sepset = sepset,
      zMin = zMin)
  if (directed) res <- udag2pdag(res)
  res
}


setMethod("summary", "pcAlgo",
          function(object) {
 	    cat("\nObject of class 'pcAlgo', from Call: \n",
                deparse(object@call),
 		"\n\nFitted skeleton based on ", object@n, "observations:\n")
 	    print(object@graph)
            cat("\nMax. order of algorithm: ",object@max.ord,
                "\nNumber of edgetests from m = 0 up to m =",object@max.ord,
                ": ",object@n.edgetests)
            tmp <- object@graph@edgeL
            nbrs <- rep(0,length(tmp))
            for (i in 1:length(tmp)) {
              nbrs[i] <- length(tmp[[i]]$edges)
            }
            cat("\nMax. number of neighbours: ",max(nbrs),
                "at node(s)",which(nbrs==max(nbrs)),
                "\nAvg. number of neighbours: ",mean(nbrs),"\n")
          })


setMethod("show", "pcAlgo",
	  function(object) {
	    cat("Object of class 'pcAlgo', from Call: \n", deparse(object@call),
		"\n\nFitted skeleton based on ", object@n, "observations:\n")
	    print(object@graph)
	    invisible(object)
	  })

setMethod("plot", signature(x = "pcAlgo"),
	  function(x, y, main = NULL, zvalue.lwd = FALSE, lwd.max = 7, ...)
      {
	if(is.null(main))
	    main <- deparse(x@call)
        if (zvalue.lwd & numEdges(x@graph)!=0) {
          lwd.Matrix <- x@zMin
          lwd.Matrix <- ceiling(lwd.max*lwd.Matrix/max(lwd.Matrix))
          z <- agopen(x@graph,name="lwdGraph")
          eLength <- length(z@AgEdge)
          for (i in 1:eLength) {
            x <- as.numeric(z@AgEdge[[i]]@head)
            y <- as.numeric(z@AgEdge[[i]]@tail)
            z@AgEdge[[i]]@lwd <- lwd.Matrix[x,y]
          }
          plot(z, main = main, ...)
        } else {
          plot(x@graph, main = main, ...)
        }
      })



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
 
  if (length(S) < 4) {
    r <- pcorOrder(x,y, S, C)
  } else {
    k <- solve(C[c(x,y,S),c(x,y,S)])
    r <- -k[1,2]/sqrt(k[1,1]*k[2,2])
    ##      r <- .C("parcorC",as.double(res),as.integer(x-1),as.integer(y-1),as.integer(S-1),as.integer(length(S)),as.integer(dim(C)[1]),as.double(as.vector(C)))[[1]]
  }
  
  res <- sqrt(n- length(S) - 3) * ( 0.5*log( (1+r)/(1-r) ) )
  if (is.na(res)) res <- 0
  
##VERBOSE  cat(" (",x,",",y,") | ",S," : z-Stat = ",res,"\n", sep='')
  
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
  ## - cutoff: Cutoff for Significance level for individual
  ##           partial correlation tests
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006, 17:32

  ## R Variante
  r <- pcorOrder(x,y, S, C)

  ## C Variante
  ##C res <- 0
  ##C p <- dim(C)[1]
  ##C r <- .C("parcor",as.double(res),as.integer(x-1),as.integer(y-1),as.integer(S-1),as.integer(length(S)),as.integer(p),as.double(C))[[1]]

  z <- 0.5*log( (1+r)/(1-r) )
  T <- sqrt(n-length(S)-3)*z
  ## cat(" (",x,",",y,") | ",S," : T = ",T,"\n", sep='')

  if (is.na(T))
      return(TRUE) ## <==>  T <- 0 # if problem, delete edge: be conservative
  ## else
  return(abs(T) <= cutoff) ## TRUE / FALSE
}

pcorOrder <- function(i,j, k, C, cut.at = 0.9999999) {
  ## Purpose: Compute partial correlation
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - i,j,k: Partial correlation of i and j given k
  ## - C: Correlation matrix among nodes
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler

  ord <- length(k)
##-   if(isTRUE(getOption("debug.pcalg"))) { ## "visualize" recursive algorithm
##-     nf <- sys.nframe()
##-     cat(sprintf("%s parCor() {nf %2d }: (i,j; k)= (%d,%d; %s)\n",
##-                 paste(rep(" ", nf), collapse = ''), nf, i,j,
##-                 paste(k, collapse = ",")))
##-   }
  if (ord == 0) {
    r <- C[i,j]
  }
  else if (ord == 1) {
    r <- (C[i,j] - C[i,k]*C[j,k]) / sqrt( (1-C[j,k]^2)*(1-C[i,k]^2) )
  }
  else { ## ord >= 2
    s <- k[ord]
    k <- k[-ord]
    w1 <- pcorOrder(i,j, k, C)
    w2 <- pcorOrder(i,s, k, C)
    w3 <- pcorOrder(j,s, k, C)
    r <- (w1 - w2*w3) / sqrt((1 - w2^2)*(1 - w3^2))
  }
  if (is.na(r)) r=0
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
                  c("standard", "Qn", "QnStable", "ogkScaleTau2", "ogkQn"))
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
	   Qn_mu <- function(x, mu.too = FALSE, ...)
	     c(if(mu.too) median(x), Qn(x, ...))
	   cov2cor(covOGK(dm, n.iter = 2, sigmamu = Qn_mu,
			  weight.fn = hard.rejection)$cov)
	 },
	 "standard" = cor(dm))
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
  ## - amat: Adjacency matrix of DAG
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
  ordered.nodes <- topOrder(amat) ##parents before children
  edge.df <- make.edge.df(amat)
  
  eOrder <- 0
  while(any(is.na(edge.df$order))) {
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
        xRes <- rep(0,length(nbr.nodes))
        for(i in 1:length(nbr.nodes)) {
          x <- nbr.nodes[i]
          ## is edge edge x-y unlabeled?
          xRes[i] <- length(intersect(which(edge.df$xmin==min(node,x) &
                                            edge.df$xmax==max(node,x)),
                                      which(is.na(edge.df$order))))>0
        }
        ## choose unlabeled edge with highest order node
        if(sum(xRes)>0) {
          nbr.unlab <- nbr.nodes[which(xRes==TRUE)] #nbrnodes w. unlabeled edges
          tmp <- ordered.nodes[which(ordered.nodes %in% nbr.unlab)]
          y <- tmp[length(tmp)]
          ## y <- last(ordered.nodes[which(ordered.nodes %in% nbr.unlab)])
          edge.df$order[which(edge.df$xmin==min(node,y)
                              & edge.df$xmax==max(node,y))] <- eOrder
          eOrder <- eOrder+1
          found <- TRUE
        }
      }
    }
  }
  edge.df
}


labelEdges <- function(amat) {
  ## Purpose: Label the edges in a DAG with "compelled" and "reversible"
  ## (for extension to a CPDAG)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: Adjacency matrix of DAG
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:40

  ## label=TRUE -> compelled
  ## label=FALSe -> reversible
  edge.df <- orderEdges(amat)
  edge.df$label <- rep(NA,dim(edge.df)[1])
  edge.df <- edge.df[order(edge.df$order),]

  while(sum(is.na(edge.df$label))>0) {
    x.y <- which(is.na(edge.df$label))[1]
    x <- edge.df$tail[x.y]
    y <- edge.df$head[x.y]
    e1 <- which(edge.df$head==x & edge.df$label==TRUE)
    if(length(e1)>0) {
      i <- 1
      go.on <- TRUE
      while(i<=length(e1) & go.on==TRUE) {
        w <- edge.df$tail[e1[i]]
        if(length(which(edge.df$tail==w & edge.df$head==y))==0) {
          edge.df$label[which(edge.df$head==y)] <- TRUE
          go.on <- FALSE
        }
        else {
          edge.df$label[which(edge.df$head==y & edge.df$tail==w)] <- TRUE
          i <- i+1
        }
      }
    }
    ## edges going to y not starting from x
    cand <- which(edge.df$head==y & edge.df$tail!=x)
    if (length(cand)>0) {
      valid.cand <- rep(FALSE,length(cand))
      for (iz in 1:length(cand)) {
        z <- edge.df$tail[cand[iz]]
        NOT.parent.of.x <- (length(which(edge.df$tail==z & edge.df$head==x))==0)
        if(NOT.parent.of.x) valid.cand[iz] <- TRUE
      }
      cand <- cand[valid.cand]
    }
    if(length(cand)>0) {
      edge.df$label[which(edge.df$head==y & is.na(edge.df$label))] <- TRUE
    } else {
      edge.df$label[which(edge.df$head==y & is.na(edge.df$label))] <- FALSE
    }
  }
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
  ## transform DAG to adjacency matrix
  if (numEdges(dag)==0) {
    cpdag <- matrix(0,p,p)
  } else {
    dag <- as(dag,"matrix")
    dag[dag!=0] <- 1
    
    ## dag is adjacency matrix
    e.df <- labelEdges(dag)
    cpdag <- matrix(rep(0,p*p),nrow=p,ncol=p)
    for (i in 1:dim(e.df)[1]) {
      if (e.df$label[i]) {
        cpdag[e.df$tail[i],e.df$head[i]] <- 1
      } else {
        cpdag[e.df$tail[i],e.df$head[i]] <- cpdag[e.df$head[i],e.df$tail[i]] <- 1
      }
    }
  }
  rownames(cpdag) <- colnames(cpdag) <- as.character(seq(1,p))
  as(cpdag,"graphNEL")
}

##################################################
## pdag2dag
##################################################
find.sink <- function(gm) {
  ## Purpose: Find sink of an adj matrix; return numeric(0) if there is none
  ## a sink my have incident indirected edges, but no directed ones
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: Adjacency matrix
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:28

  ## treat undirected edges
  undir.nbrs <- which((gm==t(gm) & gm==1),arr.ind=TRUE)
  gm[undir.nbrs] <- 0
  undir.nbrs <- unique(as.vector(undir.nbrs))
  ## treat directed edges
  res <- union(which(apply(gm,1,sum)==0),undir.nbrs)
  res
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
    for (i in 1:length(un)) {
      y <- un[i]
      adj.x <- setdiff(which(gm[x,]==1 | gm[,x]==1),y)
      adj.y <- setdiff(which(gm[y,]==1 | gm[,y]==1),x)
      axINay <- all(adj.x %in% adj.y)
      res <- (res & axINay)
    }
  }
  res
}

pdag2dag <- function(g) {
  ## Purpose: Generate a consistent extension of a PDAG to a DAG; if this
  ## is not possible, a random extension of the skeleton is returned and
  ## a warning is issued.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - g: PDAG (graph object)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Sep 2006, 15:21

  if (numEdges(g)==0) {
    res <- g
  } else {
    gm <- as(g,"matrix")
    gm[which(gm>0 & gm!=1)] <- 1
    p <- dim(gm)[1]
    
    gm2 <- gm
    a <- gm
    go.on <- TRUE
    go.on2 <- FALSE
    while(length(a)>1 & sum(a)>0 & go.on) {
      go.on <- FALSE
      go.on2 <- TRUE
      sinks <- find.sink(a)
      if (length(sinks)>0) {
        counter <- 1
        while(counter<=length(sinks) & go.on2==TRUE) {
          x <- sinks[counter]
          if (adj.check(a,x)) {
            go.on2 <- FALSE
            ## orient edges
            inc.to.x <- which(a[,x]==1 & a[x,]==1) ## undirected
            if (length(inc.to.x)>0) {
              real.inc.to.x <- as.numeric(row.names(a)[inc.to.x])
              real.x <- as.numeric(row.names(a)[x])
              gm2[real.inc.to.x,real.x] <- rep(1,length(inc.to.x))
              gm2[real.x,real.inc.to.x] <- rep(0,length(inc.to.x))
            }
            ## remove x and all edges connected to it
            a <- a[-x,-x]
          }
          counter <- counter+1
        }
      }
      go.on <- !go.on2
    }
    if (go.on2==TRUE) {
      res <- as(amat2dag(gm),"graphNEL")
      warning("PDAG not extendible: Random DAG on skeleton drawn")
    } else {
      res <- as(gm2,"graphNEL")
    }
  }
  res
}

amat2dag <- function(amat) {
  ## Purpose: Transform an adjacency matrix to a DAG (graph object)
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
  res.dag <- skel[i.ord,i.ord]
  res.dag
}

##################################################
## udag2pdag
##################################################
udag2pdag <- function(gInput) {
  ## Purpose: Transform the Skeleton of a pcAlgo-object to a PDAG using
  ## the rules of Pearl. The output is again a pcAlgo-object. 
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gInput: pcAlgo object
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Sep 2006, 15:03

  res <- gInput
  if (numEdges(gInput@graph)>0) {
    g <- as(gInput@graph,"matrix")
    p <- dim(g)[1]
    pdag <- g
    ind <- which(g==1,arr.ind=TRUE)
    
    ## Create minimal pattern
    for (i in 1:dim(ind)[1]) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(g[y,]==1),x)
      
      if (length(allZ)>0) {
        for (j in 1:length(allZ)) {
          z <- allZ[j]
          if ((g[x,z]==0) & !((y %in% gInput@sepset[[x]][[z]]) |
                  (y %in% gInput@sepset[[z]][[x]]))) {
            ## cat("\n",x,"->",y,"<-",z,"\n")
            ## cat("Sxz=",gInput@sepset[[z]][[x]],"Szx=",gInput@sepset[[x]][[z]])
            pdag[x,y] <- pdag[z,y] <- 1
            pdag[y,x] <- pdag[y,z] <- 0
          }
        }
      }
    }
    
    ## Test whether this pdag allows a consistent extension
    res2 <- pdag2dag(as(pdag,"graphNEL"))
    
    if (class(res2)=="graphNEL") {
      ## Convert to complete pattern: use rules by Pearl
      old_pdag <- matrix(rep(0,p^2),nrow=p,ncol=p)
      while (sum(!(old_pdag==pdag))>0) {
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
              ## cat("\nRule 1:",a,"->",b," und ",b,"-",indC,": ",b,"->",indC)
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
              ## cat("\nRule 2:",a,"->",b)
            }
          }
        }
        ## x11()
        ## plot(as(pdag,"graphNEL"), main="After Rule2")
        
        ## rule 3
        ind <- which((pdag==1 & t(pdag)==1), arr.ind=TRUE) ## a -> b
        if (length(ind)>0) {
          for (i in 1:dim(ind)[1]) {
            a <- ind[i,1]
            b <- ind[i,2]
            indC <- which( (pdag[a,]==1 & pdag[,a]==1) & (pdag[,b]==1 & pdag[b,]==0))
            g2 <- pdag[indC,indC]
            if (length(g2)==1) {
              g2 <- 0
            } else {
              diag(g2) <- rep(0,length(indC))
            }
            if (any(g2==0)) {
              pdag[a,b] <- 1
              pdag[b,a] <- 0
              ## cat("\nRule 3:",a,"->",b)
            }
          }
        }
        ## x11()
        ## plot(as(pdag,"graphNEL"), main="After Rule3")

        ## rule 4
        ind <- which((pdag==1 & t(pdag)==1), arr.ind=TRUE) ## a - b
        if (length(ind)>0) {
          for (i in 1:dim(ind)[1]) {
            a <- ind[i,1]
            b <- ind[i,2]
            indC <- which( (pdag[a,]==1 & pdag[,a]==1) & (pdag[,b]==0 & pdag[b,]==0))
            l.indC <- length(indC)
            if (l.indC>0) {
              found <- FALSE
              ic <- 0
              while(!found & (ic < l.indC)) {
                ic <- ic + 1
                c <- indC[ic]
                indD <- which( (pdag[c,]==1 & pdag[,c]==0) & (pdag[,b]==1 & pdag[b,]==0))
                if (length(indD)>0) {
                  found <- TRUE
                  pdag[b,a] = 0
                }
              }
            }
          }
        }

      }
      res@graph <- as(pdag,"graphNEL")
    }
  }
  return(res)
}

##################################################
## udag2cpdag
##################################################
udag2cpdag <- function(pc)
{
  ## Purpose: Transform the result of the undirected part of the
  ## PC-algorithm to a CPDAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - pc: pcAlgo-Object returned from undirected part of PC-algoritm
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 13 Sep 2006, 17:16
  res <- pc
  pc.directed <- udag2pdag(pc)
  pc.dag <- pdag2dag(pc.directed@graph)
  ## pc.cpdag <- dag2cpdag(as(pc.dag,"matrix"))
  ## rownames(pc.cpdag) <- colnames(pc.cpdag) <- as.character(seq(1,p))
  ## res@graph <- as(pc.cpdag,"graphNEL")
  res@graph <- dag2cpdag(pc.dag)
  res
}

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
  if (class(g1)=="pcAlgo") g1 <- g1@graph
  if (class(g2)=="pcAlgo") g2 <- g2@graph
  
  if (class(g1)=="graphNEL") {
    m1 <- t(wgtMatrix(g1))
    m1[m1 != 0] <- rep(1, sum(m1 != 0))
  }
  if (class(g2)=="graphNEL") {
    m2 <- t(wgtMatrix(g2))
    m2[m2 != 0] <- rep(1, sum(m2 != 0))
  }

  p <- dim(m1)[2]
  shd <- 0
  
  # Remove superfluous edges from g1
  s1 <- m1 + t(m1)
  s2 <- m2 + t(m2)
  s1[which(s1==2)] <- 1
  s2[which(s2==2)] <- 1
  ds <- s1-s2
  inds <- which(ds>0)
  m1[inds] <- 0
  shd <- shd + sum(ds>0)/2

  # Add missing edges to g1
  ind <- which(ds<0)
  m1[ind] <- m2[ind]
  shd <- shd + length(ind)/2

  # Compare Orientation
  d <- abs(m1-m2)
  shd <- shd + sum((d + t(d))>0)/2

  shd
}
