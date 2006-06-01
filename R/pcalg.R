### -*- mode: R; kept-new-versions: 30; kept-old-versions: 20; ess-indent-level: 2 -*-

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

rmvDAG <- function(n, dag, errDist = c("normal", "cauchy", "mix", "t4"),
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
              

pcAlgo <- function(dm, alpha, corMethod = "standard", verbose = FALSE) {
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
  new("pcAlgo",
      graph = Gobject,
      call = cl, n = n, max.ord = as.integer(ord-1),
      n.edgetests = n.edgetests, sepset = sepset,
      zMin = zMin)
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
                  c("standard", "Qn", "ogkScaleTau2", "ogkQn"))
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


##- glmCIformula <- function(x,y,s,v)
##- {
##-   ## Purpose: Generate a formula to be used in a GLM for testing CI(x,y|s)
##-   ##          using given variable names. The response variable is called
##-   ##          "Freq"
##-   ## ----------------------------------------------------------------------
##-   ## Arguments: - v: vector containing the variable names of the data frame
##-   ##            - x,y: integers specifying the variables of the data frame
##-   ##            - s: vector containing the indices of the variables
##-   ## ----------------------------------------------------------------------
##-   ## Value: - ci: Formula Object for fitting the conditional independence model
##-   ##        - dep: Formula Object for fitting the dependence model
##-   ## ----------------------------------------------------------------------
##-   ## Author: Markus Kalisch, Date: 23 Mar 2006, 12:21
##-   c1 <- v[x]
##-   c2 <- v[y]
##-   d <- paste(v[x],"*",v[y],sep="")
##-   for (i in 1:length(s)) {
##-     c1 <- paste(c1,"*",v[s[i]])
##-     c2 <- paste(c2,"*",v[s[i]])
##-     d <- paste(d,"*",v[s[i]])
##-   }
##-   ci <- as.formula(paste("Freq~",c1,"+",c2,sep=""))
##-   dep <- as.formula(paste("Freq~",d,sep=""))
##-   list(ci=ci,dep=dep)
##- }
##-


