trueCov <- function(g) {
  if (class(g)=="graphNEL") {
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
        ecov[j,i] = sum(w[i,pa[[i]]]*ecov[j,pa[[i]]])
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
                 w[i,pa[[i]],drop=FALSE]%*%xcov[pa[[i]],pa[[j]]]%*%t(w[j,pa[[j]],drop=FALSE]) + ecov[j,i]
             })
    }
  }
  xcov
}

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
  if (rbinom(1,1, prob) == 1) {## then generate a last edge
    edL[[n-1]] <- list(edges = n,
                       weights = runif(1, min = lB, max = uB))
  } else {
    edL[[n-1]] <- list(edges = integer(0), weights = numeric(0))
  }
  edL[[n]] <- list(edges = integer(0), weights = numeric(0))

  new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
}

## A version of this is also in	/u/maechler/R/MM/Pkg-ex/graph/weightmatrix.R

wgtMatrix <- function(g)
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

    if(!"weight" %in% (edNms <- names(edgeDataDefaults(g))))
        edgeDataDefaults(g, "weight") <- 1:1
    ew <- edgeData(g, attr = "weight")
    w <- unlist(ew)
    ## we need the *transposed* matrix anyway:
    tm <- t(as(g, "matrix"))
    ## now is a 0/1 - matrix (instead of 0/wgts) with the 'graph' bug
    if(any(w != 1)) ## fix it
        tm[tm != 0] <- w
    ## tm_[i,j]==1 iff i->j
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
              

pcSelect <- function(y,dm, alpha, corMethod = "standard", verbose = 0, directed=FALSE) {
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
  ## backward compatibility
  if (verbose==FALSE) verbose <- 0
  if (verbose==TRUE) verbose <- 1

  stopifnot((n <- nrow(dm)) >= 1,
            (p <- ncol(dm)) >= 1)
  cl <- match.call()
  sepset <- vector("list",p)
  zMin <- c(0,rep(Inf,p))
  C <- mcor(cbind(y,dm), method = corMethod)
  cutoff <- qnorm(1 - alpha/2)
  n.edgetests <- numeric(1)# final length = max { ord}
  ## G := complete graph :
  G <- c(FALSE,rep(TRUE,p))
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
      if(verbose>=1) { if(i%%100==0) cat("|i=",i,"|iMax=",nrow(ind),"\n") }
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
            if (verbose==2) cat(paste("x:",colnames(dm)[x-1],"y:",(ytmp <- round((ncol(dm)+1)/2)),"S:"),c(ytmp,colnames(dm))[nbrs[S]],paste("z:",z,"\n"))
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
          }

        } } ## end if(!done)

    } ## end for(i ..)
    ord <- ord+1
  } ## end while
  Gres <- G[-1]
  names(Gres) <- colnames(dm)
  res <- list(G=Gres,zMin=zMin[-1])
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
          function(x, y, main = NULL, zvalue.lwd = FALSE, lwd.max = 7,
                   labels = NULL, ...)
      {
	if(is.null(main))
	    main <- deparse(x@call)
        attrs <- list()
        nodeAttrs <- list()
        if (!is.null(labels)) {
          attrs$node <- list(shape = "ellipse", fixedsize = FALSE)
          names(labels) <- nodes(x@graph)
          nodeAttrs$label <- labels
        }

        if (zvalue.lwd & numEdges(x@graph)!=0) {
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
          plot(z, main = main, ...)
        } else {
          plot(x@graph, nodeAttrs = nodeAttrs, main = main,
               attrs = attrs, ...)
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
 
##   if (length(S) < 4) {
    r <- pcorOrder(x,y, S, C)
##  } else {
##    k <- solve(C[c(x,y,S),c(x,y,S)])
##    r <- -k[1,2]/sqrt(k[1,1]*k[2,2])
    ##      r <- .C("parcorC",as.double(res),as.integer(x-1),as.integer(y-1),as.integer(S-1),as.integer(length(S)),as.integer(dim(C)[1]),as.double(as.vector(C)))[[1]]
##  }
  
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
  ml[ml != 0] <- rep(1,sum(ml != 0))

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
	   Qn_mu <- function(x, mu.too = FALSE, ...)
	     c(if(mu.too) median(x), Qn(x, ...))
	   cov2cor(covOGK(dm, n.iter = 2, sigmamu = Qn_mu,
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
  ## Purpose: Find sink of an adj matrix; return numeric(0) if there is none;
  ## a sink may have incident undirected edges, but no directed ones
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: Adjacency matrix (gm_i_j is edge from j to i)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:28

  ## treat undirected edges
  undir.nbrs <- which((gm==t(gm) & gm==1),arr.ind=TRUE)
  gm[undir.nbrs] <- 0
  ## treat directed edges
  res <- which(apply(gm,2,sum)==0)
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
  res.dag <- skel[i.ord,i.ord]
  res.dag
}

##################################################
## udag2pdag
##################################################
udag2pdag <- function(gInput,verbose=0) {
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
      
      if (length(allZ)>0) {
        for (j in 1:length(allZ)) {
          z <- allZ[j]
          if ((g[x,z]==0) & !((y %in% gInput@sepset[[x]][[z]]) |
                  (y %in% gInput@sepset[[z]][[x]]))) {
            if (verbose==1) {
              cat("\n",x,"->",y,"<-",z,"\n") 
              cat("Sxz=",gInput@sepset[[z]][[x]],"Szx=",gInput@sepset[[x]][[z]])
            }
            pdag[x,y] <- pdag[z,y] <- 1
            pdag[y,x] <- pdag[y,z] <- 0
          }
        }
      }
    }
    
    ## Test whether this pdag allows a consistent extension
    res2 <- pdag2dag(as(pdag,"graphNEL"))
    
    if (res2$success) {
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
              if (verbose==1) cat("\nRule 1:",a,"->",b," und ",b,"-",indC," wobei ",a," und ",indC," nicht verbunden: ",b,"->",indC,"\n")
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
              if (verbose==1) cat("\nRule 2: Kette ",a,"->",indC,"->",
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
                if (verbose==1) cat("\nRule 3:",a,"->",b,"\n")
              }
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
##-                   if (verbose==1) cat("Rule 4 applied \n")
##-                 }
##-               }
##-             }
##-           }
##-         }

      }
      res@graph <- as(pdag,"graphNEL")
    } else {
      res@graph <- res2$graph
    }
  }
  return(res)
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

decHeur <- function(dat,gam=0.05,sim.method="t",est.method="o",n.sim=100,two.sided=FALSE,verbose=FALSE)
{
  ## Purpose: Test wether data could come from a N or a N+t3 mix
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## dat: data (col=variables, rows=samples)
  ## gam: significance level for test
  ## sim.method: "n" Normal, "t" N + 10% t3
  ## est.method: "s" standard, "o" ogkQn
  ## n.sim: Number of samples for simulation
  ## two.sided: should two sided test be used?
  ## ----------------------------------------------------------------------
  ## Value:
  ## tvec: simulated values of test statistics under H0
  ## tval: value of test statistics for real data
  ## outlier: is robust method suggested? (TRUE=YES)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Feb 2008, 17:34

  p <- ncol(dat)
  n <- nrow(dat)
  mu <- rep(0,p)
  tstat <- function(dat) {mean(apply(dat,2,sd)/apply(dat,2,Qn))}
  is.outlier <- function(a,tvec,tval,two.sided=FALSE)
    {
      ## Purpose: Decide whether the test statistics lies in the rejection
      ## region of the distribution sampled by tvec using a test
      ## on significance level a
      ## ----------------------------------------------------------------------
      ## Arguments:
      ## a: significance level
      ## tvec: sampled values from null distribution of test statistics
      ## tval: actual value of test statistics
      ## two.sided: if FALSE, only for big tval rejected
      ## ----------------------------------------------------------------------
      ## Value: use.rob: bool, TRUE if tval is in rejection region
      ## ----------------------------------------------------------------------
      ## Author: Markus Kalisch, Date: 27 Feb 2008, 10:22

      n <- length(tvec)
      tvec.sorted <- sort(tvec)


      if (two.sided) {
        ind.low <- 1+floor((n-1)*a/2)
        ind.high <- 1+ceiling((n-1)*(1-a/2))
        c.low <- tvec.sorted[ind.low]
        c.high <- tvec.sorted[ind.high]
        c.int <- c(c.low,c.high)
        
        if ((tval<c.low)|(tval>c.high)) {
          use.rob <- TRUE
        } else {
          use.rob <- FALSE
        }
      } else {
        ind.high.1 <- 1+ceiling((n-1)*(1-a))
        c.high.1 <- tvec.sorted[ind.high.1]
        c.int <- c.high.1
        
        if (tval>c.high.1) {
          use.rob <- TRUE
        } else {
          use.rob <- FALSE
        }
      }
      ##  list(use.rob=use.rob,c.int=c.int)
      use.rob
    }

  t.start <- proc.time()[1]
  ## estimate correlation matrix of data
  if (est.method=="o") {
    mc <- mcor(dat,method="ogkQn")
  } else {
    mc <- mcor(dat,method="standard")
  }

  ## run simulations
  tvec <- rep(0,n.sim)
  for (i in 1:n.sim) {
    if (verbose) cat("Run ",i," out of ",n.sim,"\n")
    if (sim.method=="t") {
      d.sim <- rbind(mvrnorm(round(9*n/10),mu,mc),
                     rmt(n-round(9*n/10),mu,mc,3))
    } else {
      d.sim <- mvrnorm(n,mu,mc)
    }
    tvec[i] <- tstat(d.sim)
  }

  ## compute actual value of test statistics
  tval <- tstat(dat)

  ## is tval in rejection region of sig.level gam?
  use.rob <- is.outlier(gam,tvec,tval,two.sided=two.sided)
  rtime <- proc.time()[1]-t.start
  
  ## return results
  list(tvec=tvec,tval=tval,use.rob=use.rob)
    
}

################################################################################
## New in V8
## uses also library(vcd)
################################################################################
ci.test <- function(x,y,S=NULL,dm.df) {
  stopifnot(class(dm.df)=="data.frame",ncol(dm.df)>1)
  tab <- table(dm.df[,c(x,y,S)])
  if ((ncol(tab) < 2) | (nrow(tab)<2)) {
    res <- 1
  } else {
    if (length(S)==0) {
      res <- fisher.test(tab,simulate.p.value=TRUE)$p.value
    } else {
      res <- coindep_test(tab,3:(length(S)+2))$p.value
    }
  }
  res
}

pcAlgo <- function(dm = NA, C = NA, n=NA, alpha, corMethod =
                   "standard", verbose = FALSE, directed=FALSE,
                   G=NULL, datatype='continuous',NAdelete=TRUE,
                   m.max=Inf,u2pd="rand",psepset=FALSE) {
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
  
  if (any(is.na(dm))) {
    stopifnot(all(!is.na(C)),!is.na(n), (p <- ncol(C))>0)
  } else {
    n <- nrow(dm)
    p <- ncol(dm)
  }
  n <- as.integer(n)
  
  cl <- match.call()
  sepset <- vector("list",p)
  n.edgetests <- numeric(1)# final length = max { ord}
  if (is.null(G)) {
    ## G := complete graph :
    G <- matrix(rep(TRUE,p*p), nrow = p, ncol = p)
    diag(G) <- FALSE
  } else {
    if (!(identical(dim(G),c(p,p)))) {
      stop("Dimensions of the dataset and G do not agree.")
    }
  }
  seq_p <- 1:p
  for (iList in 1:p) sepset[[iList]] <- vector("list",p)
  zMin <- matrix(rep(Inf,p*p),nrow=p,ncol=p)

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
                if (is.na(prob)) prob <- ifelse(NAdelete,1,0)
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
    ##Orient colliders
    for (i in 1:dim(ind)[1]) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(amat[y,]==1),x) ## x-y-z
      
      if (length(allZ)>0) {
        for (j in 1:length(allZ)) {
          z <- allZ[j]
          if ((amat[x,z]==0) & !((y %in% sepset[[x]][[z]]) |
                     (y %in% sepset[[z]][[x]]))) {
            if (verbose >= 2) {
              cat("\n",x,"*->",y,"<-*",z,"\n") 
              cat("Sxz=",sepset[[z]][[x]],"and","Szx=",sepset[[x]][[z]],"\n")
            }

            ##x o-> y <-o z
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
          ##tf = possible_d_sep(amat,x,y)
          tf <- setdiff(tf1,y)
          ## test
          if (length(tf)>0) {
            z <- abs(zStat(x,y,tf,C,n))
            if (z < zMin[x,y]) zMin[x,y] <- z
            if (z <= cutoff) {
              ##delete x-y 
              amat[x, y] <- amat[y, x] <- 0
              ##save pos d-sepset in sepset
              sepset[[x]][[y]] <- tf
            }
            if (verbose>=2) {
              cat("Possible-D-Sep of", x, "and", y, "is", tf, " - z = ",z,"\n")
            }
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
  if (length(ind)>0) {
    for (i in 1:nrow(ind)) {
      x <- ind[i,]
      res[x[1],x[2]] <- amat[x[2],x[1]]
      res[x[2],x[1]] <- amat[x[1],x[2]]
    }
  }
  res
}

pdag2dag <- function(g,keepVstruct=TRUE) {
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
    succ <- TRUE
    res <- g
  } else {
    gm <- wgtMatrix(g) ## gm_i_j is edge from j to i
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
          if (!keepVstruct | adj.check(a,x)) {
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
    }
    if (go.on2==TRUE) {
      res <- as(amat2dag(gm),"graphNEL")
      warning("PDAG not extendible: Random DAG on skeleton drawn")
      succ <- FALSE
    } else {
      res <- as(t(gm2),"graphNEL")
      succ <- TRUE
    }
  }
  list(graph=res,success=succ)
}

udag2pdagSpecial <- function(gInput,verbose=0,n.max=100) {
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
      
      if (length(allZ)>0) {
        for (j in 1:length(allZ)) {
          z <- allZ[j]
          if ((g[x,z]==0) & !((y %in% gInput@sepset[[x]][[z]]) |
                  (y %in% gInput@sepset[[z]][[x]]))) {
            if (verbose==1) {
              cat("\n",x,"->",y,"<-",z,"\n") 
              cat("Sxz=",gInput@sepset[[z]][[x]],"Szx=",gInput@sepset[[x]][[z]])
            }
            ## check if already in other direction directed
            if (pdag[x,y]==0 & pdag[y,x]==1) { 
              evisit[x,y] <- evisit[x,y] + 1
              evisit[y,x] <- evisit[y,x] + 1
            }
            if (pdag[z,y]==0 & pdag[y,z]==1) {
              evisit[z,y] <- evisit[z,y] + 1
              evisit[y,z] <- evisit[y,z] + 1
            }
            pdag[x,y] <- pdag[z,y] <- 1
            pdag[y,x] <- pdag[y,z] <- 0
          } ## if
        } ## for 
      } ## if
    } ## for

    amat0 <- pdag
    ## Test whether this pdag allows a consistent extension
    res2 <- pdag2dag(as(pdag,"graphNEL"))
    xtbl <- res2$success
    xtbl.orig <- xtbl
    
    if (!xtbl & (max(evisit)>0)) {
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
        indBase[1,1:length(dgBase)] <- dgBase
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
              if (verbose==1) cat("\nRule 1:",a,"->",b," und ",b,"-",indC," wobei ",a," und ",indC," nicht verbunden: ",b,"->",indC,"\n")
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
              if (verbose==1) cat("\nRule 2: Kette ",a,"->",indC,"->",
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
                if (verbose==1) cat("\nRule 3:",a,"->",b,"\n")
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
  return(list(pcObj=res,evisit=evisit,xtbl=xtbl,xtbl.orig=xtbl.orig,amat0=amat0,amat1=amat1,status=status,counter=counter))
}

udag2pdagRelaxed <- function(gInput,verbose=0) {
  ## Purpose: Transform the Skeleton of a pcAlgo-object to a PDAG using
  ## the rules of Pearl. The output is again a pcAlgo-object. There is
  ## NO CHECK whether the resulting PDAG is really extendable.
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
      
      if (length(allZ)>0) {
        for (j in 1:length(allZ)) {
          z <- allZ[j]
          if ((g[x,z]==0) & !((y %in% gInput@sepset[[x]][[z]]) |
                  (y %in% gInput@sepset[[z]][[x]]))) {
            if (verbose==1) {
              cat("\n",x,"->",y,"<-",z,"\n") 
              cat("Sxz=",gInput@sepset[[z]][[x]],"Szx=",gInput@sepset[[x]][[z]])
            }
            pdag[x,y] <- pdag[z,y] <- 1
            pdag[y,x] <- pdag[y,z] <- 0
          }
        }
      }
    }
    
    ## Test whether this pdag allows a consistent extension
    ## res2 <- pdag2dag(as(pdag,"graphNEL"))
    
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
            if (verbose==1) cat("\nRule 1:",a,"->",b," und ",b,"-",indC," wobei ",a," und ",indC," nicht verbunden: ",b,"->",indC,"\n")
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
            if (verbose==1) cat("\nRule 2: Kette ",a,"->",indC,"->",
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
              if (verbose==1) cat("\nRule 3:",a,"->",b,"\n")
            }
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
      ##-                   if (verbose==1) cat("Rule 4 applied \n")
      ##-                 }
      ##-               }
      ##-             }
      ##-           }
      ##-         }

    }
    res@graph <- as(pdag,"graphNEL")

  }
  return(res)
}

beta.special <- function(dat=NA,x.pos,y.pos,verbose=0,a=0.01,myDAG=NA,myplot=FALSE,perfect=FALSE,method="local",collTest=TRUE,pcObj=NA,all.dags=NA,u2pd="rand")
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
  
  tmpColl <- FALSE

  ## Covariance matrix: Perfect case / standard case
  if (perfect) {
    if (class(myDAG)!="graphNEL") stop("For perfect-option the true DAG is needed!")
    mcov <- trueCov(myDAG)
    mcor <- cov2cor(mcov)
  } else {
    mcov <- cov(dat)
  }

  ## estimate skeleton and CPDAG of given data
  if (class(pcObj)!="pcAlgo") {
    if (perfect) {
      res <- pcAlgo.Perfect(mcor, corMethod = "standard",directed=TRUE,u2pd=u2pd)
    } else {
      res <- pcAlgo(dat, alpha = a, corMethod = "standard",directed=TRUE,u2pd=u2pd)
    }
  } else {
    res <- pcObj
  }

  ## prepare adjMatrix and skeleton
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
        for (i2 in 1:length(pa2)) {
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
        rev.pth <- sp.between(gDag,as.character(y.pos),
                              as.character(x.pos))[[1]]$path
        if (length(rev.pth)>1) {
          beta.hat[i] <- 0
        } else {
          pth <- sp.between(gDag,as.character(x.pos),
                            as.character(y.pos))[[1]]$path
          if (length(pth)<2) {
            beta.hat[i] <- 0
          } else {
            wgt.unique <- t(matrix(ad[i,],p,p)) ## wgt.est is wgtMatrix of DAG
            pa1 <- which(wgt.unique[x.pos,]!=0)
            if (y.pos %in% pa1) {
              cat("Y in Parents: ",y.pos," in ",pa1,"\n")
              beta.hat[i] <- 0
            } else {
              beta.hat[i] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
            }
            if (verbose==2) {cat("Fit - y:",y.pos,"x:",c(x.pos,pa1),
                  "|b.hat=",beta.hat,"\n")}
          } ## if length(pth)
        } ## if rev.pth
      } ## for n.dags
    } ## if n.dags
  } ## if method
  beta.hat
}


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

  if (is.na(amat) | is.na(amatSkel) | is.na(t.amat)) {
  ## Code for computing precomputable variables
    ## prepare adjMatrix and skeleton
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
      for (i2 in 1:length(pa2)) {
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

lm.cov <- function(C,y,x) {
  ## C: covariance matrix
  ## y: column of response
  ## x: columns of expl. vars
  sig <- C[x,x]
  beta <- solve(sig)%*%C[x,y,drop=FALSE]
  beta[1]
}

causalEffect <- function(g,y,x) {
  ## Compute true causal effect of x on y in g
  wmat <- wgtMatrix(g)
  p <- ncol(wmat)
  vec <- matrix(0,p,1)
  vec[x] <- 1
  if ((y-x)>1) {
    for (i in (x+1):y) vec[i] <- wmat[i,]%*%vec
    beta.true <- vec[y]
  } else {
    beta.true <- wmat[y,x]
  }
  beta.true
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
  if ((length(pa2.t)>0) & any(!is.na(pa2.t))) {
    ## check whether all pa1 and all pa2.t are connected;
    ## if no, there is a new collider
    if ((length(pa1)>0) & any(!is.na(pa1))) {
      res <- (min(amatSkel[pa1,pa2.t])==0) ## TRUE if new collider
    }
    ## in addition, all pa2.t have to be connected
    if ((length(pa2.t)>1) & (!res)) {
      tmp <- amatSkel[pa2.t,pa2.t]
      diag(tmp) <- 1
      res2 <- (min(tmp)==0) ## TRUE if new collider
      res <- (res|res2)
    }
  }
  if (!res & ((length(pa2.f)>0) & any(!is.na(pa2.f)))) {
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
    if (length(papa)==0) {
      res3 <- FALSE
    } else {
      res3 <- (min(amatSkel[x,papa])==0) ## TRUE if new collider
    }
    res <- (res|res3)
  }
  res
}
    
allDags <- function(gm,a,tmp,verb=FALSE)
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
    if (verb) {
      cat("Last Call - Final Graph: \n")
      print(gm)
      cat("####################\n")
    }
    tmp2 <- rbind(tmp,c(t(gm)))
    if (all(!duplicated(tmp2))) tmp <- tmp2
  } else {
    sinks <- find.sink(a)
    if (verb) {
      cat("Main Call: ################## \n")
        print(gm)
      print(a)
      cat("Sinks: ",sinks,"\n")
    }
    n.sinks <- length(sinks)
    if (n.sinks > 0) {
      for (i in 1:n.sinks) {
        if (verb) cat("Try removing",sinks[i]," in a.\n")
        gm2 <- gm
        a2 <- a
        x <- sinks[i]
        if (adj.check(a,x)) {
          inc.to.x <- which(a[, x] == 1 & a[x, ] == 1)
          if (length(inc.to.x) > 0) {
            real.inc.to.x <- as.numeric(row.names(a)[inc.to.x])
            real.x <- as.numeric(row.names(a)[x])
            gm2[real.x, real.inc.to.x] <- rep(1, length(inc.to.x))
            gm2[real.inc.to.x, real.x] <- rep(0, length(inc.to.x))
          }
          a2 <- a[-x,-x]
          if (verb) {
            cat("Removed sink",as.numeric(row.names(a)[x]),"in g (",
                sinks[i],"in a).\n")
            cat("New graphs: \n")
            print(gm2)
            print(a)
          }
          tmp <- allDags(gm2,a2,tmp,verb)
        }
      }
    }
  }
  tmp
}

pcAlgo.Perfect <- function(C, cutoff=0.00000001, corMethod = "standard", verbose = 0, directed=FALSE,u2pd="rand",psepset=FALSE) {
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
  if (verbose==FALSE) verbose <- 0
  if (verbose==TRUE) verbose <- 1
  p <- nrow(C)
  cl <- match.call()
  sepset <- vector("list",p)
  pcMin <- matrix(rep(Inf,p*p),nrow=p,ncol=p)
  diag(pcMin) <- rep(0,p)
  for (iList in 1:p) sepset[[iList]] <- vector("list",p)
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
          S <- seq(length = ord)

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
    ##Orient colliders
    for (i in 1:dim(ind)[1]) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(amat[y,]==1),x) ## x-y-z
      
      if (length(allZ)>0) {
        for (j in 1:length(allZ)) {
          z <- allZ[j]
          if ((amat[x,z]==0) & !((y %in% sepset[[x]][[z]]) |(y %in% sepset[[z]][[x]]))) {
            if (verbose == 2) {
              cat("\n",x,"*->",y,"<-*",z,"\n") 
              cat("Sxz=",sepset[[z]][[x]],"and","Szx=",sepset[[x]][[z]],"\n")
            }
            
            ##x o-> y <-o z
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
          ##tf = possible_d_sep(amat,x,y)
          tf <- setdiff(tf1,y)
          ## test
          if (length(tf)>0) {
            pc.val <- pcorOrder(x,y,tf,C)
            if (abs(pc.val)<pcMin[x,y]){
              pcMin[x,y] <- abs(pc.val)
            }
            if (abs(pc.val) <= cutoff) {
              ##delete x-y
              amat[x,y] <- amat[y,x] <- 0
              ##save pos d-sepset in sepset
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

  if (directed) {
    res <- switch (u2pd,
                   "rand" = udag2pdag(res),
                   "retry" = udag2pdagSpecial(res)$pcObj,
                   "relaxed" = udag2pdagRelaxed(res))
  }
  res
}

##Function that computes the Possible d-sepset, done by Spirtes
reach <- function(a,b,c,adjacency)
{
                                        #reachable      set of vertices;
                                        #edgeslist      array[1..maxvertex] of list of edges
                                        #more           Boolean
                                        #reachable      list of vertices
                                        #numvertex      integer
                                        #labeled        array (by depth) of list of edges that have been labeled

  makeedge = function(x,y)(list(list(x,y)))      

  legal <- function(t,...) {UseMethod("legal")}

  legal.possibledsep = function(t,u,v,r,s) {
    if (((adjacency[r[[1]],r[[2]]] == 2) && 
         (adjacency[s,r[[2]]] == 2) && 
         (r[[1]] != s)) ||
        ((adjacency[r[[1]],s] != 0) && 
         (r[[1]] != s))){
      edgeslist[[r[[2]]]]  <<- setdiff(edgeslist[[r[[2]]]],s)
      makeedge(r[[2]],s)}}

  legal.discriminating = function(t,u,v,r,s) {
    if ((length(intersect(s,c(t,u,v))) == 0) &&    # s not in the triangle t,u,v
        (adjacency[s,r[[2]]] == 2) &&            # s collides with r edge at r[[2]]
        (r[[1]] != s) &&                           # s is not on the r edge
        (adjacency[u,r[[2]]] == 3) &&
        (adjacency[r[[2]],u] == 2)){
      edgeslist[[r[[2]]]]  <<- setdiff(edgeslist[[r[[2]]]],s)
      makeedge(r[[2]],s)}}

  initialize <- function(x,...) {UseMethod("initialize")}

  initialize.possibledsep <- function(x,y,z) {mapply(makeedge,x=a,y=edgeslist[[a]])}

  initialize.discriminating <- function(x,y,z) {mapply(makeedge,x=a,y=setdiff(which(adjacency[,a] == 2),c(b,c)))}
  
  labeled = list()
  numvertex = dim(adjacency)[1]
  edgeslist = list()
  for (i in 1:numvertex) edgeslist = c(edgeslist,list(which(adjacency[,i] != 0)))
  labeled[[1]] = initialize(a,b,c)
  edgeslist[[a]] = list()
  depth = 2
  repeat
    {more = FALSE
     labeled[[depth]] = list()
     for (i in 1:length(labeled[[depth-1]])) {
       edgestemp = edgeslist[[labeled[[depth-1]][[i]][[2]]]]
       if (length(edgestemp) == 0) break
       for (j in 1:length(edgestemp))
         {labeled[[depth]] = 
            union(legal(a,b,c,labeled[[depth-1]][[i]],edgestemp[[j]]),labeled[[depth]])}}
     if (length(labeled[[depth]]) != 0){
       more = TRUE
       depth = depth  + 1} else break}
  return(unique(unlist(labeled)))
}

##Useful for R 4

discr.path <- function(path=NA, n=NA, pag=NA ,gInput=NA ,verb=NA)
{
  ## Purpose: find a discriminating path and if it exists orient
  ## the edge like in Rule 4. We start with a, b and c like
  ## in R4 and we go recursively to the left up to we find a d,
  ## such that there is a discriminating path between d and c for b
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - path: path[1]=c=gamma, path[2]=b=beta, path[3]=a=alpha: like in R 4
  ## - n: length of path
  ## - pag: adjacencies matrix
  ## - gInput: pc object (same input as FCI)
  ## - verb: 0 no comments, 1 detailed decription
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date:  6 Mar 2009, 11:10

  if (n<=dim(pag)[1]){
    res <- pag
    c <- path[1]
    b <- path[2]
    a <- path[3]
    first.pos <- path[n]
    del.pos <- path[n-1]
    ##all d with d *--> first.pos
    indD <- which (res[first.pos,]!=0 & res[,first.pos]==2)
    ##delete del.pos from indD
    indD <- setdiff(indD,del.pos)
    if (length(indD)>0){
      for (k in 1:length(indD)){
        d <- indD[k]
        ##check if we have already oriented the edge b o--* c
        ##res[c,b]==1 means that we have made no orientation, but if
        ##res[c,b]!=1 means that we have oriented and we don't look
        ##at other possible paths
        if ((res[c,b]==1) & !(d %in% path)){
          ##check if c and d are adjacent
          ##if no there exists a discr. path
          if (res[c,d]==0 & res[d,c]==0){
            ##check whether b is in Sepset(c,d)      
            if ((b %in% gInput@sepset[[c]][[d]]) |
                (b %in% gInput@sepset[[d]][[c]])){
              if(verb==1){
                cat("Rule 4: There is a discriminating path between:",d,"and",c,"for",b,"and",b, "is in Sepset of",c,"and",d,":",b,"->",c,"\n")
              }
              ##b --> c
              res[b,c] <- 2
              res[c,b] <- 3
            }
            else {
              if(verb==1){
                cat("Rule 4: There is a discriminating path between:",d,"and",c,"for",b,"and",b, "is not in Sepset of",c,"and",d,":",a,"<->",b,"<->",c,"\n")
              }
              ##a <--> b <--> c
              res[a,b] <- res[b,c] <- res[c,b] <- 2
              ##res[b,a]==2 from definition of a and b
            }
          }
          ##else if d and c are adjacent
          else {
            ##if d is a collider on the path
            ##and a parent of c
            if (res[first.pos,d]==2 & res[d,c]==2 & res[c,d]==3){
              new.path <- c(path,d)
              m <- length(new.path)
              res <- discr.path(path=new.path, n=m, pag=res, gInput=gInput, verb=verb)
            }
            ##else there exists no discr. path for this d
            else {
              res <- res
            }
          }## else
        }## if
      }## for k
    }## if length(indTheta)
  }## if
  return(res)
}



##Useful for R 5
ucp <- function(path=NA, pag=NA, n=NA, verb=NA)
{
  ## Purpose: find recursively if there exists an uncovered circle path
  ## and if TRUE orient the graph like in R5
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - path:a subgraph of pag for which we will find an ucp
  ## - pag: matrix representing the connections between tha whole graph
  ## - n: number of elements in path
  ## - verb: 0 no comments, 1 detailed decription
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date:  9 Feb 2009, 10:09

  
  if (n<=dim(pag)[1]){
    res <- pag
    ##all x with path[n-2] o-o x
    indX <- which((res[path[n-2],]==1 & res[,path[n-2]]==1))
    if (length(indX)>0){
      for (k in 1:length(indX)){
        x <- indX[k]
        ##res[path[1],path[n]]!=3 control to check if we have
        ##already directed this path
        if (!(x %in% path) & res[path[1],path[n]]!=3){
          new.path <- c(path[c(1:(n-2))],x,path[n-1],path[n])
          if (res[x,path[n-1]]==1 & res[path[n-1],x]==1){
            ##check uncovered feature
            check.uncov <- 0
            for (l in 1:(n-1)){
              if (res[new.path[l],new.path[l+2]]==0 & res[new.path[l+2],new.path[l]]==0){
                check.uncov <- check.uncov
              }
              else {
                check.uncov <- check.uncov + 1
              }
            }
            if (check.uncov==0){
              ##there exists an uncovered circle path
              res[new.path[1],new.path[n+1]] <- res[new.path[n+1],new.path[1]] <- 3  ##a -- b
              for (j in 1:n){
                res[new.path[j],new.path[j+1]] <- res[new.path[j+1],new.path[j]] <- 3 ##each edge on the path --
              }
              if(verb==1) {
                cat("Rule 5: There exists an uncovered circle path between",new.path[1],"and",new.path[n+1],":",new.path[1],"-", new.path[n+1],"and for each edge on the path",new.path,"\n")
              }
            }
            else {
              ##recursion
              ##check that alpha and new.path[3] are not connected
              if (res[new.path[1],new.path[3]]==0 & res[new.path[3],new.path[1]]==0){
                res <- ucp(path=new.path, pag=res, n=length(new.path), verb=verb)
              }
            }
          }
          else {
            res <- ucp(path=new.path, pag=res, n=length(new.path), verb=verb)
          }
        }##if delta %in% path
      }##for k
    }## if length(indDelta)
  }##if
  return(res)
}



##Useful for R 9
upd <- function(path=NA, pag=NA, n=NA, verb=NA)
{
  ## Purpose:find recursively an uncovered potentially directed path
  ## in pag and orient the edge like in R9
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - path: a subgraph of pag for which we will find an upd, like in R 9
  ## - pag: matrix representing the connections between tha whole graph
  ## - n: number of elements in path
  ## - verb: 0 no comments, 1 detailed decription
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date:  9 Feb 2009, 11:33

  if (n<=dim(pag)[1]){
    res <- pag
    c <- path[1]
    a <- path[2]
    b <- path[n-1]
    d <- path[n]
    ##find all x such that d @--@ x (@== potentially directed)
    ##and not adjacent to b
    indX <- which((res[d,]==2 | res[d,]==1) & (res[,d]==1 | res[,d]==3) & (res[b,]==0 & res[,b]==0))
    if (length(indX)>0){
      for (j in 1:length(indX)){
        x <- indX[j]
        ##res[c,a]!=3 control to see if we have already directed
        ##the path a *--> c
        ##and check that x @--@ c
        if ((!(x %in% path)) & (res[c,a]!=3)){
          new.path <- c(path[c(2:n)],x,path[1])
          if ((res[x,c]==1 | res[x,c]==2) & (res[c,x]==1 | res[c,x]==3)){    
            ##check uncovered feature
            check.uncov <- 0
            for (l in 1:(n-1)){
              if (res[new.path[l],new.path[l+2]]==0 & res[new.path[l+2],new.path[l]]==0){
                check.uncov <- check.uncov
              }
              else {
                check.uncov <- check.uncov + 1
              }
            }
            if (check.uncov==0){
              res[c,a] <- 3
              if (verb==1) {
                cat("Rule 9: There exists an upd between",new.path,":",a,"->",c,"\n")
              }
            }
            else {
              rec.path <- c(path,x)
              res <- upd(path=rec.path, pag=res, n=length(rec.path), verb=verb)
            }
          }## if
          else {
            rec.path <- c(path,x)
            res <- upd(path=rec.path, pag=res, n=length(rec.path), verb=verb)
          }
        }## if
      }## for j
    }## if length(indDelta)
  }## if
  return(res)
}




##Useful for R 10
find.upd <- function(path=NA, a=NA, n=NA, pag=NA, verb=NA)
{
  ## Purpose: find if there exists an uncovered potentially
  ## directed path between the first and the last element in path
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - path: like in R 10, that is first element adjacent to alpha and
  ##   last element beta or theta
  ## - a: alpha
  ## - n: length of pag
  ## - pag: adjacencies matrix
  ## - verb: 0 no comments, 1 detailed decription
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date:  9 Feb 2009, 14:00

  if (n<=dim(pag)[1]){
    final <- path[n]
    mittle <- path[n-1]
    uncov.path <- NA
    res <- FALSE
    if (n==2){
      if(path[1]==path[n]){
        uncov.path <- path
        res <- TRUE
      }
      else {
        if ((pag[path[1],path[2]]==1 | pag[path[1],path[2]]==2) & (pag[path[2],path[1]]==1 | pag[path[2],path[1]]==3)){
          uncov.path <- path
          res <- TRUE
        }
      }
    }
    else {
      ##find all b s.t. mittle @-@ b (@ potentially directed)
      indB <- which((pag[mittle,]==1 | pag[mittle,]==2) & (pag[,mittle]==1 | pag[,mittle]==3))
      indB <- setdiff(indB,a)
      if (length(indB)>0){
        counter <- 0
        while ((!res) & (counter < length(indB))){
          counter <- counter + 1
          b <- indB[counter]
          if (!(b %in% path)){
            new.path <- c(path[1:(n-1)],b,path[n])
            check.uncov <- 0
            for (l in 1:(n-1)){
              if (pag[new.path[l],new.path[l+2]]==0 & pag[new.path[l+2],new.path[l]]==0){
                check.uncov <- check.uncov
              }
              else {
                check.uncov <- check.uncov + 1
              }
            }
            if (check.uncov==0){
              ##check if b @-@ final
              if ((pag[b,final]==1 | pag[b,final]==2) & (pag[final,b]==1 | pag[final,b]==3)){
                uncov.path <- new.path
                res <- TRUE
              }
              else {
                tmp <- find.upd(path=new.path, n=length(new.path), pag=pag, verb=verb)
                uncov.path <- tmp[[2]]
                res <- tmp[[1]]
              }
            }
            else {
              tmp <- find.upd(path=new.path, n=length(new.path), pag=pag, verb=verb)
              uncov.path <- tmp[[2]]
              res <- tmp[[1]]
            }
          }##if
        }##while
      }##if
    }##else
  }##if
  return(list(res=res,uncov.path=uncov.path))
}




udag2pag <- function(gInput, rules=rep(TRUE,10), verbose=TRUE) {
  ## Purpose:Transform the Skeleton of a pcAlgo-object to a PAG using
  ## the rules of Zhang. The output is an adjacency matrix.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gInput: pcAlgo object
  ## - rules: array of length 10 wich contains TRUE or FALSE corrsponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - verbose: 0 - no output, 1 - detailed output
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date:  6 Mar 2009, 11:27


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
  
  ##counter <- 0
  ##res <- gInput
  ##status <- 0
  ##xtbl <- xtbl.orig <- TRUE
  if (numEdges(gInput@graph)>0) {
    g <- as(gInput@graph,"matrix") ## g_ij==1 if i *-o j
    p <- dim(g)[1]
    evisit <- amat0 <- amat1 <- matrix(0,p,p)
    pag <- g
    ind <- which(g==1, arr.ind=TRUE)
    ## ind <- unique(t(apply(ind,1,sort)))

    
    ## Create minimal pattern
    ## rule 0
    for (i in 1:dim(ind)[1]) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(g[y,]==1),x) ## x-y-z
      
      if (length(allZ)>0) {
        for (j in 1:length(allZ)) {
          z <- allZ[j]
          if ((g[x,z]==0) & !((y %in% gInput@sepset[[x]][[z]]) |
                  (y %in% gInput@sepset[[z]][[x]]))) {
            if (verbose==1) {
              cat("\n",x,"*->",y,"<-*",z,"\n") 
              cat("Sxz=",gInput@sepset[[z]][[x]],"and","Szx=",gInput@sepset[[x]][[z]],"\n")
            }

            ##x o-> y <-o z
            pag[x,y] <- pag[z,y] <- 2

          } ## if
        } ## for 
      } ## if
    } ## for
    

    ## first while for R1-R4
    old_pag1 <- matrix(rep(0,p^2),nrow=p,ncol=p)
    while (sum(!(old_pag1==pag))>0) {
      old_pag1 <- pag

    
      ## rule 1
      if (rules[1]){
        ind <- which((pag==2 & t(pag)!=0), arr.ind=TRUE) ## a *--> b
        if (length(ind)>0) {
          for (i in 1:dim(ind)[1]) {
            a <- ind[i,1]
            b <- ind[i,2]
            ##all c with b o--* c und a not adjacent to c
            indC <- which((pag[b,]!=0 & pag[,b]==1) & (pag[a,]==0 & pag[,a]==0))
            indC <- setdiff(indC,a)
            if (length(indC)>0) {
              pag[b,indC] <- 2
              pag[indC,b] <- 3
              if (verbose==1) {
                cat("Rule 1:",a,"*->",b,"o-*",indC,
                    " where ",a," and ",indC," not connected: ",b,"->",indC,"\n")
              }
            }## if
          }## for
        }## if
      }
      
      
      
      ## rule 2
      if (rules[2]){
        ind <- which((pag==1 & t(pag)!=0), arr.ind=TRUE) ## a *--o c
        if (length(ind)>0) {
          for (i in 1:dim(ind)[1]) {
            a <- ind[i,1]
            c <- ind[i,2]
            ##all b with a --> b *--> c or a *--> b --> c
            indB <- which( ((pag[a,]==2 & pag[,a]==3) & (pag[c,]!=0 & pag[,c]==2)) | ((pag[a,]==2 & pag[,a]!=0) & (pag[c,]==3 & pag[,c]==2)))
            if (length(indB)>0) {
              pag[a,c] <- 2
              if (verbose==1) {
                cat("Rule 2:",a,"->",indB,"*->",
                    c,"or",a,"*->",indB,"->",c,"and",a,"*-o",c,":",a,"*->",c,"\n")
              }
            }## if
          }## for
        }## if
      }
      
      
      
      ## rule 3
      if (rules[3]){
        ind <- which((pag!=0 & t(pag)==1), arr.ind=TRUE) ## b o--* d
        if (length(ind)>0) {
          for (i in 1:dim(ind)[1]) {
            b <- ind[i,1]
            d <- ind[i,2]
            ##all a and c with a *--> b and a *--o d, c *--> b and c *--o d
            indAC <- which( (pag[b,]!=0 & pag[,b]==2) & (pag[,d]==1 & pag[d,]!=0))
            ##check if there are a and c which are not adjacent
            if (length(indAC)>=2){
              counter <- 0
              while ((counter<(length(indAC)-1)) & (pag[d,b]!=2)){
                counter <- counter + 1
                ii <- counter
                while ((ii<length(indAC)) & (pag[d,b]!=2)){
                  ii <- ii + 1
                  if (pag[indAC[counter],indAC[ii]]==0 & pag[indAC[ii],indAC[counter]]==0){
                    pag[d,b] <- 2
                    if (verbose==1){
                      cat("Rule 3:",d,"*->",b,"\n")
                    }
                  }
                }
              }
            }
            ## if (length(indAC)>1){
            ## for(j in 1:(length(indAC)-1)){
            ## for(k in (j+1):(length(indAC))){
            ## if(pag[indAC[j],indAC[k]]==0 & pag[indAC[k],indAC[j]]==0){
            ## pag[d,b] <- 2
            ## if (verbose==1) {
            ## cat("Rule 3:",d,"*->",b,"\n")
            ## }
            ## }## if
            ##}## for
            ##}## for
            ##}## if
          }## for
        }## if
      }
      
      
      ## rule 4
      if (rules[4]){
        ind <- which((pag!=0 & t(pag)==1), arr.ind=TRUE) ## b o--* c
        if (length(ind)>0){
          for (i in 1:dim(ind)[1]){
            b <- ind[i,1]
            c <- ind[i,2]
            ##all a with a <--* b o--* c and a --> c
            indA <- which( (pag[b,]==2 & pag[,b]!=0) & (pag[c,]==3 & pag[,c]==2))
            if (length(indA)>0){
              for (j in 1:length(indA)){
                a <- indA[j]
                path <- c(c,b,a)
                n <- length(path)
                pag <- discr.path(path=path, n=n , pag=pag, gInput=gInput, verb=verbose)
              }## for j
            }## if length(indA)
          }## for i
        }## if length(ind)
      }
    }## while



    ## second while for R5-R10
    old_pag2 <- matrix(rep(0,p^2),nrow=p,ncol=p)
    while (sum(!(old_pag2==pag))>0) {
      old_pag2 <- pag
    
  
      ##rule 5
      if (rules[5]){
        ind <- which((pag==1 & t(pag)==1), arr.ind=TRUE) ## a o--o b
        if (length(ind)>0){
          for (i in 1:dim(ind)[1]){
            a <- ind[i,1]
            b <- ind[i,2]
            ##all c s.t. a o--o c and c not adjacent to b
            indC <- which((pag[a,]==1 & pag[,a]==1) & (pag[b,]==0 & pag[,b]==0))
            ##delete b from indC
            indC <- setdiff(indC,b)
            ##all d s.t. d o--o b and d not adjacent to a
            indD <- which((pag[b,]==1 & pag[,b]==1) & (pag[a,]==0 & pag[,a]==0))
            ##delete a from indD
            indD <- setdiff(indD,a)
            if (length(indC)>0 & length(indD)>0){
              for (j in 1:length(indC)){
                c <- indC[j]
                for (l in 1:length(indD)){
                  d <- indD[l]
                  if (pag[c,d]==1 & pag[d,c]==1){
                    pag[a,b] <- pag[b,a] <- 3
                    pag[a,c] <- pag[c,a] <- 3
                    pag[c,d] <- pag[d,c] <- 3
                    pag[d,b] <- pag[b,d] <- 3
                    if (verbose==1) {
                      cat("Rule 5: There exists an uncovered circle path between",a,"and",b,":",a,"-", b,"and",a,"-",c,"-",d,"-",b, "\n")
                    }
                  }
                  else {
                    path <- c(a,c,d,b)
                    pag <- ucp(path=path, pag=pag, n=length(path), verb=verbose)
                  }
                }##for l
              }##for j
            }##if
          }##for i
        }##if
      }
      
      
      
      ## rule 6
      if (rules[6]){
        ind <- which((pag!=0 & t(pag)==1), arr.ind=TRUE) ## b o--* c
        if (length(ind)>0){
          for (i in 1:dim(ind)[1]){
            b <- ind[i,1]
            c <- ind[i,2]
            ##all a with a -- b o--* c
            indA <- which(pag[b,]==3 & pag[,b]==3)
            if (length(indA)>0){
              pag[c,b] <- 3
              if (verbose==1) {
                cat("Rule 6:",a,"-",b,"o-*",c,":",b,"-*",c,"\n")
              }
            }
          }## for i
        }## if length(ind)
      }
      
      
      
      ## rule 7
      if (rules[7]){
        ind <- which((pag!=0 & t(pag)==1), arr.ind=TRUE) ## b o--* c
        if (length(ind)>0){
          for (i in 1:dim(ind)[1]){
            b <- ind[i,1]
            c <- ind[i,2]
            ##all a with a --o b and a,c not adjacent
            indA <- which((pag[b,]==3 & pag[,b]==1) & (pag[c,]==0 & pag[,c]==0))
            indA <- setdiff(indA,c)
            if (length(indA)>0){
              pag[c,b] <- 3
              if (verbose==1) {
                cat("Rule 7:",indA,"-o",b,"o-*",c,"and", indA," and", c," are not adjacent:",b,"-*",c,"\n")
              }
            }
          }## for i
        }## if length(ind)
      }
      
      
      ## rule 8
      if (rules[8]){
        ind <- which((pag==2 & t(pag)==1), arr.ind=TRUE) ## a o--> c
        if (length(ind)>0){
          for (i in 1:dim(ind)[1]){
            a <- ind[i,1]
            c <- ind[i,2]
            ##all b with a --> b --> c or a --o b --> c
            indB <- which(((pag[a,]==2 & pag[,a]==3) | (pag[a,]==1 & pag[,a]==3)) & (pag[c,]==3 & pag[,c]==2))
            if (length(indB)>0){
              pag[c,a] <- 3
              if (verbose==1) {
                cat("Rule 8:",a,"->",indB,"->",c,"or",a,"-o",indB,"->",c,"and",a,"o->",c,":",a,"->",c,"\n")
              }
            }
          }## for i
        }## if length(ind)
      }
      
      
      
      ##rule 9
      if (rules[9]){
        ind <- which((pag==2 & t(pag)==1), arr.ind=TRUE) ## a o--> c
        if (length(ind)>0){
          for (i in 1:dim(ind)[1]){
            a <- ind[i,1]
            c <- ind[i,2]
            ##find all beta s.t. a @--@ b is potentially directed
            ##and beta is not connected to c
            indB <- which((pag[a,]==2 | pag[a,]==1) & (pag[,a]==1 | pag[,a]==3) & (pag[c,]==0 & pag[,c]==0))
            indB <- setdiff(indB,c)
            if (length(indB)>0){
              for (k in 1:length(indB)){   
                b <- indB[k]
                ##find all d such that c <--o a @-@ b @-@ d and d and a not adjacent
                indD <- which((pag[b,]==2 | pag[b,]==1) & (pag[,b]==1 | pag[,b]==3) & (pag[a,]==0 & pag[,a]==0))
                indD <- setdiff(indD,a)
                if (length(indD)>0){
                  for (l in 1:length(indD)){
                    d <- indD[l]
                    ##pag[c,a]!=3 check if we have already oriented this edge
                    if (pag[c,a]!=3){
                      if ((pag[c,d]==1 | pag[c,d]==3) & (pag[d,c]==1 | pag[d,c]==2)){
                        pag[c,a] <- 3
                        if (verbose==1){
                          cat("Rule 9, There exists an upd between",a,"and",c,":",a," ->", c, "\n")
                        }
                      }
                      else {
                        path <- c(c,a,b,d)
                        pag <- upd(path=path, pag=pag, n=length(path), verb=verbose)
                      }
                    }
                  }##for l
                }## if lenth(indD)
              }##for k
            }##if lenfth(indB)
          }##for i
        }## if length(ind)
      }


      ##rule 10
      if (rules[10]){
        ind <- which((pag==2 & t(pag)==1), arr.ind=TRUE) ## a o--> c
        if (length(ind)>0){
          for (i in 1:dim(ind)[1]){
            a <- ind[i,1]
            c <- ind[i,2]
            ##find all b s.t. b --> c
            indB <- which((pag[c,]==3 & pag[,c]==2))
            if (length(indB)>=2){
              for (j in 1:length(indB)){
                b <- indB[j]
                ##chose a d s.t. d --> c
                indD <- setdiff(indB,b)
                if (length(indD)>0 & pag[c,a]!=3){
                  for (k in 1:length(indD)){
                    d <- indD[k]
                    ##case where mu=b and omega=theta and mu and omega not adj
                    if ((pag[a,b]==1 | pag[a,b]==2) & (pag[b,a]==1 | pag[b,a]==3) & (pag[a,d]==1 | pag[a,d]==2) & (pag[d,a]==1 | pag[d,a]==3) & (pag[d,b]==0 & pag[b,d]==0)){
                      pag[c,a] <- 3
                      if (verbose==1){
                        cat("Rule 10 with mu = beta = ",b,"and omega = theta =",d,":",a,"->",c,"\n")
                      }
                    }
                    else {
                      indA <- which((pag[a,]==1 | pag[a,]==2) & (pag[,a]==1 | pag[,a]==3), arr.ind=TRUE)
                      indA <- setdiff(indA,c)
                      if (length(indA>=2)){
                        for (l in 1:length(indA)){
                          first.pos <- indA[l]
                          ##if (first.pos==b){
                          ##indAA <- setdiff(indA,c(first.pos,d))
                          ##}
                          ##if (first.pos==d){
                          ##indAA <- setdiff(indA,c(first.pos,b))
                          ##}
                          ##else {
                          indAA <- setdiff(indA,first.pos)
                          if ((length(indAA)>0) & (pag[c,a]!=3)){
                            for (s in 1:length(indAA)){
                              sec.pos <- indAA[s]
                              p1 <- find.upd(path=c(first.pos,b), a=a, n=2, pag=pag, verb=verbose)
                              p2 <- find.upd(path=c(sec.pos,d), a=a, n=2, pag=pag, verb=verbose)
                              if (p1$res==TRUE & p2$res==TRUE){
                                mu <- p1$uncov.path[1]
                                omega <- p2$uncov.path[1]
                                if ((mu!=omega) & (pag[mu,omega]==0) & (pag[omega,mu]==0)){
                                  ##then orient
                                  pag[c,a] <- 3
                                  if (verbose==1) {
                                    cat("Rule 10:",a,"->",c,"\n")
                                  }
                                }
                              }
                            }## for s
                          }## if
                          ##}## else
                        }## for l
                      }## if
                    }##else
                  }##for k
                }##if
              }##for j
            }##if
          }##for i
        }##if
      }
    
    }## while
    
  }## if

  return(pag)

}## function



plotAG <- function(amat)
{
  ## Purpose: Plot ancestral graph
  ## ATTENTION: Start with M-x R-devel; uses library(Rgraphviz)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: Adjacency matrix
  ##   amat[i,j]=3 & amat[j,i]=1 iff i 1-3 j
  ##   "0": no edge; "1": circle; "2": arrow; "3": tail
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 16 Feb 2009, 18:01

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
  
  edgeRenderInfo(g) <- list(arrowhead=ah.list,arrowtail=at.list)
  renderGraph(layoutGraph(g))
}
