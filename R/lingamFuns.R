##################################################
## exported
##################################################
LINGAM <-
function(X, verbose = FALSE)
# Copyright (c) 2013 - 2015  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    res <- uselingam(t(X), verbose = verbose)
    return(list(B = res$B, Adj = t(res$B != 0)))
}

##################################################
## internal
##################################################

all.perm <- function(n) {
  p <- matrix(1, ncol = 1)
  for (i in 2:n) {
    p <- pp <- cbind(p, i)
    v <- c(1:i, 1:(i - 1))
    for (j in 2:i) {
      v <- v[-1]
      p <- rbind(p, pp[, v[1:i]])
    }
  }
  p
}

estimate <- function( X , verbose = FALSE) {
    
    # Using the fastICA R package so make sure it is loaded
    # library('fastICA')
    
    # Call the fastICA algorithm
    icares <- fastICA( t(X), nrow(X), tol=1e-14 )
    W <- t((icares$K) %*% (icares$W))
    
    # [Here, we really should perform some tests to see if the 'icasig' 
    # really are independent. If they are not very independent, we should
    # issue a warning. This is not yet implemented.]
    
    # Try to permute the rows of W so that sum(1./abs(diag(W))) is minimized
    if(verbose)
    {
        cat('Performing row permutation...\n');
    }
    dims <- nrow(X)
    if (dims <= 8) {  
        if(verbose)
        {
            cat('(Small dimensionality, using brute-force method.)\n')
        }
        temp <- nzdiagbruteforce( W )
        Wp <- temp$Wopt
        rowp <- temp$rowp
    } else
    {
        if(verbose)
        {
            cat('(Using the Hungarian algorithm.)\n')
        }
        temp <- nzdiaghungarian( W )
        Wp <- temp$Wopt
        rowp <- temp$rowp
    }
    if(verbose)
    {
        cat('Done!\n')
    }
    # Divide each row of Wp by the diagonal element
    estdisturbancestd <- 1/diag(abs(Wp))
    Wp <- Wp/diag(Wp)
    
    # Compute corresponding B
    Best <- diag(dims)-Wp
    
    # Estimate the constants c
    m <- rowMeans(X)
    dim(m) <- c(dims,1)
    cest <- Wp %*% m
    
    # Next, identically permute the rows and columns of B so as to get an
    # approximately strictly lower triangular matrix
    if(verbose)
    {
        cat('Performing permutation for causal order...\n');
    }
    if (dims <= 8) {  
        if(verbose)
        {
            cat('(Small dimensionality, using brute-force method.)\n');
        }
        temp <- sltbruteforce( Best )
        Bestcausal <- temp$Bopt
        causalperm <- temp$optperm
    } else
    {
        if(verbose)
        {
            cat('(Using pruning algorithm.)\n')
        }
        temp <- sltprune( Best )
        Bestcausal <- temp$Bopt
        causalperm <- temp$optperm
        
    }
    if(verbose)
    {
        cat('Done!\n');
    }
    
    if(verbose)
    {
        print(Bestcausal)
    }
    # Here, we report how lower triangular the result was, and in 
    # particular we issue a warning if it was not so good!
    percentinupper <- sltscore(Bestcausal)/sum(Bestcausal^2)
    if(verbose)
    {   
        if (percentinupper>0.2) 
            cat('WARNING: Causal B not really triangular at all!!\n')
        else if (percentinupper>0.05)
            cat('WARNING: Causal B only somewhat triangular!\n')
        else
            cat('Causal B nicely triangular. No problems to report here.\n')
    }
    # Set the upper triangular to zero
    Bestcausal[upper.tri(Bestcausal,diag=FALSE)] <- 0
    
    # Finally, permute 'Bestcausal' back to the original variable
    # ordering and rename all the variables to the way we defined them
    # in the function definition
    icausal <- iperm(causalperm);
    res <- list()
    res$B <- Bestcausal[icausal, icausal];
    res$stde <- estdisturbancestd
    res$ci <- cest
    res$k <- causalperm
    res$W <- W
    res$pinup <- percentinupper # added for firmgrowth conintegration analysis
    
    # Return the result
    return(res)
    
}

iperm <- function( p ) {

  q <- array(0,c(1,length(p)))
  
  for (i in 1:length(p)) {
    ind <- which(p==i)
    q[i] <- ind[1]
  }
  
  q
  
}

nzdiagbruteforce <- function( W ) {

  #--------------------------------------------------------------------------
  # Try all row permutations, find best solution
  #--------------------------------------------------------------------------

  n <- nrow(W)
  
  bestval <- Inf;
  besti <- 0;
  allperms <- all.perm(n) 
  nperms <- nrow(allperms)
  
  for (i in 1:nperms) {
    Pr <- diag(n)
    Pr <- Pr[,allperms[i,]]
    Wtilde <- Pr %*% W
    c <- nzdiagscore(Wtilde)
    if (c<bestval) {
      bestWtilde <- Wtilde
      bestval <- c
      besti <- i
    }
  }

  res <- list()
  res$Wopt <- bestWtilde
  res$rowp <- allperms[besti,]
  res$rowp <- iperm(res$rowp)
  
  res
  
}

nzdiaghungarian <-
function( W )
{
    # Jonas' quick additional hack to make lingam work in problems with large p 
    # (the number of dimensions from which this code is used is specified in "estimate.R")
    n <- nrow(W)
    S <- matrix(1,n,n)/abs(W)
    
    ####
    ###[c,T]=hungarian(S');
    ###
    c <- as.numeric(solve_LSAP(S))
    
    # Permute W to get Wopt
    Pr <- diag(n)
    Pr <- Pr[,c]
    Wopt <- Pr %*% W
    
    
    # Return the optimal permutation as well
    res <- list()
    res$Wopt <- Wopt
    res$rowp <- iperm(c)
    return(res)
    
}

nzdiagscore <- function( W ) {

  res <- sum(1/diag(abs(W)))
  res
  
}

prune <- function( X, k, verbose=FALSE ) {
    
    # ---------------------------------------------------------------------------
    # Default values for parameters
    # ---------------------------------------------------------------------------
    
    method <- 'resampling'  # the pruning method
    
    # 'prunefactor' determines how easily weak connections are pruned
    # in the simple resampling based pruning. Zero gives no pruning 
    # whatsoever. We use a default of 1, but remember that it is quite arbitrary
    # (similar to setting a significance threshold). It should probably
    # become an input parameter along with the data, but this is left
    # to some future version.
    prunefactor <- 1
    
    # ---------------------------------------------------------------------------
    # Pruning
    # ---------------------------------------------------------------------------
    
    if(verbose)
    {
        cat('Pruning the network connections...\n')
    }
    dims <- nrow(X)
    ndata <- ncol(X)
    
    if (method == 'resampling') {
        
        # -------------------------------------------------------------------------
        # Pruning based on resampling: divide the data into several equally
        # sized pieces, calculate B using covariance and QR for each piece
        # (using the estimated causal order), and then use these multiple
        # estimates of B to determine the mean and variance of each element.
        # Prune the network using these.
        
        npieces <- 10
        piecesize <- floor(ndata/npieces)
        
        Bpieces <- array(0,c(dims,dims,npieces))
        diststdpieces <- array(0,c(dims,npieces))
        cpieces <- array(0,c(dims,npieces))
        Bfinal <- matrix(0,dims,dims)
        
        for (i in 1:npieces) {
            
            # Select subset of data, and permute the variables to the causal order
            Xp <- X[k,((i-1)*piecesize+1):(i*piecesize)]
            
            # Remember to subract out the mean 
            Xpm <- rowMeans(Xp)
            Xp <- Xp - Xpm
            
            # Calculate covariance matrix
            C <- (Xp %*% t(Xp))/ncol(Xp)
            
            
            ### HACK BY JONAS PETERS 05/2013
            #if(min(eigen(C)$values) < 0)
            if(1==1)
            {
                # regularization
                C <- C + diag(10^(-10),dim(C)[1])
            }
            while(min(eigen(C)$values) < 0)
            {
                # regularization
                C <- C + diag(10^(-10),dim(C)[1])
            }
            ### HACK BY JONAS PETERS 05/2013
            
            
            # Do QL decomposition on the inverse square root of C
            res <- tridecomp(solve(sqrtm(C)),'ql');
            Q <- res$A
            L <- res$B
            
            # The estimated disturbance-stds are one over the abs of the diag of L
            newestdisturbancestd <- 1/diag(abs(L));
            
            # Normalize rows of L to unit diagonal
            L <- L/diag(L)
            
            # Calculate corresponding B
            Bnewest <- diag(dims)-L;
            
            # Also calculate constants
            cnewest <- L %*% Xpm
            
            # Permute back to original variable order
            ik <- as.vector(iperm(k));
            Bnewest <- Bnewest[ik, ik]
            newestdisturbancestd <- newestdisturbancestd[ik]
            cnewest <- cnewest[ik]
            
            # Save results
            Bpieces[,,i] <- Bnewest
            diststdpieces[,i] <- newestdisturbancestd
            cpieces[,i] <- cnewest
            
        }
        
        for (i in 1:dims) {
            for (j in 1:dims) {
                
                themean <- mean(Bpieces[i,j,])
                thestd <- sd(Bpieces[i,j,])
                if (abs(themean)<prunefactor*thestd) {	    
                    Bfinal[i,j] <- 0
                }
                else {
                    Bfinal[i,j] <- themean
                }
                
            }
        }
        
        diststdfinal <- rowMeans(diststdpieces)
        cfinal <- rowMeans(cpieces)
        
        # Finally, rename all the variables to the way we defined them
        # in the function definition
        
        Bpruned <- Bfinal
        stde <- diststdfinal
        ci <- cfinal
        
    }
    
    if (method == 'olsboot') {
        stop('Not implemented yet!')
    }
    
    if (method == 'wald') {
        stop('Not implemented yet!')
    }
    
    if (method == 'bonferroni') {
        stop('Not implemented yet!')
    }
    
    if (method == 'hochberg') {
        stop('Not implemented yet!')
    }
    
    if (method == 'modelfit') {
        stop('Not implemented yet!')
    }
    
    if(verbose)
    {
        cat('Done!\n')
    }
    
    # Return the result
    res <- list()
    res$Bpruned <- Bpruned
    res$stde <- stde
    res$ci <- ci
    return(res)
    
}

sltbruteforce <- function( B ) {

  #--------------------------------------------------------------------------
  # Try all row permutations, find best solution
  #--------------------------------------------------------------------------

  n <- nrow(B)
  
  bestval <- Inf;
  besti <- 0;
  allperms <- all.perm(n) 
  nperms <- nrow(allperms)
  
  for (i in 1:nperms) {
    Btilde <- B[allperms[i,],allperms[i,]]
    c <- sltscore(Btilde)
    if (c<bestval) {
      bestBtilde <- Btilde
      bestval <- c
      besti <- i
    }
  }

  res <- list()
  res$Bopt <- bestBtilde
  res$optperm <- allperms[besti,]
  
  res

}

sltprune <- function( B ) {

    #Hack of JONAS PETERS 2013
    n <- nrow(B)
    rr <- sort(abs(B), index.return = TRUE)
    
    #[y,ind] = sort(abs(B(:)));
    ind <- rr$ix
    
    for(i in ((n*(n+1)/2):(n*n)))
    {
        
        # Copy original B into Bi
        Bi <- B
        
        # Set 'i' smallest (in absolute value) coefficients to zero
        Bi[ind[1:i]] <- 0
        
        # Try to do permutation
        p <- slttestperm( Bi )
        
        # If we succeeded, then we're done!
        if(any(p != 0))
        {
            Bopt <- B[p,p]
            optperm <- p
            return(list(Bopt = Bopt, optperm = p))
        }
        # ...else we continue, setting one more to zero!
    }
    return(list(Bopt = Bopt, optperm = p))

}

sltscore <- function( B ) {

  s <- sum((B[upper.tri(B,diag=TRUE)])^2)
  s
  
}

slttestperm <- function( B )
{
    #Hack of JONAS PETERS 2013
    #
    # slttestperm - tests if we can permute B to strict lower triangularity
    #
    # If we can, then we return the permutation in p, otherwise p=0.
    #
    
    # Dimensionality of the problem
    n <- nrow(B)    
    
    # This will hold the permutation
    p <- c()
    
    # Remaining nodes
    remnodes <- 1:n
    
    # Remaining B, take absolute value now for convenience
    Brem <- abs(B)
    # Select nodes one-by-one
    for(ii in 1:n)
    {
        # Find the row with all zeros
        #therow = find(sum(Brem,2)<1e-12);
        if(length(Brem) > 1)
        {
            rowS <- rowSums(Brem)
        } else
        {
            rowS <- Brem
        }
        therow <- which(rowS < 1e-12)
        
        # If empty, return 0
        #if isempty(therow),
        if(length(therow) == 0)
        {
            p <- 0
            return(p)    
        }
        # If more than one, arbitrarily select the first 
        therow <- therow[1]
        
        # If we made it to the end, then great!
        if(ii==n)
        {
            p <- c(p,remnodes)
            return(p)
        }
        
        # Take out that row and that column
        inds <- which((1:(n-ii+1)) != therow)
        Brem <- Brem[inds,inds]
        ### CHECK!!!!
        
        # Update remaining nodes
        p <- c(p,remnodes[therow])
        remnodes <- remnodes[inds]
        
    }
}

sqrtm <- function( A ) {

  e <- eigen(A)
  V <- e$vectors
  B <- V %*% diag(sqrt(e$values)) %*% t(V)
  B
  
}

tridecomp <- function( W, choice='qr' ) {

  # SYNTAX:
  # res <- tridecomp( W, choice )
  # QR, RQ, QL, or LQ decomposition specified by
  # choice = 'qr', 'rq', 'ql', or 'lq' respectively
  #
  # res$A is the first matrix and res$B is the second
  #
  # Based on MATLAB code kindly provided by Matthias Bethge
  # Adapted for R by Patrik Hoyer

  m <- nrow(W)
  n <- ncol(W)
  Jm <- matrix(0,m,m)
  Jm[m:1,] <- diag(m)
  Jn <- matrix(0,n,n)
  Jn[n:1,] <- diag(n)

  res <- list()

  if (choice == 'qr') {
    r <- qr(W)
    res$A <- qr.Q(r)
    res$B <- qr.R(r)
  } else if (choice == 'lq') {
    r <- qr(t(W))
    res$A <- t(qr.R(r))
    res$B <- t(qr.Q(r))
  } else if (choice == 'ql') {
    r <- qr(Jm %*% W %*% Jn)
    res$A <- Jm %*% qr.Q(r) %*% Jm
    res$B <- Jm %*% qr.R(r) %*% Jn
  } else if (choice == 'rq') {
    r <- qr(Jn %*% t(W) %*% Jm)
    res$A <- t(Jn %*% qr.R(r) %*% Jm)
    res$B <- t(Jn %*% qr.Q(r) %*% Jn)
  }

  res
  
}

uselingam <- function( X, verbose = FALSE ) {
  temp <- estimate( X, verbose )
  res <- prune( X, temp$k, verbose )
  res
}
