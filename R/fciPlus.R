##' Find possible Dsep links, given the adjacency matrix,
##' as the edges that satisfy the pattern described in Lemma 4.
##'
##' The edge X <-> Y is a PosDsep link in G+ if there exist nodes, U,V such that
##'  U <-> X <-> Y <-> V in G+, and U and V are not adjacent and paths V.. -> X
##' U.. -> Y exist in G+ and they do not go against arrowhead
##' @title Find POSsible D-Separation LINKS
##' @param m adjacency matrix
##' @param fixedEdges matrix with edges that are guaranteed to be present
##' @return a (possibly empty) data frame with columns "x" and "y"
##'	 where  \code{x[i] *-* y[i]} is a possible Dsep link for all i.
PosDsepLinks <- function(m, fixedEdges=NULL, verbose=TRUE)
{
  stopifnot(1 <= (p <- ncol(m)))
  i.p <- 1:p
  NAAA <- vector("list", p)
  x <- y <- integer()
  ## for all pairs  1 <= j < i <= p
  for(i in i.p[-1L]) {
    ## nodes that have a  "not against arrowhead" path to i :
    NAAA[[i]] <- NAAA_i <- NAA_ancestors(i, m)
    for (j in 1:i) {
      if( is.null(fixedEdges) | !fixedEdges[i,j] ) {
        ## first find a bidirected i <-> j edge in G+
        if (m[i,j] == 2 && m[j,i] == 2) {
          u <- v <- integer()
          ## Then find bidirected edges u <-> i and j <-> v
          for (k in i.p) if (k != i && k != j)
          {
            if (m[i,k] == 2 && m[k,i] == 2)
              u <- union(u,k)
            if (m[j,k] == 2 && m[k,j] == 2)
              v <- union(v,k)
          }

          ## find nodes u (v) that have both a bidirected edge with i (j) and
          ## a not against arrowhead path to j (i)
          u <- intersect(u, NAAA[[j]]) ## as j <= i, this is already computed
          v <- intersect(v, NAAA_i)

          ## check if there is a pair of nodes in (u,v) that is not adjacent
          found.i.j <- FALSE
          for(k in seq_along(u)) {
            u.k <- u[k]
            for(r in seq_along(v))
              if (u.k != v[r] && m[u.k,v[r]] == 0) {
                found.i.j <- TRUE
                break
              }
            if(found.i.j)
              break
          }

          ## found.i.j is true iff we found a pair of nodes fulfilling all 3 conditions.
          ## In that case, the node pair (i,j)  is added to the  PosDsepLinks set:
          if (found.i.j) {
            if(verbose) cat(sprintf("Added PosDsepLink  %2d *-* %2d\n", i,j))
            x <- c(x, i)
            y <- c(y, j)
          }
        }
      }
    }
  }

  ## returns data frame with possible Dsep links x[i] *-* y[i] for all i
  data.frame(x = x, y = y)
}


## The final fciplus function, it uses all other functions and returns
## an adjacency matrix; however, it does not go through the 10 orientation rules.
##
## the input for the function is a pc fitted graph - pc.fit
## as well as a sufficient statistic and an independence test
fciplus.intern <- function(pc.fit, alpha = 0.01, suffStat, indepTest, fixedEdges = NULL, verbose=TRUE)
{
  sepsets <- pc.fit@sepset
  cpdag <- pc.fit@graph
  m <- trafoCPDmat(as(cpdag, "matrix"))

  ## first run the augment graph function to orient invariant
  ## arrowheads according to Lemma 2. (1)
  mat <- AugmentGraph(m, suffStat, sepsets, indepTest, alpha)
  ## Check if there are any possible Dsep links
  link <- PosDsepLinks(mat, fixedEdges)

  ## if there are possible Dsep links it should be checked
  ## which of these are actually edges that should be removed
  ## an edge X - Y is a Dsep link if X and Y are Dseparated by the hierarchy
  ## of the parents of the nodes X and Y (excluding X and Y)

  ## a  counter to help us go through the data frame of all possible Dsep lins
  count <- 1L
  while (count <= length(link$x)) {

    x <- link$x[count]
    y <- link$y[count]

    ## basex and basey are vectors of the nodes neighboring nodes x,y
    basex <- setdiff(which(mat[x,] != 0), y)
    basey <- setdiff(which(mat[y,] != 0), x)
    all <- union(basex,basey)

    ## to determine whether a posDsep link is an actual Dsep link
    ## it is necessarry to go through all the subsets of the sets basex and
    ## basey (since we don't know the parents) and build a hierarchy
    ## for each combination of those subsets

    ## since it is complicated to go through all the subets of both basex and basey,
    ## I decided to combine them into one set and go through subsets of that set
    ## but only building a hierarchy using a subset that contatins neighbors of
    ## both x and y

    ## bound is the number of neighbors of x and y
    bound <- length(all)

    ## flag is an indicator that will be given the value TRUE if
    ## an actual Dsep link is discovered so we can use it to break
    ## from the for cycle and return to the while cycle
    flag <- FALSE
    for(i in seq_len(bound)) {
      if (flag) break
      S <- seq_len(i)
      exit <- FALSE
      while (!exit && !flag)
      {
        a.S <- all[S]
        ## check if there are more than 1 and less than k neighbors of x and y in the subset
        if (1 <= (l.S.bx <- length(intersect(a.S, basex))) && l.S.bx <= i &&
            1 <= (l.S.by <- length(intersect(a.S, basey))) && l.S.by <= i)
        {
          ## build a hierarchy
          potential <- HIE(union(c(x,y), a.S), sepsets)
          ## remove x and y from the hierarchy
          potential <- setdiff(potential, c(x,y))

          ## if this set is actually a separating set it should be added to
          ## the list of sepsets and we should begin the process again (AugmentGraph, PosDsepLink)
          if (indepTest(x,y, potential, suffStat) > alpha)
          {
            ## find the minimal sepset
            mindsep <- MinimalDsep(x,y,potential,suffStat,indepTest)

            ## update sepsets
            if (is.null(mindsep)) {
              sepsets[[x]][y] <- list(NULL)
            } else {
              sepsets[[x]][[y]] <- mindsep
            }

            ## update matrix by removing the Dsep link
            mat[x,y] <- mat[y,x] <- 0

            cat("Found Dsep Link for x=",x,"y=",y,"sepset=",mindsep,"\n")

            ## orient invariant arrowheads
            mat <- AugmentGraph(mat,suffStat,sepsets,indepTest,alpha)

            ## search for new PosDsepLinks
            link <- PosDsepLinks(mat, fixedEdges)

            ## reset counter and set flag to TRUE
            count <- 1L
            flag <- TRUE
          }
        }

        z <- getNextSet(bound, i, S)
        S <- z$nextSet
        exit <- z$wasLast
      } ## while(!exit ..)
    } ## for(i .. )
    ## if we've exited the for cycle and flag is still false it means x-y is not
    ## a Dseplink so we should move on to the next posDseplink
    if (!flag)
      count <- count+1L
  } ## while(count <= ..)
  ## return the adjusted adjacency matrix and sepset
  list(mat = mat, sepset = sepsets)
} ## {fciplus.intern}


fciPlus <- function(suffStat, indepTest, alpha, labels, p, verbose=TRUE, selectionBias = TRUE, jci = c("0","1","12","123"), contextVars = NULL)
{
  ## Purpose: Perform FCI+-Algorithm, i.e., estimate PAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - suffStat, indepTest: info for the tests
  ## - p: number of nodes in the graph
  ## - alpha: Significance level of individual partial correlation tests
  ## - verbose: 0 - no output, 1 - detailed output
  ## - labels: names of the variables or nodes
  ## - selectionBias: allow for selection bias (default: TRUE)
  ## - jci: specifies the JCI background knowledge that is used; can be either:
  ##     "0"   no JCI background knowledge (default),
  ##     "1"   JCI assumption 1 only (i.e., no system variable causes any context variable),
  ##     "12"  JCI assumptions 1 and 2 (i.e., no system variable causes any context variable,
  ##           and no system variable is confounded with any context variable),
  ##     "123" all JCI assumptions 1, 2 and 3
  ## - contextVars: subset of variable indices that will be treated as context variables
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  7 Jul 2014, 12:08; update: Joris Mooij, 2020

  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    ## FIXME ---- make function for this and use in skeleton(), pc(), fci(), ....
    ## if(is.matrix(C <- suffStat[["C"]]) && (d <- dim(C))[1] == d[2]) {
    ##   ## we can derive 'p'  *and* 'labels' -- in 99% of cases !
    ##   p <- d[1]
    ##   if(is.null(labels <- colnames(C)))
    ##      labels <- as.character(seq_len(p))
    ## } else {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
    ## }
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }

  ## Check that jci background knowledge is valid
  jci <- match.arg(jci)
  ## Check whether contextVars is valid
  if( !is.null(contextVars) && length(contextVars) > 0 ) {
    if( !is.numeric(contextVars) )
      stop("contextVars has to be a vector of integers in {1,2,..,p}, where p is the number of variables.")
    if( !all(sapply(contextVars, function(i) i == as.integer(i))) )
      stop("contextVars has to be a vector of integers in {1,2,..,p}, where p is the number of variables.")
    if( min(contextVars) < 1 || max(contextVars) > p )
      stop("contextVars has to be a vector of integers in {1,2,..,p}, where p is the number of variables.")
  }
  ## Set fixed edges from JCI assumption 3 if asked for
  fixedEdges<-matrix(FALSE,p,p)
  if( jci == "123" && length(contextVars) > 0 ) {
    fixedEdges[contextVars,contextVars]<-TRUE
  }

  skel <- skeleton(suffStat = suffStat, indepTest = indepTest, alpha = alpha,
                   labels = labels, fixedEdges = fixedEdges, p = p)
  fit1 <- udag2pdagRelaxed(gInput = skel, orientCollider = FALSE)
  fcip <- fciplus.intern(pc.fit = fit1, alpha=alpha, suffStat=suffStat,
                         indepTest=indepTest, fixedEdges=fixedEdges, verbose=verbose)
  fcip$mat <- (fcip$mat != 0) + 0 # forget orientations in augmented graph: FCI orientation rules seem more accurate  

  rules = rep(TRUE,10)
  if( !selectionBias )
    rules[5:7] <- FALSE
  else {
    if( jci != "0" )
      stop( 'The current JCI implementation does not support selection bias (use selectionBias=FALSE instead).' )
  }
  fciplus.amat <- udag2pag(pag = fcip$mat, sepset = fcip$sepset, rules = rules,
                           orientCollider = TRUE, jci = jci, contextVars = contextVars, 
                           verbose = verbose)
  colnames(fciplus.amat) <- rownames(fciplus.amat) <- labels
  new("fciAlgo", amat = fciplus.amat, call = cl, n = integer(0),
      max.ord = integer(0),
      max.ordPDSEP = integer(0),
      n.edgetests = integer(0), n.edgetestsPDSEP = integer(0),
      sepset = list(), pMax = matrix(0,1,1), allPdsep = list())
} ## {fciPlus}
