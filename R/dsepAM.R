dsepAM <- function(X,Y, S = NULL, amat, verbose=FALSE)
{
  ## Purpose: Tests whether sets X and the set Y are m-separated given set S
  ##          in a MAG.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - X,Y,S: vectors with node indices (row/column indices, NOT node names)
  ## - amat: adjacency matrix of the MAG
  ##         amat[i,j] = 0 iff no edge btw i,j
  ##         amat[i,j] = 2 iff i *-> j
  ##         amat[i,j] = 3 iff i *-- j
  ## ----------------------------------------------------------------------
  ## Value:
  ## Boolean decision
  ## ----------------------------------------------------------------------
  ## Authors: Markus Kalisch, Joris Mooij
##   suppressMessages(require(RBGL))
##  suppressMessages(require(igraph))
##  suppressMessages(require(graph))

  ## Check whether amat is a valid MAG
  if (dim(amat)[1] != dim(amat)[2]) 
    stop("dsepAM: amat should be a square matrix")
  if (length(setdiff(unique(as.vector(amat)),c(0,2,3))) > 0) 
    stop("dsepAM: amat contains invalid edge marks")
  p <- dim(amat)[1]
  my.nodes <- colnames(amat)
  if( is.null(my.nodes) )
    my.nodes <- as.character(seq_len(p))

  ## build node union of X,Y,S
  nodeUnion <- if(length(S) > 0) unique(c(X,Y,S)) else unique(c(X,Y))

  ## take subgraph of anteriors of nodeUnion
  ant.set <- searchAM(amat,nodeUnion,"ant")
  amatS <- matrix(0,length(ant.set),length(ant.set))
  amatS[1:length(ant.set),1:length(ant.set)] <- amat[ant.set,ant.set]

  ## find all <-> edges in the subMAG
  indA <- which((amatS == 2 & t(amatS) == 2), arr.ind = TRUE)
  amatS.bidirected <- matrix(0,length(ant.set),length(ant.set))
  amatS.bidirected[indA] <- 1
  gS.bidirected <- as(amatS.bidirected,"graphNEL")
  ## Get all districts in the subMAG
  districts <- connectedComp(gS.bidirected)

  ## Moralize
  amatSM <- matrix(0,length(ant.set),length(ant.set))
  for( d in seq_len(length(districts)) ) {
    distset <- as.numeric(districts[[d]])
    paset <- searchAM(amatS,distset,"pa")
    padistset <- unique(c(paset,distset))
    # add fully connected subgraph on district + parents
    amatSM[padistset,padistset] <- 1
  }
  amatSM <- amatSM - diag(diag(amatSM))
  rownames(amatSM) <- my.nodes[ant.set]
  colnames(amatSM) <- my.nodes[ant.set]
  gSM <- as(amatSM,"graphNEL")

  # do separation test in moralized graph
  separates(my.nodes[X],my.nodes[Y],my.nodes[S],gSM)
} ## {dsepAM}


dsepAMTest <- function(x,y, S = NULL, suffStat) 
{
  ## suffStat$g: Adjacency matrix of MAG
  ## suffStat$verbose: Verbosity level
  ## Returns "p-value" P =
  ##	0: d-connected
  ##  1: d-separated

  stopifnot(length(x) == 1 && length(y) == 1)
  if( x == y || x %in% S || y %in% S) {
    0
  } else {
    stopifnot(is(g <- suffStat$g, "matrix"))
    stopifnot(dim(g)[1] == dim(g)[2])
    if( is.null(suffStat$verbose) )
      verbose=FALSE
    else
      verbose=suffStat$verbose
    as.numeric(dsepAM(x,y,S,g,verbose))
  }
} ## {dsepAMTest}
