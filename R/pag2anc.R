pag2anc <- function(P,verbose=FALSE) {
  ## Purpose: Constructs a matrix which contains identifiable ancestral and
  ## non-ancestral relations in the Markov equivalence class represented
  ## by directed partial ancestral graph P.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - P: adjacency matrix of type amat.pag
  ##        P[i,j] = 0 iff no edge btw i,j
  ##        P[i,j] = 1 iff i *-o j
  ##        P[i,j] = 2 iff i *-> j
  ##        P[i,j] = 3 iff i *-- j
  ##      P has to be directed, i.e., it should not contain undirected
  ##      (x --- y) or circle-tail (x o-- y) edges
  ## - verbose: gives more verbose output if TRUE
  ## ----------------------------------------------------------------------
  ## Value:
  ## - anc: matrix containing for each ordered pair of nodes whether an
  ##   ancestral relation, non-ancestral relation, or neither
  ##   was identified
  ## ----------------------------------------------------------------------
  ## Author: Joris Mooij, Date: Nov 2020

  p <- dim(P)[1]
  anc <- matrix(0,p,p)
  rownames(anc)<-rownames(P)
  colnames(anc)<-colnames(P)
  for( a in seq_len(p) ) {
    for( c in seq_len(p) ) {
      if( a == c ) {
        anc[a,c] <- 1
      } else { ## Proposition 4 in Mooij & Claassen (2020)
        ## check whether a is an identifiable cause of c
        if( a %in% searchAM(P,c,type="an") )
          anc[a,c] <- 1
        else {
          ## search with a breitensuche two minimal uncovered circle paths

          ## find all x s.t. a (o-)--(o>) x
          indX <- which((P[a, ] == 1 | P[a, ] == 2) &
                        (P[, a] == 1 | P[, a] == 3), arr.ind = TRUE)
          indX <- setdiff(indX, c)
          done <- FALSE
          if (length(indX >= 2)) {
            if( verbose )
              cat('Starting search for upd paths from ',a,' to ',c,' with first nodes ',indX,'\n')
            counterX1 <- 0
            while (counterX1 < length(indX) && !done) {
              counterX1 <- counterX1 + 1
              first.pos <- indX[counterX1]
              indX2 <- setdiff(indX, first.pos)
              counterX2 <- 0
              while (counterX2 < length(indX2) && !done) {
                counterX2 <- counterX2 + 1
                sec.pos <- indX2[counterX2]
                if( P[first.pos,sec.pos] == 0 ) { # first.pos and sec.pos are not adjacent
                  if( verbose )
                    cat('  pair {',first.pos,',',sec.pos,'} is promising, searching...\n')
                  t1 <- minUncovPdPath(p, P, a, first.pos, c, unfVect = c(), verbose = FALSE)
                  if( verbose )
                    cat('    first path: ',t1,'\n')
                  if (length(t1) > 1) { # otherwise, can skip next minUnc..()
                    t2 <- minUncovPdPath(p, P, a, sec.pos, c, unfVect = c(), verbose = FALSE)
                    if( verbose )
                      cat('    second path: ',t2,'\n')
                    if (length(t2) > 1 && first.pos != sec.pos) {
                      ## we found 2 uncovered pd paths
                      done <- TRUE
                    }
                  }
                }
              } #  # while ( counterX2 .. )
            }
          }
          if( done )
            anc[a,c] <- 1
        }
      } ## Proposition 5 in Mooij & Claassen (2020)
      ## check whether c is not a possible descendant of a
      if( !(c %in% searchAM(P,a,type="pde")) ) {
        if( anc[a,c] == 1 )
          stop('Inconsistency in pag2anc!')
        anc[a,c] <- -1
      }
    } # for c
  } # for a
  return(anc)
} ## {pag2anc}
