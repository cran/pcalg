pag2edge <- function(P){
  ## Purpose: Constructs a matrix which contains identifiable parental and
  ## non-parental relations in the Markov equivalence class represented
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
  ## ----------------------------------------------------------------------
  ## Value:
  ## - edge: matrix containing for each ordered pair of nodes whether a
  ##   parental relation, non-parental relation, or neither was identified
  ## ----------------------------------------------------------------------
  ## Author: Joris Mooij, Date: Nov 2020

  p<-dim(P)[1]
  edge<-matrix(0,p,p)
  rownames(edge)<-rownames(P)
  colnames(edge)<-colnames(P)
  for( i in seq_len(p) ) {
    for( j in seq_len(p) ) { 
      if( P[i,j] == 0 ) { # i  j
        ## Proposition 7 in Mooij & Claassen (2020)
        edge[i,j] = -1
        edge[j,i] = -1
      } else if( P[i,j] == 2 && P[j,i] == 3 ) { # i --> j
        ## Proposition 8 in Mooij & Claassen (2020)
        if( !visibleEdge(P,i,j) ) { # case 1: no possibly directed path from i to j that avoids edge i-->j
          # copy P but remove edge i --> j
          PP <- P
          PP[i,j] <- 0
          PP[j,i] <- 0
          if( !(j %in% searchAM(PP,i,type="pde")) ) {
            ## no DMG in PP contains a directed path i --> .. --> j
            edge[i,j] = 1
          } else {
            ## some DMG in PP may contain a directed path i --> .. --> j and therefore i --> j is possibly-indirect
            edge[i,j] = 0
          }
        } else { # case 2: i-->j is definitely visible
          edge[i,j] = 1 # optimistic start
          for( k in seq_len(p) ) {
            if( k != i && k != j ) {
              if( (P[i,k] == 1 | P[i,k] == 2) & (P[k,i] == 1 | P[k,i] == 3) ) { ## i (o-)--(o>) k
                if( (P[k,j] == 1 | P[k,j] == 2) & (P[j,k] == 1 | P[j,k] == 3) ) { ## k (o-)--(o>) j
                  # copy P but modify k o-> j into k --> j
                  Pmod<-P
                  Pmod[j,k]<-3
                  Pmod[k,j]<-2
                  if( !visibleEdge(Pmod,k,j) ) {
                    edge[i,j] = 0 # we were too optimistic
                    break
                  }
                }
              }
            }
          }
        }
        edge[j,i] = -1
      } else if( P[i,j] == 2 && P[j,i] == 2 ) { # i <-> j
        ## Proposition 7 in Mooij & Claassen (2020)
        edge[i,j] = -1
        edge[j,i] = -1
      } else if( P[i,j] == 2 && P[j,i] == 1 ) { # i o-> j
        ## Proposition 7 in Mooij & Claassen (2020)
        edge[i,j] = 0
        edge[j,i] = -1
      } else if( P[i,j] == 1 && P[j,i] == 1 ) { # i o-o j
        ## Proposition 7 in Mooij & Claassen (2020)
        edge[i,j] = 0
        edge[j,i] = 0
      }
    }
  }

  return(edge)
} ## {pag2edge}
