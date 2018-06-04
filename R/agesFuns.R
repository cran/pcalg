ages <- function(data, lambda_min = 0.5*log(nrow(data)), labels=NULL, fixedGaps = NULL, adaptive = c("none","vstructures", "triples"), maxDegree = integer(0), verbose = FALSE, ...){

  ## Initialization of the parameters
  p <- ncol(data)
  if(!is.data.frame(data)){
    if(is.null(colnames(data))){
      colnames(data) <- as.character(1:p)
    }
    data <- as.data.frame(data)
  }
  if(!is.null(labels)){
    names(data) <- labels
  }
  lambda <- lambda_min
  lambda_input <- lambda_min
  dataArray <- data
  i <- 1
  
  ## Initialize the empty graph
  nodes <- colnames(dataArray)
  in.edges <- replicate(length(nodes), integer(0))
  score <- new("GaussL0penObsScore", data=dataArray,lambda=lambda)
    startGraph <- new("EssGraph", nodes = nodes, in.edges = in.edges, score = score)

    ## Store the value for the first edge addition
  globalscore <- score$global.score(startGraph$repr())
  bool <- startGraph$greedy.step(direction = "forward", verbose=verbose, fixedGaps = fixedGaps, adaptive = adaptive, maxDegree = maxDegree, ...)
  lambda <- c(lambda, (score$global.score(startGraph$repr())-globalscore)+lambda_input)
  lengthlambda <- length(lambda)
  s <- lambda[lengthlambda]
  bool <- TRUE
  CPDAGToAdd <- as(startGraph,"matrix")
  colnames(CPDAGToAdd) <- rownames(CPDAGToAdd) <- nodes
  #CPDAGsList <- list(as(startGraph,"matrix"))
  CPDAGsList <- list(CPDAGToAdd)
  
  ## We now have a CPDAG containing one edge. This CPDAG is stored in CPDAGsList and
  ## its respective lambda value in lambda.
  ## At the end CPDAGsList and lambda should have the same length and
  ## the values in lambda should correspond to the CPDAGs found by GES
  ##(only forward and backward phase without iterations!) with these values.
  
  ## TODO: bool=FALSE should never occur now, can delete it (double check first)
  while((lambda[1] < s) & bool){
    
    ForBackwardInEdges <- startGraph$.in.edges
    globalscore <- score$global.score(startGraph$repr())
    
    bool <- startGraph$greedy.step(direction = "forward", verbose=verbose, fixedGaps = fixedGaps, adaptive = adaptive, maxDegree = maxDegree, ...)
    PotentialNextLambda <- score$global.score(startGraph$repr())-globalscore + lambda_input
    
    
    ## If the next lambda candidate is not smaller than the actual penalty parameter
    ## the forward phase is not over yet.
    ## Hence, we simply continue
    if(PotentialNextLambda-lambda[lengthlambda]< -1e-8){
      
      ## Now we have the next penalty candidate. For all values between lambda[lengthlambda] and PotentialNextLambda the forward phase returns the same CPDAG
      ## We start the backward phase with the smallest value, i.e. PotentialNextLambda and then increase it
      ## (The backward phase deletes more edges for larger penalties)
      
      scoreBackward <- new("GaussL0penObsScore", data=dataArray,lambda=PotentialNextLambda)
      BackwardGraph <- new("EssGraph", nodes = nodes, in.edges = ForBackwardInEdges, score = scoreBackward)
      Bbefore <- score$global.score(BackwardGraph$repr())
      
      boolBackward <- BackwardGraph$greedy.step(direction = "backward", verbose=verbose, fixedGaps = fixedGaps, adaptive = adaptive, maxDegree = maxDegree, ...)
      Bafter <- score$global.score(BackwardGraph$repr())
      
      
      
      while(boolBackward){
        
        Bbefore <- score$global.score(BackwardGraph$repr())
        
        boolBackward <- BackwardGraph$greedy.step(direction = "backward", verbose=verbose, fixedGaps = fixedGaps, adaptive = adaptive, maxDegree = maxDegree, ...)
        Bafter <- score$global.score(BackwardGraph$repr())
        
        
      }
      
      ## This is the end of the backward phase for the smallest lambda
      ## We then have to check if more edges will be deleted with a larger lambda (up to s)
      ## This graph should be saved with lambda Potentialnextlambda
      ## In this process we may find more lambdas and CPDAGs
      ## We temporarely save them in NewLambdas and ListOfCPDAGs, respectively
      ListOfCPDAGs <- list(as(BackwardGraph,"matrix"))
      ListOfLambdasFromPotentialNext <- list(c(0))
      
      scoreBackward <- new("GaussL0penObsScore", data=dataArray,lambda=s)
      BackwardGraph <- new("EssGraph", nodes = nodes, in.edges = BackwardGraph$.in.edges, score = scoreBackward)
      
      Bbefore <- score$global.score(BackwardGraph$repr())
      
      boolBackward <- BackwardGraph$greedy.step(direction = "backward", verbose=verbose, fixedGaps = fixedGaps, adaptive = adaptive, maxDegree = maxDegree, ...)
      Bafter <- score$global.score(BackwardGraph$repr())
      
      deltaBackward <-  s - (Bafter - Bbefore + s) + lambda_input
      
      ## If boolBackward is FALSE we did not delete anything, hence we can go on.
      ## TODO: the second part of the if statement is always true now, can delete it (but first double check)
      if(boolBackward & abs(deltaBackward - s)>1e-8){
        
        NewLambdas <- PotentialNextLambda
        
        while(boolBackward){
          
          if(!boolBackward){
            next
          }
          
          ## For each additional deletion we need to compute the lambda value to obtain it.
          ## Also, we need to check if other deletions are possible with the same lambda (previously done with --- < -1e-8 and >1e-8)
          deltaBackward <- s - (Bafter - Bbefore + s) + lambda_input
          
          s_prov <- deltaBackward
          scoreBackward_prov <- new("GaussL0penObsScore", data=dataArray,lambda=s_prov)
          BackwardGraph_prov <- new("EssGraph", nodes = nodes, in.edges = BackwardGraph$.in.edges, score = scoreBackward_prov)
          
          boolBackward_prov <- boolBackward
          
          while(boolBackward_prov){
            ## We do here the further deletions
            Bbefore <- score$global.score(BackwardGraph_prov$repr())
            
            boolBackward_prov <- BackwardGraph_prov$greedy.step(direction = "backward", verbose=verbose, fixedGaps = fixedGaps, adaptive = adaptive, maxDegree = maxDegree, ...)
            Bafter <- score$global.score(BackwardGraph_prov$repr())
          }
          
          ## We are now done, we have a new lambda and a new CPDAG. We need to save it.
          
          BackwardGraph$.in.edges <- BackwardGraph_prov$.in.edges
          
          
          NewLambdas <- c(NewLambdas, deltaBackward)
          PotentialNextLambda <- NewLambdas[1]
          ListOfCPDAGs <- append(ListOfCPDAGs,list(as(BackwardGraph,"matrix")))
          
          
          Bbefore <- score$global.score(BackwardGraph$repr())
          
          boolBackward <- BackwardGraph$greedy.step(direction = "backward", verbose=verbose, fixedGaps = fixedGaps, adaptive = adaptive, maxDegree = maxDegree, ...)
          Bafter <- score$global.score(BackwardGraph$repr())
          
        }
        
        ## Now we add ListOfCPDAGs and NewLambdas to the correspective lists
        
        if(length(NewLambdas)>1){
          for(q in length(NewLambdas):2){
            lambda <- c(lambda,NewLambdas[[q]])
            lengthlambda <- length(lambda)
            s <- lambda[lengthlambda]
            CPDAGToAdd <- ListOfCPDAGs[[q]]
            colnames(CPDAGToAdd) <- rownames(CPDAGToAdd) <- nodes
            #CPDAGsList <- list(as(startGraph,"matrix"))
            CPDAGsList <- append(CPDAGsList,list(CPDAGToAdd))
          }
        }
        lambda <- c(lambda,PotentialNextLambda)
        lengthlambda <- length(lambda)
        s <- lambda[lengthlambda]
        CPDAGToAdd <- ListOfCPDAGs[[1]]
        colnames(CPDAGToAdd) <- rownames(CPDAGToAdd) <- nodes
        CPDAGsList <- append(CPDAGsList,list(CPDAGToAdd))
        
      }else{
        if(PotentialNextLambda-lambda[lengthlambda]< -1e-8){
          lambda <- c(lambda,PotentialNextLambda)
          lengthlambda <- lengthlambda + 1
          s <- lambda[lengthlambda]
          CPDAGToAdd <- ListOfCPDAGs[[1]]
          colnames(CPDAGToAdd) <- rownames(CPDAGToAdd) <- nodes
          CPDAGsList <- append(CPDAGsList,list(CPDAGToAdd))
        }
      }
    }
  }
  lambda <- c(lambda,lambda[1])
  lambda <- lambda[-1]
  
  bool.final <- T
  while(bool.final){
  bool.final <- startGraph$greedy.step(direction = "forward", verbose=verbose, fixedGaps = fixedGaps, adaptive = adaptive, maxDegree = maxDegree, ...)
  }

  bool.final.bw <- T
  while(bool.final.bw){
    bool.final.bw <- startGraph$greedy.step(direction = "backward", verbose=verbose, fixedGaps = fixedGaps, adaptive = adaptive, maxDegree = maxDegree, ...)
  }

  CPDAGToAdd <- as(startGraph,"matrix")
  colnames(CPDAGToAdd) <- rownames(CPDAGToAdd) <- nodes
  CPDAGsList <- append(CPDAGsList,list(CPDAGToAdd))

  
  resFinal <- AggregationMethodsSAGES(CPDAGsList)[[1]]
  
  ### Construct the EssGraph object to return so that is similar to the output of ges()
  
  in.edges <- list(which(resFinal[,1]==1))
  if(p>1){
  for(i in 2:p){
    in.edges <- append(in.edges,list(which(resFinal[,i]==1))) 
  }
  }
  
  essgraph <- new("EssGraph", nodes, in.edges)
  
  repr <- pdag2dag(as(essgraph,"graphNEL"))
  stopifnot(repr$success)
  
  repr <- repr$graph
  
  # return(list(CPDAGsList,lambda,resFinal,resFinal500,resFinal1000,resFinal5000,resFinal10000,resFinal50000,ges.return,ess,re))
  Ages <- list(essgraph,repr,CPDAGsList,lambda)
  names(Ages) <- c("essgraph","repr","CPDAGsList","lambda")
  return(Ages)
  
  # return(list(essgraph,repr,CPDAGsList,lambda))
  }

AggregationMethodsSAGES <- function(ListOfGraphs, startconsideringfrom = 2){
  
  len <- length(ListOfGraphs)
  
  skeleton <- ((ListOfGraphs[[len]] + t(ListOfGraphs[[len]]))>0)*1
  dir <- ((ListOfGraphs[[len]] - t(ListOfGraphs[[len]]))>0)*1
  
  for(j in 2:(len-startconsideringfrom+1)){
    i <- len-j+1
    sk2 <- ((ListOfGraphs[[i]] + t(ListOfGraphs[[i]]))>0)*1
    if(sum((skeleton - sk2)==-1)>0){
      next
    }
    dir2 <- ((ListOfGraphs[[i]] - t(ListOfGraphs[[i]]))>0)*1
    dir.prov <- ((dir + dir2 - t(dir))>0)*1
    pdag.prov <- ((skeleton - t(dir.prov))>0)*1
    if(pdag2dag(as(pdag.prov,"graphNEL"))[[2]]){
      
      pdag <- applyOrientationRules(pdag.prov)
      dir <- ((pdag - t(pdag))>0)*1
      
    }
    
  }
  
  res <- ((skeleton - t(dir))>0)*1
  
  return(list(res))
  
}
