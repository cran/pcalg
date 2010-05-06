library(pcalg)

## geht nicht
G <- randomDAG(5, 0.001)
B <- as(G, "matrix") ## geht nicht
t(wgtMatrix(G)) ## geht

G <- randomDAG(5, 0.4)
B <- as(G, "matrix") ## geht nicht
t(wgtMatrix(G)) ## geht


## geht
G2 <- new("graphNEL", nodes = as.character(1:5))
B2 <- as(G2, "matrix")

##- Loesungsvorschlag:
##- ==================  
randomDAG <- function (n, prob, lB = 0.1, uB = 1) 
{
    stopifnot(n >= 2, is.numeric(prob), length(prob) == 1, 0 <= 
        prob, prob <= 1, is.numeric(lB), is.numeric(uB), lB <= 
        uB)
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
