library(pcalg)

val <- rep(FALSE, 9)
## Test 1
gm <- rbind(c(0,1),
            c(1,0))
colnames(gm) <- rownames(gm) <- LETTERS[1:2]
res <- pdag2allDags(gm, verbose = FALSE)
## plotAllDags(res)
## for (i in 1:2) {cat(paste("c(",paste(res$dags[i,], collapse = ","),"),",sep=""),"\n")}
res.truth <- rbind(c(0,1,0,0), 
                   c(0,0,1,0))
val[1] <- identical(res.truth, res$dags)

## Test 2
gm <- rbind(c(0,0),
            c(0,0))
colnames(gm) <- rownames(gm) <- LETTERS[1:2]
res <- pdag2allDags(gm, verbose = FALSE)
res.truth <- rbind(c(0,0,0,0))
val[2] <- identical(res.truth, res$dags)

## Test 3, non-collider
gm <- rbind(c(0,1,0),
            c(1,0,1),
            c(0,1,0))
colnames(gm) <- rownames(gm) <- LETTERS[1:ncol(gm)]
res <- pdag2allDags(gm, verbose = FALSE)
res.truth <- rbind(c(0,1,0,0,0,1,0,0,0), 
                   c(0,1,0,0,0,0,0,1,0), 
                   c(0,0,0,1,0,0,0,1,0))
val[3] <- identical(res.truth, res$dags)
## for (i in 1:3) {cat(paste("c(",paste(res$dags[i,], collapse = ","),"),",sep=""),"\n")}
## plotAllDags(res)

## Test 4, collider
gm <- rbind(c(0,0,0),
            c(1,0,1),
            c(0,0,0))
colnames(gm) <- rownames(gm) <- LETTERS[1:ncol(gm)]
res <- pdag2allDags(gm, verbose = FALSE)
res.truth <- rbind(c(0,0,0,1,0,1,0,0,0))
val[4] <- identical(res.truth, res$dags)
## plotAllDags(res)

## Test 5, trick question
gm <- rbind(c(0,0,0),
            c(1,0,1),
            c(0,1,0))
colnames(gm) <- rownames(gm) <- LETTERS[1:ncol(gm)]
res <- pdag2allDags(gm, verbose = FALSE)
res.truth <- rbind(c(0,0,0,1,0,0,0,1,0))
val[5] <- identical(res.truth, res$dags)
## plotAllDags(res)

## Test 6,complete
gm <- rbind(c(0,1,1),
            c(1,0,1),
            c(1,1,0))
colnames(gm) <- rownames(gm) <- LETTERS[1:ncol(gm)]
res <- pdag2allDags(gm, verbose = FALSE)
## for (i in 1:6) {cat(paste("c(",paste(res$dags[i,], collapse = ","),"),",sep=""),"\n")}
res.truth <- rbind(c(0,1,1,0,0,1,0,0,0), 
                   c(0,1,1,0,0,0,0,1,0), 
                   c(0,0,1,1,0,1,0,0,0), 
                   c(0,0,0,1,0,1,1,0,0), 
                   c(0,1,0,0,0,0,1,1,0), 
                   c(0,0,0,1,0,0,1,1,0))
val[6] <- identical(res.truth, res$dags)
## plotAllDags(res)

## Test 7, 4 nodes: being really mean
## No DAG possible
gm <- rbind(c(0,1,1,0),
            c(1,0,0,1),
            c(1,0,0,1),
            c(0,1,1,0))
colnames(gm) <- rownames(gm) <- LETTERS[1:ncol(gm)]
res <- pdag2allDags(gm, verbose = FALSE)
val[7] <- is.null(res$dags)
## plotAllDags(res)

## Test 8, 4 nodes
gm <- rbind(c(0,1,1,0),
            c(1,0,0,0),
            c(1,0,0,0),
            c(0,1,1,0))
colnames(gm) <- rownames(gm) <- LETTERS[1:ncol(gm)]
res <- pdag2allDags(gm, verbose = FALSE)
res.truth <- rbind(c(0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,0),
                   c(0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0),
                   c(0,1,0,0,0,0,0,0,1,0,0,0,0,1,1,0))
val[8] <- identical(res.truth, res$dags)
## plotAllDags(res)

## Test 9, 5 nodes; D -> E must always be there
gm <- rbind(c(0,1,1,0,0),
            c(1,0,0,0,0),
            c(1,0,0,0,0),
            c(0,1,1,0,1),
            c(0,0,0,1,0))

colnames(gm) <- rownames(gm) <- LETTERS[1:ncol(gm)]
res <- pdag2allDags(gm, verbose = FALSE)

res.truth <- rbind(c(0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0),
                   c(0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,0),
                   c(0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,0))
val[9] <- identical(res.truth, res$dags)

if (!all(val)) stop("Error in testing pdag2allDags!\n")
