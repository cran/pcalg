library(pcalg)

res <- rep(FALSE, 10)
set.seed(123)
g <- pcalg::randomDAG(n = 7, prob = 0.3)
## plot(g)
cpdag <- dag2cpdag(g)
## plot(cpdag)
cpdag.mat <- t(as(cpdag,"matrix")) ## has correct encoding

## test 1: using graph, valid
g1 <- addBgKnowledge(gInput = cpdag, x = 3, y = 5)
m1 <- t(as(g1, "matrix"))
res[1] <- (m1[3,5]==0 & m1[5,3]==1) ## should be true

## test 2: using matrix, valid
m2 <- addBgKnowledge(gInput = cpdag.mat, x = 3, y = 5, verbose = TRUE)
res[2] <- (m2[3,5]==0 & m2[5,3]==1) ## should be true

## test 3: using matrix, invalid
m3 <- addBgKnowledge(gInput = cpdag.mat, x = 6, y = 3, verbose = TRUE)
res[3] <- is.null(m3)

## test 4: using graph, invalid
g4 <- addBgKnowledge(gInput = cpdag, x = 6, y = 3, verbose = TRUE)
res[4] <- is.null(g4)

## test 5: using matrix, invalid
m5 <- addBgKnowledge(gInput = cpdag.mat, x = 1, y = 3, verbose = TRUE)
res[5] <- is.null(m5)

## test 6: using graph, invalid
g6 <- addBgKnowledge(gInput = cpdag, x = 1, y = 3, verbose = TRUE)
res[6] <- is.null(g6)

## test 7: empty background knowledge: Meek rule 1
m7 <- matrix(0, 3,3)
colnames(m7) <- rownames(m7) <- as.character(1:ncol(m7))
r7T <- m7
m7[2,1] <- m7[3,2] <- m7[2,3] <- 1
r7 <- addBgKnowledge(gInput = m7, x = c(), y = c(), verbose = TRUE,
                     checkInput = FALSE)
r7T[2,1] <- r7T[3,2] <- 1
res[7] <- identical(r7,r7T)

## test 8: empty background knowledge: Meek rule 2
m8 <- matrix(0, 3,3)
colnames(m8) <- rownames(m8) <- as.character(1:ncol(m8))
r8T <- m8
m8[1,2] <- m8[2,3] <- m8[3,1] <- m8[3,2] <- 1
r8 <- addBgKnowledge(gInput = m8, x = c(), y = c(), verbose = TRUE,
                     checkInput = FALSE)
r8T[1,2] <- r8T[3,1] <- r8T[3,2] <- 1
res[8] <- identical(r8,r8T)

## test 9: empty background knowledge: Meek rule 3
m9 <- matrix(0, 4,4)
colnames(m9) <- rownames(m9) <- as.character(1:ncol(m9))
r9T <- m9
m9[1,2:4] <- m9[2,1] <- m9[3,c(1,2,4)] <- m9[4,1] <- 1
r9 <- addBgKnowledge(gInput = m9, x = c(), y = c(), verbose = TRUE,
                     checkInput = FALSE)
r9T[1,c(2,4)] <- r9T[2,1] <- r9T[3,c(1,2,4)] <- r9T[4,1] <- 1
res[9] <- identical(r9,r9T)

## test 10: empty background knowledge: Meek rule 4
m10 <- matrix(0, 4,4)
colnames(m10) <- rownames(m10) <- as.character(1:ncol(m10))
r10T <- m10
m10[1,2:4] <- m10[2,c(1,3)] <- m10[3,c(1,4)] <- m10[4,1] <- 1
r10 <- addBgKnowledge(gInput = m10, x = c(), y = c(), verbose = TRUE,
                      checkInput = FALSE)
r10T[1,c(3,4)] <- r10T[2,c(1,3)] <- r10T[3,c(1,4)] <- r10T[4,1] <- 1
res[10] <- identical(r9,r9T)

## final result
stopifnot(all(res))
