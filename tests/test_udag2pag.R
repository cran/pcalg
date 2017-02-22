library(pcalg)

m <- matrix(FALSE, 3,3)
tmp <- list(NULL, NULL, NULL)

(pag0 <- udag2pag(m, sepset = rep(tmp, 3)))
stopifnot(is.numeric(pag0))

