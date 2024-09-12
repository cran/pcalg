library(pcalg)

## possDe
## a -> b -- c
amat <- matrix(c(0,1,0, 0,0,1, 0,1,0), 3,3)
colnames(amat) <- rownames(amat) <- letters[1:3]
## plot(as(t(amat), "graphNEL"))

stopifnot(
  all(possDe(m = amat, x = 1, possible = TRUE, ds = FALSE, type = "pdag") == c(1,2,3)),
  all(possDe(m = amat, x = 1, possible = FALSE, ds = FALSE, type = "pdag") == c(1,2)),
  all(possDe(m = amat, x = 1, y = 2, possible = TRUE, ds = FALSE, type = "pdag") == 1) )

## possAn
## a -- b -> c
amat <- matrix(c(0,1,0, 1,0,1, 0,0,0), 3,3)
colnames(amat) <- rownames(amat) <- letters[1:3]
## plot(as(t(amat), "graphNEL"))

stopifnot(
  all(possAn(m = amat, x = 3, possible = TRUE, ds = FALSE, type = "pdag") == c(1,2,3)),
  all(possAn(m = amat, x = 3, y = 2, possible = TRUE, ds = FALSE, type = "pdag") == 3) )
