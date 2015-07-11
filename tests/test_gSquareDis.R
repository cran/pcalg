library(pcalg)

cat("doExtras:", (doExtras <- pcalg:::doExtras()), "\n")
showProc.time <- local({
    pct <- proc.time()
    function() { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	cat('Time elapsed: ', (pct - ot)[1:3],'\n')
    }
})

load("discr100k.rda")# in directory tests/ i.e., typically *not* installed
str(dat)# 100'000 x 5 data frame of integers
sapply(dat, table)
nlev.dat <- sapply(dat, function(v) nlevels(factor(v)))
suffStat <- list(dm = dat, nlev = nlev.dat, adaptDF = TRUE)
dat. <- data.matrix(dat)

showProc.time()

## Tests for   |S| = 0  and  |S| = 1

pv.5 <- c(
    ## Collider
    ##        x y  S
    disCItest(1,2,NULL,suffStat),# 1 : Must be !=0
    disCItest(1,2, 3,  suffStat),# 2 : 0
    disCItest(5,2, 4,  suffStat),# 3 : 0
    disCItest(5,2,NULL,suffStat),# 4 : != 0
    ## Non-collider
    disCItest(4,3, 2,  suffStat),# 5 : !=0
    disCItest(4,3,NULL,suffStat),# 6 : 0
    ## Neighbours
    disCItest(4,2,NULL,suffStat),# 7 : 0
    disCItest(1,3,NULL,suffStat),# 8 : 0
    disCItest(1,5,NULL,suffStat))# 9 : !=0

stopifnot(all.equal(pv.5,
                    c(0.5150866,
                      0,
                      0,
                      0.289990516,
                      ## Non-collider
                      0.818905483,
                      0,
                      # Neighbours
                      0,
                      0,
                      0.774158263), tol = 4e-9))

showProc.time()

if(!doExtras) ## smaller data set
    dat. <- dat.[1:4096,]

##  |S| = 2 :
pv2.5 <- c(
    gSquareDis(3,4, S = c(1,2), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,4, S = c(1,3), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,3, S = c(1,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,3, S = c(1,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,4, S = c(2,3), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,3, S = c(2,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,3, S = c(2,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,2, S = c(3,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,2, S = c(3,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,2, S = c(4,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),

    gSquareDis(3,5, S = c(1,2), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,5, S = c(1,3), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,5, S = c(1,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,4, S = c(1,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,5, S = c(2,3), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,5, S = c(2,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,4, S = c(2,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,5, S = c(3,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,4, S = c(3,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,3, S = c(4,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE))

showProc.time()

pv2.5.EXP <- {
    if(doExtras)
        c(0.98726, 0,       0, 0, 0.77097, 0,       0,       0, 0, 0.98974,
          0.95964, 0.99966, 0, 0, 0.99685, 0.99999, 0.42513, 1, 0.11179, 0)
    else
        c(0.73219, 0, 0.0003349, 0, 0.90758, 0, 0, 0.98106, 0.46011, 1,
          0.93245, 0.99998, 8e-07, 0, 0.99797, 0.99933, 0.86036, 0.99988, 0.84859, 0)
}


(eq52 <- all.equal(pv2.5, pv2.5.EXP, tol = 1e-5))
if(!isTRUE(eq52)) stop("gSquareDis(|S| = 2) not ok:", eq52)


##  |S| = 3 :
pv3.5 <- c(
    gSquareDis(1,2, S = c(3,4,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,3, S = c(2,4,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,4, S = c(2,3,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,5, S = c(2,3,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,3, S = c(1,4,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,4, S = c(1,3,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,5, S = c(1,3,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(3,4, S = c(1,2,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(3,5, S = c(1,2,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(4,5, S = c(1,2,3), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE))
showProc.time()

pv3.5.EXP <- {
    if(doExtras)
        c(0,       0, 0.95861, 1, 0,       0, 0,       0.99204, 1, 0)
    else
        c(0.99999, 0, 0.84316, 1, 0.15977, 0, 0.83432, 0.66329, 1, 0)
}

eq53 <- all.equal(pv3.5, pv3.5.EXP, tol=6e-6)
if(!isTRUE(eq53)) stop("gSquareDis(|S| = 3) not ok:", eq53)

## FIXME: Extend data, and do (few) tests with  |S| in {4, 5, 6, 7} !!
