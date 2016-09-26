library(pcalg)

xx <- TRUE
##################################################
## DAG / CPDAG
##################################################
## CPDAG 1: Paper Fig 1
mFig1 <- matrix(c(0,1,1,0,0,0, 1,0,1,1,1,0, 0,0,0,0,0,1,
                  0,1,1,0,1,1, 0,1,0,1,0,1, 0,0,0,0,0,0), 6,6)
type <- "cpdag"
x <- 3; y <- 6
## FIXME: test more than just $gac
## Ver.1: Let gac() return an S3 class, say "GACfit" or "gacFit", with a print() method
##        and (auto)print(.) everywhere below, save *.Rout.save -> output compared: Is ok, as all "discrete"

xx <- xx &  gac(mFig1,x,y, z=c(2,4), type)$gac
xx <- xx &  gac(mFig1,x,y, z=c(4,5), type)$gac
xx <- xx &  gac(mFig1,x,y, z=c(4,2,1), type)$gac
xx <- xx &  gac(mFig1,x,y, z=c(4,5,1), type)$gac
xx <- xx &  gac(mFig1,x,y, z=c(4,2,5), type)$gac
xx <- xx &  gac(mFig1,x,y, z=c(4,2,5,1), type)$gac
xx <- xx & !gac(mFig1,x,y, z= 2,    type)$gac
xx <- xx & !gac(mFig1,x,y, z= NULL, type)$gac

## CPDAG 2: Paper Fig 5a
mFig5a <- matrix(c(0,1,0,0,0, 1,0,1,0,0, 0,0,0,0,1, 0,0,1,0,0, 0,0,0,0,0), 5,5)
type <- "cpdag"
x <- c(1,5); y <- 4
xx <- xx &  gac(mFig5a, x,y, z=c(2,3), type)$gac
xx <- xx & !gac(mFig5a, x,y, z= 2,     type)$gac

## DAG 1 from Marloes' Talk
mMMd1 <- matrix(c(0,1,0,1,0,0, 0,0,1,0,1,0, 0,0,0,0,0,1,
                  0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0),6,6)
type <- "dag"
x <- 1; y <- 3
xx <- xx &  gac(mMMd1, x,y, z=NULL, type)$gac
xx <- xx & !gac(mMMd1, x,y, z= 2, type)$gac
xx <- xx &  gac(mMMd1, x,y, z= 4, type)$gac
xx <- xx & !gac(mMMd1, x,y, z= 5, type)$gac
xx <- xx & !gac(mMMd1, x,y, z= 6, type)$gac
xx <- xx & !gac(mMMd1, x,y, z=c(4,5), type)$gac

## DAG 2 from Marloes' Talk
mMMd2 <- matrix(c(0,1,0,1,0,0, 0,0,0,0,0,0, 0,1,0,0,1,0,
                  0,0,0,0,1,0, 0,0,0,0,0,1, 0,0,0,0,0,0), 6,6)
type <- "dag"
x <- 4; y <- 6
xx <- xx &  gac(mMMd2, x,y, z= 1, type)$gac
xx <- xx &  gac(mMMd2, x,y, z= 3, type)$gac
xx <- xx & !gac(mMMd2, x,y, z= 5, type)$gac
xx <- xx & !gac(mMMd2, x,y, z=c(1,5), type)$gac
xx <- xx &  gac(mMMd2, x,y, z=c(1,2), type)$gac
xx <- xx &  gac(mMMd2, x,y, z=c(1,3), type)$gac
xx <- xx & !gac(mMMd2, x,y, z= 2, type)$gac

##################################################
## PAG
##################################################
mFig3a <- matrix(c(0,1,0,0, 1,0,1,1, 0,1,0,1, 0,1,1,0), 4,4)
mFig3b <- matrix(c(0,2,0,0, 3,0,3,3, 0,2,0,3, 0,2,2,0), 4,4)
mFig3c <- matrix(c(0,3,0,0, 2,0,3,3, 0,2,0,3, 0,2,2,0), 4,4)
mFig4a <- matrix(c(0,0,1,0,0,0, 0,0,1,0,0,0, 2,2,0,3,3,2,
                   0,0,2,0,2,2, 0,0,2,1,0,2, 0,0,1,3,3,0), 6,6)
mFig4b <- matrix(c(0,0,1,0,0,0, 0,0,1,0,0,0, 2,2,0,0,3,2,
                   0,0,0,0,2,2, 0,0,2,3,0,2, 0,0,2,3,2,0), 6,6)
mFig5b <- matrix(c(0,1,0,0,0,0,0, 2,0,2,3,0,3,0, 0,1,0,0,0,0,0, 0,2,0,0,3,0,0,
                   0,0,0,2,0,2,3, 0,2,0,0,2,0,0, 0,0,0,0,2,0,0), 7,7)
type <- "pag"
xx <- xx & !gac(mFig3a, x=2,      y=4, z=NULL,   type)$gac
xx <- xx & !gac(mFig3b, x=2,      y=4, z=NULL,   type)$gac
xx <- xx &  gac(mFig3c, x=2,      y=4, z=NULL,   type)$gac
xx <- xx & !gac(mFig4a, x=3,      y=4, z=NULL,   type)$gac
xx <- xx &  gac(mFig4a, x=3,      y=4, z= 6,     type)$gac
xx <- xx &  gac(mFig4a, x=3,      y=4, z=c(1,6), type)$gac
xx <- xx &  gac(mFig4a, x=3,      y=4, z=c(2,6), type)$gac
xx <- xx &  gac(mFig4a, x=3,      y=4, z=c(1,2,6), type)$gac
xx <- xx & !gac(mFig4b, x=3,      y=4, z=NULL,   type)$gac
xx <- xx & !gac(mFig4b, x=3,      y=4, z= 6,     type)$gac
xx <- xx & !gac(mFig4b, x=3,      y=4, z=c(5,6), type)$gac
xx <- xx & !gac(mFig5b, x=c(2,7), y=6, z=NULL,   type)$gac
xx <- xx &  gac(mFig5b, x=c(2,7), y=6, z=c(4,5), type)$gac
xx <- xx &  gac(mFig5b, x=c(2,7), y=6, z=c(4,5,1), type)$gac
xx <- xx &  gac(mFig5b, x=c(2,7), y=6, z=c(4,5,3), type)$gac
xx <- xx &  gac(mFig5b, x=c(2,7), y=6, z=c(1,3,4,5), type)$gac

## PAG in Marloes' talk
mMMp <- matrix(c(0,0,0,3,2,0,0, 0,0,0,0,1,0,0, 0,0,0,0,1,0,0, 2,0,0,0,0,3,2,
                 3,2,2,0,0,0,3, 0,0,0,2,0,0,0, 0,0,0,2,2,0,0), 7,7)
x <- c(5,6); y <- 7
xx <- xx & !gac(mMMp, x,y, z=NULL, type)$gac
xx <- xx & !gac(mMMp, x,y, z= 1,   type)$gac
xx <- xx & !gac(mMMp, x,y, z= 4,   type)$gac
xx <- xx & !gac(mMMp, x,y, z= 2,   type)$gac
xx <- xx & !gac(mMMp, x,y, z= 3,   type)$gac
xx <- xx & !gac(mMMp, x,y, z=c(2,3), type)$gac
xx <- xx &  gac(mMMp, x,y, z=c(1,4), type)$gac
xx <- xx &  gac(mMMp, x,y, z=c(1,4,2), type)$gac
xx <- xx &  gac(mMMp, x,y, z=c(1,4,3), type)$gac
xx <- xx &  gac(mMMp, x,y, z=c(1,4,2,3), type)$gac

##################################################
## type = "pag" -- Tests from Ema
##################################################
type <- "pag"
pag.m <- readRDS(system.file("external/gac-pags.rds", package="pcalg"))
m1 <- pag.m[["m1"]]
x <- 6; y <- 9
xx <- xx & !gac(m1,x,y, z=NULL, type)$gac
xx <- xx & !gac(m1,x,y, z=1, type)$gac
xx <- xx & !gac(m1,x,y, z=2, type)$gac
xx <- xx & !gac(m1,x,y, z=3, type)$gac
xx <- xx & !gac(m1,x,y, z=4, type)$gac
xx <- xx & !gac(m1,x,y, z=c(2,3), type)$gac
xx <- xx &  gac(m1,x,y, z=c(2,3,8), type)$gac
xx <- xx &  gac(m1,x,y, z=c(2,3,7,8), type)$gac
xx <- xx & !gac(m1,x,y, z=c(2,3,5,8), type)$gac
xx <- xx & !gac(m1,x,y, z=c(2,3,5,7,8), type)$gac

x <- c(6,8); y <- 9
xx <- xx & !gac(m1,x,y, z=NULL, type)$gac
xx <- xx & !gac(m1,x,y, z=1, type)$gac
xx <- xx & !gac(m1,x,y, z=2, type)$gac
xx <- xx & !gac(m1,x,y, z=3, type)$gac
xx <- xx & !gac(m1,x,y, z=4, type)$gac
xx <- xx &  gac(m1,x,y, z=c(2,3), type)$gac
xx <- xx &  gac(m1,x,y, z=c(2,3,4), type)$gac
xx <- xx &  gac(m1,x,y, z=c(2,3,7), type)$gac
xx <- xx & !gac(m1,x,y, z=c(2,3,5), type)$gac
xx <- xx & !gac(m1,x,y, z=c(2,3,5,7), type)$gac

x <- 3; y <- 1
xx <- xx & !gac(m1,x,y, z=NULL, type)$gac
xx <- xx & !gac(m1,x,y, z=2, type)$gac
xx <- xx & !gac(m1,x,y, z=4, type)$gac
xx <- xx & !gac(m1,x,y, z=5, type)$gac
xx <- xx & !gac(m1,x,y, z=6, type)$gac
xx <- xx & !gac(m1,x,y, z=c(2,6), type)$gac
xx <- xx & !gac(m1,x,y, z=c(2,8), type)$gac
xx <- xx & !gac(m1,x,y, z=c(2,7,8), type)$gac
xx <- xx & !gac(m1,x,y, z=c(2,5,8), type)$gac
xx <- xx & !gac(m1,x,y, z=c(2,5,7,8), type)$gac

m2 <- pag.m[["m2"]]
x <- 3; y <-1
xx <- xx & !gac(m2,x,y, z=NULL, type)$gac
xx <- xx &  gac(m2,x,y, z=2, type)$gac
xx <- xx & !gac(m2,x,y, z=4, type)$gac
xx <- xx & !gac(m2,x,y, z=c(2,8), type)$gac
xx <- xx & !gac(m2,x,y, z=8, type)$gac
xx <- xx & !gac(m2,x,y, z=9, type)$gac
xx <- xx & !gac(m2,x,y, z=c(2,8,9), type)$gac
xx <- xx &  gac(m2,x,y, z=c(2,5), type)$gac

x <- c(3,9); y <- 1
xx <- xx & !gac(m2,x,y, z=NULL, type)$gac
xx <- xx & !gac(m2,x,y, z=2, type)$gac
xx <- xx & !gac(m2,x,y, z=4, type)$gac
xx <- xx & !gac(m2,x,y, z=c(2,8), type)$gac
xx <- xx & !gac(m2,x,y, z=8, type)$gac
xx <- xx & !gac(m2,x,y, z=9, type)$gac
xx <- xx & !gac(m2,x,y, z=c(2,8,9), type)$gac
xx <- xx & !gac(m2,x,y, z=c(2,5), type)$gac

m3 <- pag.m[["m3"]]
x <- 1; y <- 9
xx <- xx & !gac(m3,x,y, z=NULL, type)$gac
xx <- xx & !gac(m3,x,y, z=2, type)$gac
xx <- xx & !gac(m3,x,y, z=3, type)$gac
xx <- xx & !gac(m3,x,y, z=5, type)$gac
xx <- xx & !gac(m3,x,y, z=7, type)$gac
xx <- xx & !gac(m3,x,y, z=8, type)$gac
xx <- xx &  gac(m3,x,y, z=c(2,3), type)$gac
xx <- xx &  gac(m3,x,y, z=c(5,7), type)$gac

x <- 1; y <- 8
xx <- xx & !gac(m3,x,y, z=NULL, type)$gac
xx <- xx & !gac(m3,x,y, z=2, type)$gac
xx <- xx & !gac(m3,x,y, z=3, type)$gac
xx <- xx & !gac(m3,x,y, z=5, type)$gac
xx <- xx &  gac(m3,x,y, z=7, type)$gac
xx <- xx & !gac(m3,x,y, z=9, type)$gac
xx <- xx &  gac(m3,x,y, z=c(2,3), type)$gac
xx <- xx & !gac(m3,x,y, z=c(5,9), type)$gac

if (!xx) stop("Problem when testing function gac.")
