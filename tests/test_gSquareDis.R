library(pcalg)

load("discr100k.rda")# in directory tests/ i.e., typically *not* installed
dim(dat)# 100'000 x 5
suffStat <- list(dm = dat, nlev = c(3,2,3,4,2), adaptDF = TRUE)
.pt <- proc.time()

stopifnot(## Collider
	  disCItest(1,2,NULL,suffStat) > 0.05, ## Must be !=0
	  disCItest(1,2,3,   suffStat) < 0.05, ## 0
	  disCItest(5,2,4,   suffStat) < 0.05, ## 0
	  disCItest(5,2,NULL,suffStat) > 0.05, ## != 0
	  ## Non-collider
	  disCItest(4,3, 2,  suffStat) > 0.05, ## !=0
	  disCItest(4,3,NULL,suffStat) < 0.05, ## 0
	  ## Neighbours
	  disCItest(4,2,NULL,suffStat) < 0.05, ## 0
	  disCItest(1,3,NULL,suffStat) < 0.05, ## 0
	  disCItest(1,5,NULL,suffStat) > 0.05) ## !=0

cat('Time elapsed: ', proc.time() - .pt,"\n")

