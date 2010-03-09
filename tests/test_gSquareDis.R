library(pcalg)

load(file = "discreteData.rda")
suffStat <- list(dm = dat, nlev = c(3,2,3,4,2), adaptDF = TRUE)
## Collider
ok <- all(disCItest(1,2,NULL,suffStat) > 0.05, ## Must be !=0
    disCItest(1,2,3,suffStat) < 0.05, ## 0
    disCItest(5,2,4,suffStat) < 0.05, ## 0
    disCItest(5,2,NULL,suffStat) > 0.05, ## != 0
    ## Non-collider
    disCItest(4,3,2,suffStat) > 0.05, ## !=0
    disCItest(4,3,NULL,suffStat) < 0.05, ## 0
    ## Neighbours
    disCItest(4,2,NULL,suffStat) < 0.05, ## 0
    disCItest(1,3,NULL,suffStat) < 0.05,## 0
    disCItest(1,5,NULL,suffStat) >0.05) ## !=0

if (!ok) stop("Test gSquareDis wrong: Some dependence was not estimated correctly!")
