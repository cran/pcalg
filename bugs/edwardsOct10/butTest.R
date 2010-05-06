library(pcalg)
library(Rgraphviz)
load("kalisch.Rdata")

cpdag.carc <- pc(suffStat,gaussCItest,p=ncol(carcass),alpha=0.05, verbose = TRUE)
## nodes(cpdag.carc@graph)<-names(carcass)
## plot(cpdag.carc@graph,"neato")
plot(cpdag.carc)

write.table(carcass, file = "kalisch.dat", row.names = FALSE)

### Test
skel.carc <- skeleton(suffStat,gaussCItest,p=ncol(carcass),alpha=0.05, verbose = TRUE)
cpdag.carc <- udag2pdagSpecial(skel.carc, verbose = 1)
plot(cpdag.carc)
