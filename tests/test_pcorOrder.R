library(pcalg)

set.seed(123)

g <- randomDAG(4,0.6)
dat <- rmvDAG(1000,g,errDist="normal")

res1 <- round(pcorOrder(3,4,c(2,1),cor(dat)),8)
res2 <- round(pcorOrder(3,4,2,cor(dat)),8)

if((res1!=-0.00474175) | (res2!=0.00912709)) {
  stop("Test of pcorOrder: Theoretical result not matched!")
}
