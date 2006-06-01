library(pcalg)

n <- 1000
p <- 10

set.seed(101)
dag <- randomDAG(p, 0.2)

checkNode <- sample(p, 1)
##Bug conMat <- as(dag,"matrix")
##Bug parents <- conMat[,checkNode]
parents <- wgtMatrix(dag)[checkNode, ]

data <- rmvDAG(n, dag, errDist="normal", mix=0.5)
realVal <- checkVal <- numeric(n)
for (i in 1:n) {
  realVal[i] <- data[i,checkNode]
  checkVal[i] <- sum(data[i,]*parents)
}
sampleVal <- realVal - checkVal
qqnorm(sampleVal)
shapiro.test(sampleVal)


