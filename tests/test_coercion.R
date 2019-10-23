# Tests for coercion functions.

library(pcalg)

# Tests for coercion functions of EssGraph objects.
set.seed(123)

for (i in 1:10) {
  # Test that a CPDAG can be coerced to an EssGraph object and back.
  dag <- randomDAG(10, 0.3)
  cpdag <- dag2cpdag(dag)
  essgraph <- as(cpdag, "EssGraph")
  cpdag2 <- as(essgraph, "graphNEL")

  if (!all.equal(cpdag, cpdag2)) {
    stop("Original and coerced CPDAG are not equal.")
  }
}
