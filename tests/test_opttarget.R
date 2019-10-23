# Tests for optimal intervention targets.

library(pcalg)

# Test graph: essential graph of Figure 4 in [1]. It has an unoriented component
# with 3 vertices (1, 2, 3) and one with 2 vertices (4, 5).
#
# [1] A. Hauser, P. Bühlmann: Two optimal strategies for active learning of
# causal models from interventions. Proceedings of the 6th European Workshop on
# Probabilistic Graphical Models (PGM-2012), pp. 123 to 130, 2012.
essgraph <- new("EssGraph",
                nodes = as.character(1:5),
                in.edges = list(2, c(1, 3), 2, c(2, 3, 5), c(2, 3, 4)),
                targets = list(integer(0), 1:3))

# Optimal single vertex intervention: vertex 2 (==> all edges between nodes 1,
# 2, 3 become orientable).
stopifnot(opt.target(essgraph, max.size = 1, use.node.names = FALSE) == 2)
stopifnot(opt.target(essgraph, max.size = 1) == "2")

# Optimal intervention of arbitrary size: vertices 1, 3, 4 (makes all edges
# orientable).
stopifnot(all.equal(opt.target(essgraph, max.size = 5, use.node.names = FALSE),
                    c(1, 3, 4)))
stopifnot(all.equal(opt.target(essgraph), c("1", "3", "4")))
