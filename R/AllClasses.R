##################################################
## Classes
##################################################
setClass("gAlgo",
         representation(call = "call",
                        n	   = "integer",
                        max.ord = "integer",
                        n.edgetests= "numeric",
                        sepset= "list",
                        pMax= "matrix"), "VIRTUAL")


setClass("fciAlgo",
         representation(amat = "matrix", allPdsep = "list",
                        n.edgetestsPDSEP = "numeric", max.ordPDSEP = "integer"),
         contains = "gAlgo")

setClass("pcAlgo",
         representation(graph = "graph", zMin = "matrix"), ## zMin for compatibility
         contains = "gAlgo")
