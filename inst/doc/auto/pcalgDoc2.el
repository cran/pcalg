(TeX-add-style-hook "pcalgDoc2"
 (lambda ()
    (LaTeX-add-bibliographies
     "Mybib")
    (LaTeX-add-labels
     "fig:intro1"
     "sec:gm"
     "pc"
     "algo:pc"
     "sec:skel"
     "fig:skelExpl"
     "fig:skel2"
     "sec:pc"
     "fig:pcFit1"
     "sec:fci"
     "fig:fci"
     "sec:ida"
     "fig:ida"
     "fig:allDags"
     "fig:userSpec")
    (TeX-run-style-hooks
     "algorithm"
     "algorithmic"
     "latex2e"
     "jss10"
     "jss"
     "nojss"
     "article")))

