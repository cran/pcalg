(TeX-add-style-hook "pcalg_overview"
 (function
  (lambda ()
    (LaTeX-add-labels
     "fig:trueDAG"
     "fig:corGraph"
     "fig:skel"
     "fig:skelLWD")
    (TeX-run-style-hooks
     "subfigure"
     "texab"
     "graphicx"
     "latex2e"
     "art12"
     "article"
     "a4paper"
     "12pt"))))

