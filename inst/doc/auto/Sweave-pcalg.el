(TeX-add-style-hook "Sweave-pcalg"
 (function
  (lambda ()
    (LaTeX-add-labels
     "fig:trueDAG"
     "fig:corGraph"
     "fig:skel"
     "fig:skelLWD")
    (TeX-run-style-hooks
     "a4wide"
     "latex2e"
     "art10"
     "article"
     "a4paper"))))

