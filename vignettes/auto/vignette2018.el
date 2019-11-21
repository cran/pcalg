(TeX-add-style-hook
 "vignette2018"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "amsmath"
    "amssymb"
    "natbib"
    "xcolor"
    "hyperref")
   (TeX-add-symbols
    '("pkg" 1))
   (LaTeX-add-labels
    "sec:intro-structure-learning"
    "sec:skeleton"
    "fig:intro1"
    "sec:fci"
    "fig:intro2"
    "sec:ges"
    "eqn:score-definition"
    "fig:gesFit"
    "sec:gies"
    "fig:giesFit"
    "sec:arges"
    "sec:addBK"
    "fig:addBgKnowledge"
    "sec:slAssumptions"
    "lab"
    "defadjustment"
    "eqn:totalEffect"
    "fig:ida"
    "fig:jointIda"
    "fig:backdoor"
    "fig:gac"
    "sec:assumptions"
    "sec:amat"
    "fig:rgraphviz"
    "fig:igraph"
    "fig:simFig4"
    "fig:simFig5")
   (LaTeX-add-bibliographies
    "Mybibliography")
   (LaTeX-add-xcolor-definecolors
    "Red"
    "Blue"))
 :latex)

