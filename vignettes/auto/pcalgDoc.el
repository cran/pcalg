(TeX-add-style-hook
 "pcalgDoc"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("jss" "nojss" "shortnames" "article")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8")))
   (add-to-list 'LaTeX-verbatim-environments-local "alltt")
   (TeX-run-style-hooks
    "latex2e"
    "jss"
    "jss10"
    "graphicx"
    "xcolor"
    "framed"
    "alltt"
    "inputenc"
    "algorithmic"
    "algorithm"
    "amsmath"
    "tikz"
    "xstring"
    "upquote")
   (TeX-add-symbols
    '("gredge" 1)
    '("AUT" 1)
    '("hlkwd" 1)
    '("hlkwc" 1)
    '("hlkwb" 1)
    '("hlkwa" 1)
    '("hlstd" 1)
    '("hlopt" 1)
    '("hlcom" 1)
    '("hlstr" 1)
    '("hlnum" 1)
    "maxwidth"
    "hlipl"
    "FrameCommand"
    "substarrow")
   (LaTeX-add-labels
    "sec:introduction"
    "fig:intro1"
    "sec:gm"
    "pc"
    "algo:pc"
    "sec:bounds"
    "sec:skel"
    "fig:skelExpl"
    "fig:skel2"
    "sec:pc"
    "fig:pcFit1"
    "sec:ges"
    "fig:gesFit"
    "sec:fci"
    "fig:fci"
    "sec:rfci"
    "sec:gies"
    "fig:giesFit"
    "sec:ida"
    "fig:ida"
    "fig:allDags"
    "sec:backdoor"
    "fig:backdoor"
    "fig:userSpec"
    "fig:yeast"
    "fig:yeast.new")
   (LaTeX-add-environments
    "kframe"
    "knitrout")
   (LaTeX-add-bibliographies
    "Mybib")
   (LaTeX-add-lengths
    "edgelength")
   (LaTeX-add-color-definecolors
    "fgcolor"
    "shadecolor"
    "messagecolor"
    "warningcolor"
    "errorcolor"))
 :latex)

