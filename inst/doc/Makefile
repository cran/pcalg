###--- Makefile mainly from Douglas Bates

## need to get 'texclean' command

DRIVERS = pcalgDoc
#R_LIB = /usr/local/lib/R
R_LIB = /usr/local/app/R/current
#R_PROG = R-devel
R_PROG = R
#R_CALL = ${R_PROG} --vanilla --silent
#                  ^^^^^^^^^ also loses R_LIBS
R_CALL = ${R_PROG} --no-init-file --no-save --silent
# unfortunately, R CMD Sweave is no longer useful (as it has --vanilla in R 2.12.0)
# R_SWEAVE = ${R_PROG} CMD Sweave

TPUTS = ${TEXINPUTS}:${R_LIB}/share/texmf//

## "Writing R Extensions" requires: 1st target must also clean up
doc: all
	${MAKE} clean
## but this is a pain --->   use  'make all'  which does *not* delete!
all: $(DRIVERS:=.tex) $(DRIVERS:=.pdf)

%.tex: %.Snw
	echo "utils::Sweave(\"$<\")" | ${R_PROG} --no-restore --no-save --slave
#	${R_SWEAVE} $<

%.R: %.Rnw
	${R_PROG} CMD Stangle $<

%.pdf : %.tex *.bib
	TEXINPUTS=${TPUTS} texi2pdf $<

%.dvi : %.tex *.bib
	TEXINPUTS=${TPUTS} texi2dvi $<


purge:
	$(RM) Rplots.pdf $(DRIVERS:=.out) *.aux *.log *.blg
clean: purge
	$(RM) $(DRIVERS:=.tex) $(DRIVERS:=.bbl)
	$(RM) $(DRIVERS:=-*.pdf) $(DRIVERS:=-*.eps)
	$(RM) $(DRIVERS:=-*.pdf) $(DRIVERS:=-*.eps)

realclean: clean
	$(RM) $(DRIVERS:=.pdf)
