LATEX       = pdflatex
CHECK_RERUN = grep "Rerun to get" $*.log

all: openclusters_mkness_A.pdf

%.pdf: %.tex %.bbl
	${LATEX} $<
	( ${CHECK_RERUN} && ${LATEX} $< ) || echo "Done."
	( ${CHECK_RERUN} && ${LATEX} $< ) || echo "Done."
