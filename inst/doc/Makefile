
## There is an old bug in texidvi that makes it not swallow the ~
## marker used to denote whitespace. This is actually due to fixing
## another bug they whereby you could not run texidvi on directory
## names containing a tilde (as we happen to do for Debian builds of R
## alpha/beta/rc releases). The 'tilde' bug will go away as it
## reportedly has been squashed upstream but I am still bitten by it
## on Ubuntu so for Dirk will insist on pdflatex and this variable helps
whoami=$(shell whoami)

all: clean predictionet.pdf

pdfclean:
	rm -f predictionet.pdf

clean:
	if test -d cache; then rm -fr cache; fi
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.tikz *.toc *.map *.sh *.dvi
	rm -f predictionet-*.pdf priors_ras_from_webapp.csv biblio.bib
	cp -f vignette-src/predictionet-fake.Rnw predictionet.Rnw

predictionet.pdf: vignette-src/predictionet.Rnw
	cp -f vignette-src/predictionet.Rnw .
	cp -f vignette-src/biblio.bib .
	cp -f vignette-src/priors_ras_from_webapp.csv .
	Rscript --vanilla -e "library(pgfSweave); pgfSweave('predictionet.Rnw');"
	pdflatex predictionet
	bibtex predictionet
	pdflatex predictionet
	pdflatex predictionet
	if test -d cache; then rm -fr cache; fi
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.tikz *.toc *.map *.sh *.dvi
	rm -f predictionet-*.pdf priors_ras_from_webapp.csv biblio.bib
	cp -f vignette-src/predictionet-fake.Rnw predictionet.Rnw
