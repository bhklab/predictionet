whoami=$(shell whoami)

all: clean predictionet

predictionet:
	cp -f vignette-src/predictionet.Rnw .
	cp -f vignette-src/biblio.bib .
	cp -f vignette-src/priors_ras_from_webapp.csv .
	"${R_HOME}/bin/R" --vanilla -e "library(pgfSweave); pgfSweave('predictionet.Rnw');"
	pdflatex predictionet
	bibtex predictionet
	pdflatex predictionet
	pdflatex predictionet
	if test -d cache; then rm -fr cache; fi
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.tikz *.toc *.map *.sh *.dvi
	rm -f predictionet-*.pdf priors_ras_from_webapp.csv biblio.bib
	rm -f Rplots.pdf
	cp -f vignette-src/predictionet-fake.Rnw predictionet.Rnw

clean:
	if test -d cache; then rm -fr cache; fi
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.tikz *.toc *.map *.sh *.dvi
	rm -f predictionet-*.pdf priors_ras_from_webapp.csv biblio.bib
	rm -f Rplots.pdf
	cp -f vignette-src/predictionet-fake.Rnw predictionet.Rnw

pdfclean:
	rm -f predictionet.pdf