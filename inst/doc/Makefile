whoami=$(shell whoami)

all: clean pdf
	
pdf:
	cp -f vignette-src/predictionet.Rnw .
	cp -f vignette-src/biblio.bib .
	cp -f vignette-src/predictionet-regrnet_design.pdf .
	cp -f vignette-src/predictionet-cytoscape.pdf .
	cp -f vignette-src/predictionet-pn_webapp_ras.pdf .
	"${R_HOME}/bin/R" --vanilla -e "library(pgfSweave); Stangle('predictionet.Rnw'); pgfSweave('predictionet.Rnw');"
	pdflatex predictionet
	bibtex predictionet
	pdflatex predictionet
	pdflatex predictionet
	if test -d cache; then rm -fr cache; fi
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.tikz *.toc *.map *.sh *.dvi *.bib
	rm -f predictionet-*.pdf Rplots.pdf
	cp -f predictionet.Rnw predictionet2.Rnw
	cat predictionet2.Rnw | perl  -ne 's/>>=$$/>>=\nif(dorun) {/;print $$_' > tmp && mv -f tmp predictionet2.Rnw
	cat predictionet2.Rnw | perl  -ne 's/\@$$/}\n\@/;print $$_' > tmp && mv -f tmp predictionet2.Rnw
	echo "<<run>>=\ndorun <- FALSE\n@\n\n" | cat - predictionet2.Rnw > predictionet.Rnw
	rm -f predictionet2.Rnw
	if test -d ../../build; then rm -fr ../../build; fi
	rm -f predictionet.R

clean:
	if test -d ../../build; then rm -fr ../../build; fi
	if test -d cache; then rm -fr cache; fi
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.tikz *.toc *.map *.sh *.dvi *.bib

