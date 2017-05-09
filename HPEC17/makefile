#
#  USAGE:
#     make       ... to build the GrB_futures_hpec17.pdf file
#     make tar   ... to create a tar file of the paper's source
#     make clean ... to remove LaTeX working files
#

LATEX = pdflatex
BIBTEX = bibtex

BIBSRC = GrB_futures_hpec17.tex 

SRC  = $(BIBSRC) GrBDescrip.tex ack.tex conclusion.tex dag.tex \
       inspecExec.tex intro.tex kronProd.tex makefile masks.tex \
       stddefs.tex usrTypes.tex 

TARSRC = GrB.bib GrBDescrip.tex GrB_futures_hpec17.tex ack.tex \
         conclusion.tex dag.tex figures inspecExec.tex intro.tex \
         kronProd.tex makefile masks.tex stddefs.tex usrTypes.tex 

GrB_futures_hpec17.pdf: $(SRC)
	$(LATEX) $(SRC)
	$(BIBTEX) $(BIBSRC:.tex=)
	$(LATEX) $(SRC)
	$(LATEX) $(SRC)

tar: $(TARSRC)
	rm -rf GrB_futures_hpec17.tar
	tar cvf GrB_futures_hpec17.tar $(TARSRC) 

clean:
	rm -f *.dvi *.log *.aux *.bbl *.blg

cleanall: clean
	rm -f GrB_futures_hpec17.pdf
