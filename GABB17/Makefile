pdf:
	pdflatex GABB17.tex
	bibtex GABB17
	pdflatex GABB17.tex
	pdflatex GABB17.tex

GabbBC4M.o: GabbBC4M.c
	gcc -c --std=c11 GabbBC4M.c

clean:
	rm -f GABB17.aux GABB17.log GABB17.out GABB17.toc GABB17.bbl GABB17.blg GABB17.pdf GabbBC4M.o
