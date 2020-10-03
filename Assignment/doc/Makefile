all: instructions.pdf

%.pdf: %.tex
	pdflatex $*
	pdflatex $*

clean:
	rm *.pdf *.aux *.log *.out

