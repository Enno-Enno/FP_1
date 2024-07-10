all: build/V01.pdf

build/plot1.pdf: plot1.py data1.txt | build
	python plot1.py

build/V01.pdf: V01.tex *.tex lit.bib programme.bib | build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V01.tex
	BIBINPUTS=..:
	biber build/V01.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V01.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean