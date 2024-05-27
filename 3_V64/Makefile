all: build/V64.pdf

build/plot1.pdf: plot1.py data1.txt | build
	python plot1.py

build/V64.pdf: V64.tex *.tex lit.bib programme.bib | build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V64.tex
	BIBINPUTS=..:
	biber build/V64.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V64.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean