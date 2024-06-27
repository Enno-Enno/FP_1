all: build/V46.pdf

build/plot1.pdf: plot1.py data1.txt | build
	python plot1.py

build/V46.pdf: build/plot1.pdf V46.tex *.tex lit.bib programme.bib | build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V46.tex
	BIBINPUTS=..:
	biber build/V46.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V46.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean