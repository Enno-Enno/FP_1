all: build/v46.pdf

build/plot1.pdf: plot1.py data1.txt | build
	python plot1.py

build/v46.pdf: build/plot1.pdf v46.tex lit.bib programme.bib | build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error v46.tex
	BIBINPUTS=..:
	biber build/v46.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error v46.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean