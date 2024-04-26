all: build/V60.pdf

build/plot1.pdf: plot1.py data1.txt | build
	python plot1.py

build/V60.pdf: V60.tex lit.bib| build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V60.tex
	BIBINPUTS=..:
	biber build/V60.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V60.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean