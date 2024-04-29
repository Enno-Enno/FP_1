all: build/V18.pdf

build/plot1.pdf: plot1.py data1.txt | build
	python plot1.py

build/V18.pdf: V18.tex | build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V18.tex
	BIBINPUTS=..:
	biber build/V18.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V18.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean