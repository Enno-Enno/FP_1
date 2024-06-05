all: build/V64.pdf

build/plot1.pdf: plot1.py data.py | build
	python plot1.py

build/plot2.pdf: plot2.py data.py | build
	python plot2.py

build/V64.pdf: V64.tex *.tex lit.bib programme.bib build/plot1.pdf build/plot2.pdf | build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V64.tex
	BIBINPUTS=..:
	biber build/V64.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V64.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean