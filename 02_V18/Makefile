all: build/V18.pdf

build/plt1_Eu.pdf: plot1.py 152Eu-Spektrum.txt | build
	python plot1.py

build/plt5_Cs.pdf: plot2.py 137Cs-Spektrum.txt | build
	python plot2.py

build/plt9_Un.pdf: plot3.py 000Unbekannt-Spektrum.txt | build
	python plot3.py

build/V18.pdf: build/plt1_Eu.pdf build/plt5_Cs.pdf build/plt9_Un.pdf V18.tex | build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V18.tex
	BIBINPUTS=..:
	#biber build/V18.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error V18.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean