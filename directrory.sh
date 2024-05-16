#!/bin/bash
echo "Ordnername"
read dir
echo "Name Protokolldatei ohne Endung"
read Name
mkdir -p $dir
printf "all: build/$Name.pdf

build/plot1.pdf: plot1.py data1.txt | build
	python plot1.py

build/$Name.pdf: $Name.tex lit.bib programme.bib | build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error $Name.tex
	BIBINPUTS=..:
	biber build/$Name.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error $Name.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean" > $dir/Makefile
cd $dir
cp -r ../Vxxx/Bilder .
cp ../Vxxx/plot1.py .
touch data1.txt
cp ../Vxxx/Vxxx.tex .
cp ../Vxxx/lit.bib .
cp ../Vxxx/programme.bib .
mv Vxxx.tex $Name.tex
make
cd ..