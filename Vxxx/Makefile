all: build/Vxxx.pdf

#build/plot1.pdf: plot1.py data1.txt | build
#	python plot1.py

build/Vxxx.pdf: Vxxx.tex | build
	lualatex --output-directory=build --interaction=batchmode --halt-on-error Vxxx.tex
	#biber build/Vxxx.bcf
	lualatex --output-directory=build --interaction=batchmode --halt-on-error Vxxx.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean