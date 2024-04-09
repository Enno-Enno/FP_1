#!/bin/bash
echo "Ordnername"
read dir
echo "Name Protokolldatei ohne Endung"
read Name
mkdir -p $dir
echo "all: build/$Name.pdf

# hier Python-Skripte:
build/plot.pdf: plot.py ../matplotlibrc ../header-matplotlib.tex | build 
        # so that matplotlib can find the tex header when running        
        # LaTeX in the tmp directory
        # and set the matplotlibrc
        TEXINPUTS=$$(pwd)/..: MATPLOTLIBRC=../matplotlibrc python plot.py

# hier weitere Abhängigkeiten für build/$Name.pdf deklarieren:
build/$Name.pdf: build/plot.pdf

build/$Name.pdf: FORCE | build
        # to find header and bib files in the main directory
        TEXINPUTS=..: \\
        BIBINPUTS=..: \\
        max_print_line=1048576 \\
        latexmk \\
          --lualatex \\
          --output-directory=build \\
          --interaction=nonstopmode \\
          --halt-on-error \\
        $Name.tex

build:
        mkdir -p build

clean:
        rm -rf build

FORCE:

.PHONY: all clean" > $dir/Makefile
echo "\\input{header.tex}      

\\subject{VERSUCH NUMMER}
\\title{TITEL}
\\date{%
  Durchführung: DATUM
  \\hspace{3em}
  Abgabe: DATUM
}

\\begin{document}

\\maketitle
\\thispagestyle{empty}
\\tableofcontents
\\newpage

\\input{content/theorie.tex}
\\input{content/durchfuehrung.tex}
\\input{content/auswertung.tex}
\\input{content/diskussion.tex}

\\printbibliography{}

\\end{document}" > $Name
echo "import matplotlib.pyplot as plt
import numpy as np" > plot.py