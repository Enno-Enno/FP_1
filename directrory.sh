#!/bin/bash
echo "Ordnername"
read dir
echo "Name Protokolldatei ohne Endung"
read Name
mkdir -p $dir
printf "all: build/$Name.pdf

# hier Python-Skripte:
build/plot.pdf: plot.py ../matplotlibrc ../header-matplotlib.tex | build 
        # so that matplotlib can find the tex header when running        
        # LaTeX in the tmp directory
        # and set the matplotlibrc
        TEXINPUTS=\$\$(pwd)/..: MATPLOTLIBRC=../matplotlibrc python plot.py

# hier weitere Abhängigkeiten für build/$Name.pdf deklarieren:
build/$Name.pdf: build/plot.pdf

build/$Name.pdf: FORCE | build
\t\t# to find header and bib files in the main directory
\t\tTEXINPUTS=..: \\
\t\tBIBINPUTS=..: \\
\t\tmax_print_line=1048576 \\
\t\tlatexmk \\
\t\t  --lualatex \\
\t\t  --output-directory=build \\
\t\t  --interaction=nonstopmode \\
\t\t  --halt-on-error \\
\t\t$Name.tex

build:
\t\tmkdir -p build

clean:
\t\trm -rf build

FORCE:

.PHONY: all clean" > $dir/Makefile
echo "\\input{header.tex}      

\subject{VERSUCH NUMMER}
\\title{TITEL}
\\date{%
  Durchführung: DATUM
  \hspace{3em}
  Abgabe: DATUM
}

\\begin{document}

\\maketitle
\\thispagestyle{empty}
\\tableofcontents
\\newpage

\\section{Theorie}

\\printbibliography{}

\\end{document}" > $dir/$Name.tex
echo "import matplotlib.pyplot as plt
import numpy as np" > $dir/plot.py