all: build/V60_Benedikt_Lütke_Lanfer_Enno_Wellmann.pdf

# hier Python-Skripte:
build/plot.pdf: plot.py ../matplotlibrc ../header-matplotlib.tex | build 
		# so that matplotlib can find the tex header when running        
		# LaTeX in the tmp directory
		# and set the matplotlibrc
		TEXINPUTS=37055(pwd)/..: MATPLOTLIBRC=../matplotlibrc python plot.py

# hier weitere Abhängigkeiten für build/V60_Benedikt_Lütke_Lanfer_Enno_Wellmann.pdf deklarieren:
build/V60_Benedikt_Lütke_Lanfer_Enno_Wellmann.pdf: build/plot.pdf

build/V60_Benedikt_Lütke_Lanfer_Enno_Wellmann.pdf: FORCE | build
		# to find header and bib files in the main directory
		TEXINPUTS=..: \
		BIBINPUTS=..: \
		max_print_line=1048576 \
		latexmk \
		  --lualatex \
		  --output-directory=build \
		  --interaction=nonstopmode \
		  --halt-on-error \
		V60_Benedikt_Lütke_Lanfer_Enno_Wellmann.tex

build:
		mkdir -p build

clean:
		rm -rf build

FORCE:

.PHONY: all clean
