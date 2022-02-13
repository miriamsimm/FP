all: build/main.pdf

figures/T_1.pdf: T_1.py matplotlibrc header-matplotlib.tex 
	python T_1.py

figures/T_2.pdf: T_2.py matplotlibrc header-matplotlib.tex 
	python T_2.py

figures/echo.pdf figures/fourier_plot.pdf figures/fourier_trafo.pdf: diffusion.py matplotlibrc header-matplotlib.tex 
	python diffusion.py

build/main.pdf: figures/T_1.pdf figures/T_2.pdf figures/echo.pdf figures/fourier_plot.pdf figures/fourier_trafo.pdf

build/main.pdf: FORCE | build
	  TEXINPUTS=build: \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	main.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
