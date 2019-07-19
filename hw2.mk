
##GeneradPdf: GenerarDatos resultados.tex
##	pdflatex resultados.tex
##	find . -type f -not -name '*py' -not -name '*tex' -not -name '*mk' -not -name ##'*png' -not -name 'resultados.pdf' -print0 | xargs -0 rm --	


prueba.tex: plots.png GenerarDatos ## prueba engaño.. cabeza engaño

plots.png: data_euler.dat data_leap.dat data_runge.dat
	python Plots_hw2.py

%.dat : a.out
	./a.out

a.out: ODEs.cpp
	g++ -o a.out ODEs.cpp
       
GenerarDatos: Fourier.py
	python3 Fourier.py
    

