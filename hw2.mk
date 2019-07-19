prueba.tex: plots.png GenerarDatos

plots.png: data_euler.dat data_leap.dat data_runge.dat
	python Plots_hw2.py

%.dat : a.out
	./a.out

a.out: ODEs.cpp
	g++ -o a.out ODEs.cpp
       
GenerarDatos: fourier.py
	python3 fourier.py
    

