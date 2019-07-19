prueba.tex: plots.png GenerarDatos

plots.png: data_euler.dat data_leap.dat data_runge.dat
	python Plots_hw2L.py

%.dat : a.out
	./a.out

a.out: ODEsL.cpp
	g++ -o a.out ODEsL.cpp
       
GenerarDatos: fourier.py
	python3 fourier.py
    

