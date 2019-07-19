prueba.tex: plots.png GenerarDatos ## prueba.tex engaño

plots.png: data_euler.dat data_leap.dat data_runge.dat
	python Plots_hw2.py

%.dat : a.out
	./a.out

a.out: ODEs.cpp
	g++ -o a.out ODEs.cpp
       
GenerarDatos: Fourier.py
	python3 Fourier.py
    

