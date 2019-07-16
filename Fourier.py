import numpy as np
import matplotlib.pyplot as plt
# almacenar los datos 
imagen1 =plt.imread("cara_03_grisesMF.png")
imagen2 =plt.imread("cara_02_grisesMF.png")

def Trasformada (Fev,N): # entran puntos y numero de datos 
    TransFourier = [] # lista vacia para que añada los valores evaluados 
    N= len(Fev) # tamaño de mis datos un valor 
    for i in range (N): # revisar los picos en el espacio de fourier osea los n
        contadorFourier= 0 # contador para que 
        for k in range ( N ): # revisa las posiciones k 
            contadorFourier+= Fev[k]*np.exp(1j*2*np.pi*(i/N)*k)# sumatoria de fourier
            TransFourier.append(contadorFourier)
    return contadorFourier