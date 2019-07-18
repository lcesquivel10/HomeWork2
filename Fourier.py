import numpy as np
from numpy.fft import fft2, ifft2, fftshift, ifftshift
from scipy import misc,ndimage
import matplotlib.pyplot as plt


# almacenar los datos 
imagen1 =plt.imread("cara_03_grisesMF.png")
imagen2 =plt.imread("cara_02_grisesMF.png")


import numpy as np
from numpy.fft import fft2, ifft2, fftshift, ifftshift
from scipy import misc,ndimage

#implementacion propia transformada de fourier 2D

def propiaFFT2d(arr):
    transformada=[[] for i in range(len(arr))]
    for i1 in range(len(arr)):
        for j1 in range(len(arr[i])):
            Xnm=0
            for i2 in range(len(arr)):
                for j2 in range(len(arr[i1])):
                    Xnm+=arr[i2][j2]*np.exp(-2j*np.pi*(i1*i2/len(arr)+j1*j2/len(arr[i1])))
            transformada[i1].append(((Xnm.real**2+Xnm.imag**2)**0.5))
        transformada[i1]=np.array(transformada[i1])
    transformada=np.array(transformada)
    return transformada