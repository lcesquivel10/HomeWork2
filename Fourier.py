import numpy as np
from numpy.fft import fft2, ifft2, fftshift, ifftshift
from scipy import misc,ndimage

#ixtransformada de fourier 2D
def propiaFFT2d(arr):
    transformada=[[] for i in range(len(arr))]
    for i1 in range(len(arr)):
        for j1 in range(len(arr[i])):
            Xnm=0
            for i2 in range(len(arr)):
                for j2 in range(len(arr[i1])):
                    Xnm+=arr[i2][j2]*np.exp(-2j*np.pi*(i1*i2/len(arr)+j1*j2/len(arr[i1]))) ## formula 2d
            transformada[i1].append(((Xnm.real**2+Xnm.imag**2)**0.5))
        transformada[i1]=np.array(transformada[i1])
    transformada=np.array(transformada)
    return transformada

#Matriz  filtro 

def CrearMatrizFiltroFFT(cerca,coeficienteClaridad):
    matriz=[[-1 for i in range(170)] for j in range(254)]
    for i in range(254):
        for j in range(170):
            #Se calcula el complemento si es la imagen 2
            if cerca:
                matriz[i][j]=1-np.exp(-(pow(i-127,2)+pow(j-85,2))/coeficienteClaridad)
            
            else:
                matriz[i][j]=np.exp(-(pow(i-127,2)+pow(j-85,2))/coeficienteClaridad)
                
    return matriz


def crearImagenHibrida(nombreImagen1,nombreIMagen2): # funcion para imagen... usando unicamente paquetes
    # pasar las imagenes a una matriz .... https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.ndimage.imread.html
    imagen1 = ndimage.imread(nombreImagen1, flatten=True)
    imagen2 = ndimage.imread(nombreIMagen2, flatten=True)
    
    
    #calcular la transformada de fourier de ambas imagenes
    #como es una matriz de 2 dimensiones se utiliza la funcion fft2...https://docs.scipy.org/doc/numpy-1.15.1/reference/generated/numpy.fft.fftshift.html
    #O la implementacion propia
    
    # misc.imsave... scipy.. guarda arreglo como imagen
    transformada1= fftshift(fft2(imagen1))
    transformada2=fftshift(fft2(imagen2))
    misc.imsave("transformadaImagen1.pdf", np.real(transformada1))
    misc.imsave("transformadaImagen2.pdf", np.real(transformada2))
    
    
    # FILTROS
    #Crear una matriz filtro la cual al multiplicarla por la fft nos permita cambial la claridad de la imagen
    filtro1=CrearMatrizFiltroFFT(True,1000)
    filtro2=CrearMatrizFiltroFFT(False,300)
    misc.imsave("filtro1.pdf", np.real(filtro1))
    misc.imsave("filtro2.pdf", np.real(filtro2))  
    
    #multiplicar estas matrices por las transformadas
    imagenFiltrada1=filtro1*transformada1
    imagenFiltrada2=filtro2*transformada2
    misc.imsave("ImagenFiltrada1.pdf", np.real(imagenFiltrada1))
    misc.imsave("ImagenFiltrada2.pdf", np.real(imagenFiltrada2))  

    #INVERSA para ver los cambios que hizo el filtro... https://www.programcreek.com/python/example/61036/numpy.fft.ifftshift
    imagenNueva1=ifft2(ifftshift(imagenFiltrada1))
    imagenNueva2=ifft2(ifftshift(imagenFiltrada2))
    misc.imsave("ImagenNueva1.pdf", np.real(imagenNueva1))
    misc.imsave("ImagenNueva2.pdf", np.real(imagenNueva2))
 
    #COMBINAR.. para tener una sola imagen
    imagenFinal=imagenNueva1+imagenNueva2
    misc.imsave("ImHybrid.pdf", np.real(imagenFinal))
    

crearImagenHibrida("cara_03_grisesMF.png","cara_02_grisesMF.png")
