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
#Matriz la cual nos dará el filtro que deseamos
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


def crearImagenHibrida(nombreImagen1,nombreIMagen2):
    #primer paso pasar las imagenes a una matriz
    imagen1 = ndimage.imread(nombreImagen1, flatten=True)
    imagen2 = ndimage.imread(nombreIMagen2, flatten=True)
    misc.imsave("sonriendo.pdf", np.real(imagen1))
    misc.imsave("normal.pdf", np.real(imagen2))
    #segundo paso
    #calcular la transformada de fourier de ambas imagenes
    #como es una matriz de 2 dimensiones se utiliza la funcion fft2
    #O la implementacion propia
    transformada1= fftshift(fft2(imagen1))
    transformada2=fftshift(fft2(imagen2))
    misc.imsave("transformadaImagen1.pdf", np.real(transformada1))
    misc.imsave("transformadaImagen2.pdf", np.real(transformada2))
    #tercer paso
    #Crear una matriz filtro la cual al multiplicarla por la fft nos permita cambial la claridad de la imagen
    filtro1=CrearMatrizFiltroFFT(True,1000)
    filtro2=CrearMatrizFiltroFFT(False,300)
    misc.imsave("filtro1.pdf", np.real(filtro1))
    misc.imsave("filtro2.pdf", np.real(filtro2))  
    #cuarto paso
    #multiplicar estas matrices por las transformadas
    imagenFiltrada1=filtro1*transformada1
    imagenFiltrada2=filtro2*transformada2
    misc.imsave("ImagenFiltrada1.pdf", np.real(imagenFiltrada1))
    misc.imsave("ImagenFiltrada2.pdf", np.real(imagenFiltrada2))  
    #quinto paso
    #hacer la fft inversa a cada imagen para ver los cambios que hizo el filtro
    imagenNueva1=ifft2(ifftshift(imagenFiltrada1))
    imagenNueva2=ifft2(ifftshift(imagenFiltrada2))
    misc.imsave("ImagenNueva1.pdf", np.real(imagenNueva1))
    misc.imsave("ImagenNueva2.pdf", np.real(imagenNueva2))
    #Sexto paso
    #Ya teniendo estas dos imagenes con sus filtros podemos pasar a combinarlas a una sola imagen
    imagenFinal=imagenNueva1+imagenNueva2
    misc.imsave("ImagenFinal.pdf", np.real(imagenFinal))
    

crearImagenHibrida("cara_03_grisesMF.png","cara_02_grisesMF.png")
