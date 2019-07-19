import numpy as np
from numpy.fft import fft2, ifft2, fftshift, ifftshift
from scipy import misc,ndimage
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm

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

## para graficar las imagenes 

def graficarImagenes(img1,img2,nombreGuardado,titulo):

    f = plt.figure()
    f.add_subplot(1,2, 1)
    plt.imshow(np.real(img1), cmap=cm.gray)
    f.add_subplot(1,2, 2)
    plt.imshow(np.real(img2), cmap=cm.gray)
    f.suptitle(titulo, fontsize=13)
    f.savefig(nombreGuardado)
    
def crearImagenHibrida(nombreImagen1,nombreIMagen2): # funcion para imagen... usando unicamente paquetes
    # pasar las imagenes a una matriz .... https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.ndimage.imread.html
    imagen1 = ndimage.imread(nombreImagen1, flatten=True)
    imagen2 = ndimage.imread(nombreIMagen2, flatten=False)
    #graficarImagenes(imagen1,imagen2,"imagenInicial.pdf","Imagenes Iniciales")
    
    #calcular la transformada de fourier de ambas imagenes
    #como es una matriz de 2 dimensiones se utiliza la funcion fft2...https://docs.scipy.org/doc/numpy-1.15.1/reference/generated/numpy.fft.fftshift.html
    #O la implementacion propia
    
    # misc.imsave... scipy.. guarda arreglo como imagen
    transformada1= fftshift(fft2(imagen1))
    transformada2=fftshift(fft2(imagen2))
    graficarImagenes(transformada1,transformada2,"FFtIm.pdf","Transformadas de las imagenes")
    
    
    # FILTROS
    #Crear una matriz filtro la cual al multiplicarla por la fft nos permita cambial la claridad de la imagen
    filtro1=CrearMatrizFiltroFFT(True,1000)
    filtro2=CrearMatrizFiltroFFT(False,300)
    graficarImagenes(filtro1,filtro2,"Improceso.pdf","Filtros") 
    
    #multiplicar estas matrices por las transformadas
    imagenFiltrada1=filtro1*transformada1
    imagenFiltrada2=filtro2*transformada2
    graficarImagenes(imagenFiltrada1,imagenFiltrada2,"Improceso2.pdf","Imagenes con el filtro")
    #INVERSA para ver los cambios que hizo el filtro... https://www.programcreek.com/python/example/61036/numpy.fft.ifftshift
    imagenNueva1=ifft2(ifftshift(imagenFiltrada1))
    imagenNueva2=ifft2(ifftshift(imagenFiltrada2))
    graficarImagenes(imagenNueva1,imagenNueva2,"Improceso3.pdf","Nuevas imagenes")
    #COMBINAR.. para tener una sola imagen
    imagenFinal=imagenNueva1+imagenNueva2
    misc.imsave("ImHybrid.pdf", np.real(imagenFinal))
      
## Imagen final
crearImagenHibrida("cara_03_grisesMF.png","cara_02_grisesMF.png")
