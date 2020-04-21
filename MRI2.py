from random import seed, random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import fft 

import matplotlib.pyplot as plt
import cv2
import pylab 
import time

dT = 1  
T = 500
N = int(np.ceil(T/dT)+1)
T1 = 600    
T2 = 100
color = ['red','orange','blue']
xdata = [[],[],[]]
ydata = [[],[],[]]
B=[]
dataxPlot = [[],[],[]]
datayPlot=[[],[],[]]
z = np.linspace(0, 5, N)
img = cv2.imread('2d-mri.png',0)
larmor_frequency= None
B =[]
zlam=np.arange(0,2000,1)

def LarmerFreq():
    seed(1)
    zlam=np.arange(0,2000,1)
    for _ in range(0,2000):
        B.append(random()*0.5+3)
    larmor_frequency = 42.6*np.array(B)*2*np.pi
    
    return(larmor_frequency)
  
def zrot(phi):   
    Rz = [[np.cos(phi) ,-np.sin(phi), 0],[np.sin(phi) ,np.cos(phi), 0],[ 0, 0 ,1]]
    return(Rz)

def freeprecess(dT ,T1 ,T2 , df):
    phi = df
    E1 = np.exp(-dT/T1) 
    E2 = np.exp(-dT/T2) 
    Afp = np.dot([[E2,0,0],[0,E2,0],[0,0,E1]],zrot(phi))      
    Bfp = [0 ,0 ,1-E1]   
    return(Afp,Bfp)

def Magnetization(A,B):
    M = np.empty((N,3))    
    M[0,:] =np.array([1,0,0])
    for k in range (N-1):
        M[k+1,:] = np.dot(A,M[k,:])+B
    return(M[:,0],M[:,1])


def blochEquation():  
    df = LarmerFreq()
    for i in range(3):
        A,B = freeprecess(dT,T1,T2,df[i+2])
        xdata[i],ydata[i]=Magnetization(A,B)
    plot(xdata,ydata)
  

def plot(xdata,ydata):
    ax = pylab.gca(projection='3d')   
    ax.set_xlim(min(xdata[0]), max(xdata[0]))
    ax.set_ylim(min(ydata[0]),max(ydata[0]))
    ax.set_zlim(0, 5)
    for i in range(1,50):
        for j in range(3):
            dataxPlot[j].append(xdata[j][i-1])
            datayPlot[j].append(ydata[j][i-1]) 
            pylab.plot(dataxPlot[j],datayPlot[j],z[:i],color =color[j],linewidth=1.5)

        pylab.draw()
        pylab.pause(1e-117)        
     
    plt.show() 



def transform_image_to_kspace(img, dim=None, k_shape=None):
  
    if not dim:
        dim = range(img.ndim)
    k = fft.fftshift(fft.fftn(fft.ifftshift(img, axes=dim), s=None, axes=dim), axes=dim)
    k /= np.sqrt(np.prod(np.take(img.shape, dim)))
    k = np.real(k)
    k = 20*np.log(np.abs(k)+1)
    return(k)

def KSpace():
    magnitudeSpectrum = transform_image_to_kspace(img)
    plt.subplot(121),plt.imshow(img,cmap='gray')
    plt.title('Original ')
    plt.subplot(122),plt.imshow(magnitudeSpectrum, cmap = 'gray',interpolation='nearest')
    plt.title('K sapce')
    plt.show()

def LarmerFreqPlot():
    df = LarmerFreq()
    plt.subplot(211)
    plt.suptitle('B')
    plt.plot(zlam,B)
    plt.xlabel('Z')
    plt.ylabel('B')
    axes = pylab.gca()
    axes.set_ylim(0,10)
    plt.subplot(212)
    plt.plot(zlam,df)
    plt.xlabel('z')
    plt.ylabel('Larmer Freq')
    plt.show()
    
# blochEquation()
KSpace()
# LarmerFreqPlot()
