# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

"""
DATE: 04/07/17
OBJECTIF: Test Ellipsoïde contour lignes magnétiques 
"""

def polynomel(x):
    
    w = -6.0*np.abs(x)*np.abs(x)*np.abs(x)*np.abs(x)*np.abs(x) + 15.0*x*x*x*x - 10.0*np.abs(x)*np.abs(x)*np.abs(x) + 1
    
    y = ( w if (np.abs(x) <= 1.0)  else 0.0)
    
    return y
#    y = ( x if (np.abs(x) <= 1.0)  else 0.0)
#    w = -float(int(10*y))/10 + y
#    z = (int(w*100) if (0 != int(w*100)) else 0.0)
#
#    return z

#    y = ( x if (np.abs(x) <= 1.0)  else 0.0)
#
#    if y == 1.0 :
#        w = 1.0
#    elif y==0.5 :
#        w = 0.5
#    else :
#        w =0.0
#
#    return w

def coefpol4(x,y,a,b):
    a4 = 1
    a3 = 2*a +2*b 
    a2 = a**2 + 4*a*b + b**2 - x**2 - y**2
    a1 = 2*b*a**2 + 2*a*b**2 - 2*b*x**2 - 2*a*y**2
    a0 = a*a*b*b - x*x*b*b - y*y*a*a

    return a4, a3, a2, a1, a0

def equations(p,*kargs):
    t, de = p
    
    a,b,x,y = kargs
    
    return (x-(b+de)*np.cos(t) , 
            y-(a+de)*np.sin(t) )


def equations2(p,*kargs):
    al, be = p
    
    a,b,x,y = kargs
    
    return (-be * x/(b*b) + al*y/(a*a)  - al* be * (1/(a*a) - 1/(b*b)) , 
            al *al/(a*a) + be*be/(b*b) - 1 )
def equation(de,*kargs):
        
    a,b,x,y = kargs
    
    return x**2/((a+de)**2) + y**2/((b+de)**2) -1 

def ellipseB(x,y,x0,y0,a,b,R,B0):
# calcul la valeur des composantes du champ B par rapport à sa position
# x0,y0 : centre ellipse 
# a, b  : semi-major & semi minor
# R :     half width of the B ribbon 
    wx = x - x0
    wy = y - y0
    
    #modulus of arg
    wr = np.sqrt(wx**2 + wy**2)

    #angle of arg
    ar = (0 if wr == 0 else np.arccos(wy/wr))

    we = a*b/np.sqrt((b*np.cos(ar))**2 + (a*np.sin(ar))**2)
    if wx == 0.0 and wy ==0.0 :
    	print we
    wa = (wr - we)/R

    ux = (0 if wr == 0 else wy/wr)
    uy = (0 if wr == 0 else wx/wr)

    Bx = B0 * polynomel(wa)*ux
    By = B0 * polynomel(wa)*uy

    return Bx,By



def ellipseBis(x,y,x0,y0,a,b,R,B0):
# calcul la valeur des composantes du champ B par rapport à sa position
# x0,y0 : centre ellipse 
# a, b  : semi-major & semi minor
# R :     half width of the B ribbon 
    wx = x - x0
    wy = y - y0

    co = wx**2 / (a*a) +wy**2 /(b*b)
    wa = np.abs(1-np.sqrt(co)) * a/R

    #modulus of arg
    wr = np.sqrt(wx**2 + wy**2)

    ux = (0 if wr == 0 else wy/wr)
    uy = (0 if wr == 0 else wx/wr)

    Bx = B0 * polynomel(wa)*ux
    By = B0 * polynomel(wa)*uy

    return Bx,By


def ellipseTer(x,y,x0,y0,a,b,R,B0):
# calcul la valeur des composantes du champ B par rapport à sa position
# x0,y0 : centre ellipse 
# a, b  : semi-major & semi minor
# R :     half width of the B ribbon 
    wx = x - x0
    wy = y - y0
    
    #modulus of arg
    wr = np.sqrt(wx**2 + wy**2)
    
    #angle of arg
    ar = (0 if wr == 0 else np.arccos(wx/wr))
    
    args = (a,b,wx,wy)
    #delta =  fsolve(equation, np.sqrt(a**2+b**2)-wr,args=args)
    coefs = coefpol4(wx,wy,a,b)
    delta = np.roots(coefs)
    tmp_d = np.min(np.abs(delta))
    wa = tmp_d / R

    #modulus of vect
    wv = np.sqrt((wy/((b+tmp_d)**2))**2 + (wx/((a+tmp_d)**2))**2)

    ux = (0 if wr == 0 else wy/(wv*(b+tmp_d)**2))
    uy = (0 if wr == 0 else wx/(wv*(a+tmp_d)**2))

    Bx = B0 * polynomel(wa)*ux
    By = B0 * polynomel(wa)*uy

    return Bx,By

def ellipseTerDisplay(x,y,x0,y0,a,b,R,B0):
# calcul la valeur des composantes du champ B par rapport à sa position
# x0,y0 : centre ellipse 
# a, b  : semi-major & semi minor
# R :     half width of the B ribbon 
    wx = x - x0
    wy = y - y0
    
    #modulus of arg
    wr = np.sqrt(wx**2 + wy**2)
    
    #angle of arg
    ar = (0 if wr == 0 else np.arccos(wx/wr))
    
    args = (a,b,wx,wy)
    #delta =  fsolve(equation, np.sqrt(a**2+b**2)-wr,args=args)
    coefs = coefpol4(wx,wy,a,b)
    delta = np.roots(coefs)
    #print ' DELTA : ', delta
    tmp_d = np.min(np.abs(delta))
    #print ' RETENU : -> ',tmp_d 
    wa = tmp_d / R

    #modulus of vect
    wv = np.sqrt((wy/((b+tmp_d)**2))**2 + (wx/((a+tmp_d)**2))**2)

    ux = (0 if wr == 0 else wy/(wv*(b+tmp_d)**2))
    uy = (0 if wr == 0 else wx/(wv*(a+tmp_d)**2))

    Bx = B0 * polynomel(wa)*ux
    By = B0 * polynomel(wa)*uy

    return Bx,By,ux,uy,tmp_d

def fourierFlux(ncells,L,Bx,By):
      
    import pyfftw
    from pylab import meshgrid
    
    numOfThreads = 2
    dim = [ncells+1,ncells+1]
    
    # arrays have to be aligned... for pyfftw to work
    bx = pyfftw.empty_aligned(dim, dtype = 'complex128')
    by = pyfftw.empty_aligned(dim, dtype = 'complex128')
    BX = pyfftw.empty_aligned(dim, dtype = 'complex128')
    BY = pyfftw.empty_aligned(dim, dtype = 'complex128')
    AZ = pyfftw.empty_aligned(dim, dtype = 'complex128')
    az = pyfftw.empty_aligned(dim, dtype = 'complex128')
    
    # 2d fftw plans
    fftForx = pyfftw.FFTW(bx, BX, axes=(0, 1), direction='FFTW_FORWARD', flags=('FFTW_MEASURE', ))
    fftFory = pyfftw.FFTW(by, BY, axes=(0, 1), direction='FFTW_FORWARD', flags=('FFTW_MEASURE', ))
    fftBack = pyfftw.FFTW(AZ, az, axes=(0, 1), direction='FFTW_BACKWARD', flags=('FFTW_MEASURE', ))
    
    # fill the input arrays
    bx.real[...] = Bx
    by.real[...] = By
    bx.imag = np.zeros(dim, dtype = float)
    by.imag = np.zeros(dim, dtype = float)
    
    # update needed before execute plan
    fftForx.update_arrays(bx, BX)
    fftFory.update_arrays(by, BY)
    pyfftw.FFTW.execute(fftForx)
    pyfftw.FFTW.execute(fftFory)
    
    # this intend to build the k values, first half is ascending
    # positives values, last one being the descending negative values
    
    axis = 0
    N = dim[axis]
    k = np.linspace(0, N-1, N)
    kplus = k[:N/2+1]
    kminus =-k[(N-1)/2:0:-1]
    kx = np.concatenate((kplus, kminus))*2*np.pi/L
    axis = 1
    N = dim[axis]
    k = np.linspace(0, N-1, N)
    kplus = k[:N/2+1]
    kminus =-k[(N-1)/2:0:-1]
    ky = np.concatenate((kplus, kminus))*2*np.pi/L
    
    # create 2d values to avoid nested "for" loops
    wy, wx = meshgrid(ky, kx)
    k2 = np.square(wx)+np.square(wy)
    k2[0][0] = 1.0
    AZ = np.divide(1j*(wx*BY-wy*BX), k2)
    
    # remove dc component while k=0 mode is not defined
    AZ[0][0] = 0.0
    fftBack.update_arrays(AZ, az)
    pyfftw.FFTW.execute(fftBack)
    _flux = az.real
    
    return _flux

# -------------   Initialisation des paramètres  ----------#
B0 = 1.0 
R  = 2
a  = 4.0
b  = 6.0
x0 = 10.0
y0 = 10.0

L0 = 20.0  # length of the B box
X0 = 120   # # of cells

l  = np.linspace(0.0, L0, num = X0 + 1)
dl = l[1] - l[0]

vect = range(20,60)

Bx = np.zeros(len(vect))
By = np.zeros(len(vect))
ux = np.zeros(len(vect))
uy = np.zeros(len(vect))
d  = np.zeros(len(vect))
for index,i in enumerate(vect):
    [Bx[index],By[index],ux[index],uy[index],d[index]] = ellipseTerDisplay(i*dl,i*dl,x0,y0,a,b,R,B0)
        


Bx = -Bx
ux = -ux
#Az = fourierFlux(X0,L0,Bx,By)
