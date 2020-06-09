# -*- coding: utf-8 -*-
"""
Comment on the script
...
Get several image and display it on a grid. Since they have the same dimension 
The axis are shared in order to simplify the readiness.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
from matplotlib import rc
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import ImageGrid
sys.path.append("/home/sbolanos/Documents/code/Python/")
import lectILZ as ilz
sys.path.append("/home/sbolanos/Documents/code/mypywi/")
import runs.heckle as heckle 
#.. use latex fonts
rc('text',usetex = True)
rc('font',size = 16,family ='sans-serif')
rc('axes', labelsize = 'larger')
rc('mathtext',default = 'regular')

def GetIndex(vect,val):
	#Donne l'indice de la valeur sur le vecteur ou l'indice inf si non exact
	return np.where(vect<=val)[0][-1]


# constant :
b0   = 100 # Tesla
n0   = 1.0e27 #m-3
A    = 63
Z    = 19
mu   = 4*np.pi*1.0e-7
kb   = 1.38e-23
mp   = 1.67e-27
qe   = 1.6e-19
mi   = A*mp

va = b0/np.sqrt(mu*n0*mi/Z)       # alfven velocity
dp = np.sqrt(mi/(mu*n0*Z*qe*qe))*1e3  # inertial length [mm]
tc = mi/(Z*qe*b0)*1e9                 # gyroperiod [ns]

G = 10.0 #Grossissement

paths    =  [['/media/sbolanos/BatDRIVE/HECKLE/PRO/anti60_T1/400x600_pro60_T1/',
			  '/media/sbolanos/BatDRIVE/HECKLE/PRO/T1/400x600_T1/',
			  '/media/sbolanos/BatDRIVE/HECKLE/PRO/pro60_T1/400x600_pro-60_T1/'],
			['/media/sbolanos/BatDRIVE/HECKLE/PRO/anti60_T10/400x600_pro60_T10/',
			 '/media/sbolanos/BatDRIVE/HECKLE/PRO/T10/400x600_T10/',
			 '/media/sbolanos/BatDRIVE/HECKLE/PRO/pro60_T10/400x600_pro-60_T10/']]
			 

time     = [[13.4, 9.6, 7.8],[ 19.4, 7.0, 4.8]]
#time     = [[30., 20., 18.],[ 23.0,15.0,11.]]
titles   = ['$anti-Hall$','$coplanaire$','$pro-Hall$']
ylabel   = ['$\\beta\ =\ 1$\n$y c / \omega_{pi}$','$\\beta\ =\ 10$\n$y c / \omega_{pi}$','$\\beta\ =\ 100$\n$y c / \omega_{pi}$']
subnot   = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)']
F = plt.figure(1,(11,6))
grid = ImageGrid(F, 111,
                      nrows_ncols= np.shape(paths),
                      direction="row",
                      axes_pad=0.2,
                      add_all=True,
                      label_mode="L",
                      share_all=False,
                      cbar_location="right",
                      cbar_mode="single",
                      cbar_size="2%",
                      cbar_pad=0.05,
                      )
n_beta = np.shape(paths)[0]
n_conf = np.shape(paths)[1]

index = 0
name = '400x600'
vmaxh = 0.4
vminh = - vmaxh 

x_max =  30.
x_min = -30.
y_max =  20.
y_min = -20.




for i in range(n_beta):
	
	for j in range(n_conf):

		if i == 0 :
			grid[index].set_title(titles[index])
		
		path = paths[i][j]
		t    = time[i][j]

		run  = heckle.Heckle(path, name)

		dataAz = run.fourierFlux(t)
		#dataEz = run.GetE(t)[...,2]
		dataBz = run.GetB(t)[...,2]

		if index == 0 : 

			X = run.GetCoords(axis=0)
			X = X - X[-1]/2
			Y = run.GetCoords(axis=1)
			Y = Y - Y[-1]/2

			ix_max = GetIndex(X,x_max)
			ix_min = GetIndex(X,x_min)
			iy_max = GetIndex(Y,y_max)
			iy_min = GetIndex(Y,y_min)

			extent = [X[ix_min],X[ix_max],Y[iy_min],Y[iy_max]]

		im = grid[index].imshow(np.transpose(dataBz[ix_min:ix_max,iy_min:iy_max]), #vmin = vminh, vmax = vmaxh, 
									  origin = "lower",
									  extent = extent,
									  vmax = vmaxh,
									  vmin = vminh,
									  cmap   = 'bwr')

		co = grid[index].contour(X[ix_min:ix_max],Y[iy_min:iy_max],np.transpose(dataAz[ix_min:ix_max,iy_min:iy_max]),4,
								 colors = ('k'),
								 origin = 'lower',
								 extent = extent,
								 linestyles = 'solid',
								 linewidths = 1)

		
		if j == n_conf - 1 :
			cb = grid[index].cax.colorbar(im)
			cb.set_label_text('$B_{z}$')


		#if index == 7 :
		#	cb1 = grid[index].cax.colorbar(im,ticks = [-1.0,0.0,1.0])
		#	tick_locator = ticker.MaxNLocator(nbins=1)
		#	cb1.set_label_text('$\\Delta N / N_{0}$')
		#	cb1.locator = tick_locator
			
		
		grid[index ].set_ylabel(ylabel[i])
		grid[index ].set_xlabel('$x c / \omega_{pi}$')

		index = index + 1 

## This affects all axes with option label_mode="L",
#grid.axes_llc.set_xticks([int(10.0*3.0*X[0,0]/(4.0*G))/10.0, 0.0 , int(10.0*3.0*X[-1,-1]/(4.0*G))/10.0])
#grid.axes_llc.set_yticks([int(10.0*3.0*Y[0,0]/(4.0*G))/10.0, 0.0 , int(10.0*3.0*Y[-1,-1]/(4.0*G))/10.0])
#
#grid[4].cax.colorbar(im)
#grid[3].axvline(x=6.0, ls = '--',color = 'k' )
plt.show()