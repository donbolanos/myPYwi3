
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
from matplotlib import rc
from mpl_toolkits.axes_grid1 import ImageGrid
sys.path.append("/home/sbolanos/Documents/code/Python/")

import runs.heckle  as heckle

#.. use latex fonts
rc('text',usetex = False)
rc('font',size = 16,family ='sans-serif')
rc('axes', labelsize = 'larger')
rc('mathtext',default = 'regular')

titles   = ['$Initialisation$','$Onset$', '$Reconnection$', '$End\ of\ MR$']
#time     = [[0.0,10,42.2,50.0],[0.0,33.0,34.8,42.2],[0.0,30.0,43.6,50.0]]
#ylabel   = ['T = 1\ny [mm]','T = 10\ny [mm]','T = 100\ny [mm]']

time = [3.0,10.,20.,30.,40.]
F = plt.figure(1,(12,3.))
grid = ImageGrid(F, 111,
                      nrows_ncols=(1,len(time)),
                      direction="row",
                      axes_pad= 0.12,#0.85, #0.12 normalement
                      add_all=True,
                      label_mode="L",
                      share_all=True,
                      cbar_location="right",
                      cbar_mode="edge",
                      cbar_size="5%",
                      cbar_pad=0.05,
                      )

path  = '/media/sbolanos/BatDRIVE/HECKLE/LMJ/nb3_T1/'
name  = '400x600'
vmaxh = 0.4
vminh = -vmaxh

#
zoomx = [100,400]
zoomy = [150,450]
  
run  = heckle.Heckle(path , name)

goodtime,grouptime = run.getTimeGroups([0,75])

#.. passage d'unité adimensionnée à physique 
# constants :
A  = 1.
Z  = 1.
mu = 4*np.pi*1.0e-7
kb = 1.38e-23
mp = 1.67e-27 
qe = 1.6e-19
mi = A*mp 
b0 = 200  #express in T to denormalized
n0 = 0.5e26 #express in cm-3 to denormalized
  
va = b0/np.sqrt(mu*n0*mi/Z)
dp = np.sqrt(mi/(mu*n0*Z*qe*qe))
tc = mi/(Z*qe*b0)

for t in goodtime:
	if t in time :
		index = np.where(t == time)[0][0]
		
		dataJz = run.GetB(t)[...,2]
		dataAz = run.fourierFlux(t)
		
		X = run.GetCoords(axis=0)
		#X = X - X[-1]/2
		X = X - 44.4 #with 3 NBeams
		X = X*dp*1e3 #express in mm

		Y = run.GetCoords(axis=1)
		Y = Y - Y[-1]/2
		Y = Y*dp*1e3    #express in mm
		
		extent = [X[zoomx[0]],X[zoomx[-1]],Y[zoomy[0]],Y[zoomy[-1]]]
		
		# PLOT 1st LINE Bx:
		im = grid[index].imshow(np.transpose(dataJz[zoomx[0]:zoomx[-1],zoomy[0]:zoomy[-1]]), vmin = vminh, vmax = vmaxh, 
								origin = "lower",
								extent = extent,
								cmap   = 'bwr')

		co = grid[index].contour(X[zoomx[0]:zoomx[-1]],Y[zoomy[0]:zoomy[-1]],np.transpose(dataAz[zoomx[0]:zoomx[-1],zoomy[0]:zoomy[-1]]),5,
						 		 colors = ('k'),
						 		 origin = 'lower',
						 		 extent = extent,
						 		 linestyles = 'solid',
						 		 linewidths = 1)
		title = '$t\Omega_{ci} = $'+'${:.1f}$'.format(t)
		grid[index].set_title(title)
		#grid[index].set_ylabel('$y c / \omega_{pi}$')
		grid[index].set_ylabel('y [mm]')

		
		
		#grid[index].set_xlabel('$x c / \omega_{pi}$')
		grid[index].set_xlabel('x [mm]')
		grid[index].set_xticks([-0.6,0.0,0.6,1.2])
		if index == 0 :
			cb1 = grid[0].cax.colorbar(im)
			
			cb1.set_label_text('$B_{z}$')
			
plt.show()