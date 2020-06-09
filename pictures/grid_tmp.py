
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
rc('text',usetex = True)
rc('font',size = 16,family ='sans-serif')
rc('axes', labelsize = 'larger')
rc('mathtext',default = 'regular')

titles   = ['$Initialisation$','$Onset$', '$Reconnection$', '$End\ of\ MR$']
#time     = [[0.0,10,42.2,50.0],[0.0,33.0,34.8,42.2],[0.0,30.0,43.6,50.0]]
#ylabel   = ['T = 1\ny [mm]','T = 10\ny [mm]','T = 100\ny [mm]']

time = [0.0,10.0,30.0,50.0]
F = plt.figure(1,(11,13))
grid = ImageGrid(F, 111,
                      nrows_ncols=(2,len(time)),
                      direction="row",
                      axes_pad=0.12,
                      add_all=True,
                      label_mode="L",
                      share_all=True,
                      cbar_location="right",
                      cbar_mode="edge",
                      cbar_size="5%",
                      cbar_pad=0.05,
                      )

path  = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T1/'
name  = 'laser'
vmaxh = 1.0
vminh = -vmaxh

liste =os.listdir(path)

for folder in liste :

  if folder[0:len(name)] == name :
    
	run  = heckle.Heckle(path + folder, name)
	goodtime,grouptime = run.getTimeGroups([0,75])
	for t in goodtime:
		if t in time :
			index = np.where(t == time)[0][0]
			
			dataBz = run.GetB(t)[...,2]
			dataBy = run.GetB(t)[...,1]
			dataBx = run.GetB(t)[...,0]
			dataAz = run.fourierFlux(t)
			
			X = run.GetCoords(axis=0)
			X = X - X[-1]/2
			Y = run.GetCoords(axis=1)
			Y = Y - Y[-1]/2
			extent = [X[0],X[-1],Y[0],Y[-1]]
			
			# PLOT 1st LINE Bx:
			im = grid[index].imshow(np.transpose(np.sqrt(dataBx**2 + dataBy**2)), vmin = 0.0, vmax = vmaxh, 
										  origin = "lower",
										  extent = extent,
										  cmap   = 'Reds')

			co = grid[index].contour(X,Y,np.transpose(dataAz),3,
							 colors = ('k'),
							 origin = 'lower',
							 extent = extent,
							 linestyles = 'solid',
							 linewidths = 1)

			grid[index].set_title(titles[index])
			grid[index].set_ylabel('$y c / \omega_{pi}$')
			
			# PLOT 2nd LINE Bz:
			im1 = grid[index + len(time)].imshow(np.transpose(dataBz), vmin = -0.25, vmax = 0.25, 
										  origin = "lower",
										  extent = extent,
										  cmap   = 'bwr')
			
			co = grid[index + len(time)].contour(X,Y,np.transpose(dataAz),3,
							 colors = ('k'),
							 origin = 'lower',
							 extent = extent,
							 linestyles = 'solid',
							 linewidths = 1)

			grid[index + len(time)].set_ylabel('$y c / \omega_{pi}$')
			
			
			grid[index + len(time)].set_ylabel('$y c / \omega_{pi}$')
			grid[index + len(time)].set_xlabel('$x c / \omega_{pi}$')
			if index == 0 :
				cb1 = grid[0].cax.colorbar(im)
				cb2 = grid[  len(time)].cax.colorbar(im1)


				cb1.set_label_text('$B_{xy}$')
				cb2.set_label_text('$B_{z}$')



plt.show()