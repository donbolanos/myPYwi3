
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
rc('font',size = 12,family ='sans-serif')
rc('axes', labelsize = 'larger')
rc('mathtext',default = 'regular')

paths    =  ['/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T10/ILZ/n_1e+27_B_100.0/']
			 
titles   = ['A \nInitialisation','B \n Onset', 'C \n During MR', 'D \n End of MR']
#time     = [[0.0,10,42.2,50.0],[0.0,33.0,34.8,42.2],[0.0,30.0,43.6,50.0]]
#ylabel   = ['T = 1\ny [mm]','T = 10\ny [mm]','T = 100\ny [mm]']

time = [0.0,33.0,34.8,42.2]
F = plt.figure(1,(9,4))
grid = ImageGrid(F, 111,
                      nrows_ncols=(1,len(time)),
                      direction="row",
                      axes_pad=0.5,
                      add_all=True,
                      label_mode="L",
                      share_all=True,
                      cbar_location="right",
                      cbar_mode="each",
                      cbar_size="5%",
                      cbar_pad=0.05,
                      )

path  = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T10/'
name  = '400x600'
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
			dataNi = run.GetN(t,"p")[...]
			#dataBz = run.GetB(t)[...,2]
			#dataBx = run.GetB(t)[...,0]
			#dataJz = run.GetJ(t)[...,2]
			
			X = run.GetCoords(axis=0)
			X = X - X[-1]/2
			Y = run.GetCoords(axis=1)
			Y = Y - Y[-1]/2
			extent = [X[0],X[-1],Y[0],Y[-1]]
			
			# PLOT 1st LINE Bx:
			im = grid[index].imshow(np.transpose(dataNi), vmin = 0.0, 
										  origin = "lower",
										  extent = extent,
										  cmap   = 'jet')
			grid[index].set_title(titles[index])
			grid[index].set_ylabel('$n_{i}$\ny-axis')
			
			grid[index].set_xlabel('z-axis')
			grid[index].cax.colorbar(im)


plt.show()