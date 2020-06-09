# -*- coding: utf-8 -*-
"""
GRILLE DE COMPARAISON ENTRE cas coplanaire et 45deg ! 
POUR UNE SIMU reproduisant la reconnexion LULI2017 i.e avec plasma de cuivre
...
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
from matplotlib import rc
from mpl_toolkits.axes_grid1 import ImageGrid

# .. import the modulus
import runs.heckle  as heckle
import shapes.field as field
import draws.colp   as colp

#paramètre physique
mp = 1.67e-27 # exprimé en kg
qe = 1.6e-19  # exprimé en C
b0 = 200      # exprimé en T

tc = mp/(qe*b0) # gyropériode d'un proton

#.. use latex fonts
rc('text',usetex = True)
rc('font',size = 16,family ='sans-serif')
rc('axes', labelsize = 'larger')
rc('mathtext',default = 'regular')

#.. init pour graphs
cmap   = 'bwr'
flines = - np.arange(1.,150.,4.0)[::-1]
vmax   = +5.0
vmin   = -0.5

name = ''

paths    =  ['/media/sbolanos/BatDRIVE/HECKLE/TEST/store/las_CU_L18_n10/',
			 '/media/sbolanos/BatDRIVE/HECKLE/TEST/store/las_CU_L18_gf_n10_ts+/']
			 #'/media/sbolanos/BatDRIVE/HECKLE/PRO/anti_T1/ILZ/n_3e+26_B_200.0/']
			# '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T100/ILZ/n_1e+27_B_100.0/']

time     = np.arange(0.,86.,17.)
#titles   = ['Initialisation','Beginning \n of MR', 'During MR', 'End of MR']
titles = ['t = {:.1f} \n {:.1f} ns'.format(tt,tt*tc*1e9) for tt in time]

ylabel   = ['coplanar\ny [mm]','GF [45] \ny [mm]']


F = plt.figure(1,(13,7))
grid = ImageGrid(F, 111,
                 nrows_ncols=[len(paths),len(time)],
                 direction="row",
                 axes_pad=0.05,
                 add_all=True,
                 label_mode="L",
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="3%",
                 cbar_pad=0.05,
                 )
index = 0
for i,path in enumerate(paths):
	
	run  = heckle.Heckle(path, name)
	
	y    = run.GetCoords(axis = 1)
	x    = run.GetCoords(axis = 0)

	extent = [x[0],x[-1],y[0],y[-1]]
		

	for t in time:
		if i == 0 :
			grid[index].set_title(titles[index])
		

		dataAz    = run.fourierFlux(t)
		dataAz    = dataAz - np.max(dataAz) + 0.1 # <<<=== APPLY A "GAUGE" 
      
		#dataField = run.GetJ(t)[...,2]
		dataField = run.GetN(t,"e")
		
		im = grid[index].imshow(np.transpose(dataField), vmin = vmin, vmax = vmax, 
									  origin = "lower",
									  extent = extent,
									  cmap   = cmap)
		
		co = grid[index].contour(x,y,np.transpose(dataAz),flines,
								 colors = ('k'),
								 origin = 'lower',
								 extent = extent,
								 linestyles = 'solid',
								 linewidths = 1)

		#grid[index].set_yticks([2.0*Y[0,0]/(3.0*G),0.0,2.0*Y[-1,-1]/(3.0*G)])
		#grid[index].set_xticks([2.0*X[0,0]/(3.0*G),0.0,2.0*X[-1,-1]/(3.0*G)])
		grid[index].set_xlabel('z [mm]')
		grid[index].set_ylabel(ylabel[i])
		index = index + 1 

grid[0].cax.colorbar(im)

plt.show()