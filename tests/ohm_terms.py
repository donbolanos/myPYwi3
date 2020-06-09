#
#
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#
# .. add the path ---> from Roch
#sys.path.append("/home/roch/code/mypywi/runs")
#sys.path.append("/home/roch/code/mypywi/shapes")
#sys.path.append("/home/roch/code/mypywi/draws")
#
# .. import the modulus
import runs.heckle   as heckle
import shapes.field  as field
import fields.fields as tf
import draws.colp    as colp
import ohm           as ohm
#
#.. use latex fonts
rc('text',usetex = False)
rc('font',size = 18,family ='sans-serif')
rc('axes', labelsize = 'larger')
rc('mathtext',default = 'regular')

def plot_lineout( run, time, savepath = None, xmin = None, xmax = None) : 

	shiftx = 0.0
	shifty = -60.0
	x    = run.GetCoords(axis = 0) + shiftx
	y    = run.GetCoords(axis = 1) + shifty
	x1   = 400/2 - 20  #x-position decalee du x-point
	x2   = 400/2       #x-position du x-point

	if xmin is None :
		xmin = np.min(y)

	if xmax is None : 
		xmax = np.max(y)

	Bz   = run.GetB(time)[...,2]
	Ez   = run.GetE(time)[...,2]
	l1ez = Ez[x1,:]
	l2ez = Ez[x2,:]
	################# CALCUL DES DIFFERENTS TERMES DE LA LOI D OHM GENERALISEE ###################

	####  -------------  grad of electron stress tensor  -------------  ####
	divpe = ohm.electron_pressure(run, time)
	l1dpe = divpe[x1,:,2] # lineout  
	l2dpe = divpe[x2,:,2] # lineout 
	####  -------------  Laplacian of j, anomalous viscousity  -----------  ####
	hprst = ohm.hyperresistivity(run, time)
	l1hrt = hprst[x1,:,2] #lineout
	l2hrt = hprst[x2,:,2] #lineout
	####  -------------  Hall term, j x B  -------------  ####
	hall   = ohm.hall(run,time)
	l1hall = hall[x1,:,2]
	l2hall = hall[x2,:,2]
	####  -------------  v x B, fluid velocity (vi all the ion)  -------------  ####
	viB   = ohm.ion_ideal(run,time)
	l1viB = viB[x1,:,2]
	l2viB = viB[x2,:,2]
	####  -------------  resistive term, eta J -------------  ####
	rst   = ohm.resistivity(run,time)
	l1rst = rst[x1,:,2]
	l2rst = rst[x2,:,2]

	l2sum = l2dpe + l2hrt + l2hall + l2viB + l2rst 

	#plt.figure()
	#vmax = np.max(np.abs(Bz))
	#plt.imshow(Bz,vmax = vmax , vmin = -vmax , cmap = 'bwr')
	#plt.axhline(y=x2, color='k', lw=1, ls = '-.')
	#plt.axhline(y=x1, color='k', lw=1, ls = '-.')
	#plt.colorbar()
	#
	plt.figure()
	#plt.title('v-linout X-point\n t = {:.1f}'.format(time))

	plt.plot(y,l2ez  ,'k:', linewidth = 2., label = 'Ex   ')
	#plt.plot(y,l2dpe , label = 'divPe')
	plt.plot(y,l2hrt , label = 'hprst')
	plt.plot(y,l2hall, label = 'Hall ')
	plt.plot(y,l2viB , label = 'vi x B  ')
	#plt.plot(y,l2rst , label = 'rst  ')
	plt.xlim([xmin,xmax])
	plt.legend(loc = 'best')
	plt.grid()
	plt.xlabel('y-axis')
	plt.tight_layout()
	if savepath is not None :
		plt.savefig(savepath+'X_ohm_t{:.1f}.png'.format(time), format='png', dpi=300)
	
	plt.figure()
	plt.title('v-linout left-side of X-point\n t = {:.1f}'.format(time))

	plt.plot(y,l1ez  , 'k:', linewidth = 2., label = 'Ez   ')
	plt.plot(y,l1dpe , label = 'divPe')
	plt.plot(y,l1hrt , label = 'hprst')
	plt.plot(y,l1hall, label = 'Hall ')
	plt.plot(y,l1viB , label = 'viB  ')
	plt.plot(y,l1rst , label = 'rst  ')
	plt.xlim([xmin,xmax])

	plt.legend(loc = 'best')
	plt.grid()
	if savepath is not None :
		plt.savefig(savepath+'L_ohm_t{:.1f}.png'.format(time), format='png', dpi=300)
	
	plt.figure()
	plt.title('Ez at X-point\n t = {:.1f}'.format(time))
	plt.plot(y,l2sum , label = 'SUM')
	plt.plot(y,l2ez  , label = 'Ez ')
	plt.xlim([xmin,xmax])

	plt.legend(loc = 'best')
	plt.grid()

	if savepath is not None :
		plt.savefig(savepath+'Ez_t{:.1f}.png'.format(time), format='png', dpi=300)
	

# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/new0_28-42/'
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T10/'
name = '400x600'
savepath = path+'ohm/'

try : 
  os.mkdir(savepath)
except OSError:
  pass

liste =os.listdir(path)

#
# .. get the desired data giving time
#tv  = np.arange(32.0,36.0,0.2)
tv  = [34.8]
t_total = []
for folder in liste :

  if folder[0:len(name)] == name :
    
        run  = heckle.Heckle(path + folder, name)
        #
        # .. get the desired data giving time
        #data = run.GetN(time, "a")
        goodtime,grouptime = run.getTimeGroups([0,75])

        for t in tv : #goodtime[::5] :
          if t in goodtime and not(t in t_total):
			#plot_lineout( run, t, savepath, xmin = 50.0, xmax = 70.0)
			plot_lineout( run, t , xmin = -5.0, xmax = 5.0)
			plt.show()
			t_total.append(t)
			#plt.close('all')			
#plt.show(block = False)
#plt.close()
