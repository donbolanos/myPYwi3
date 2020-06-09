#
'''
This script plot the density of the ions and electron and display at which stages of the MR we are 
by plotting the ez(t) and Az(t)
'''
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
import runs.heckle  as heckle
import shapes.field as field
import draws.colp   as colp
import potential    as pot 
#
#.. use latex fonts
rc('text',usetex = False)
rc('font',size = 16,family ='sans-serif')
rc('axes', labelsize = 'larger')
rc('mathtext',default = 'regular')

def readEZ(path,filename ='ez_norm.txt',skiprows = 1): 
	fname = path + filename 
	data_e = np.loadtxt(fname, skiprows = skiprows)
	
	tt = data_e[:,0]
	ez = data_e[:,1]
	az = data_e[:,2]

	return tt, ez, az

def plotFig(run,time,tt,ez,az,vmin = None , vmax=None, t_max =None, savepath = None):
	
	plt.close('all')

	dutu = run.fourierFlux(time)

	N_p  = run.GetN(time,"p") 
	N_e  = run.GetN(time,"e")

	# ------------- GRAPH PART -----------#
	fig = plt.figure(1,(14,6))
	ax1 = plt.subplot2grid((2, 3), (0, 0), rowspan=2)
	ax2 = plt.subplot2grid((2, 3), (0, 1), rowspan=2)
	ax3 = plt.subplot2grid((2, 3), (0, 2))
	ax4 = plt.subplot2grid((2, 3), (1, 2), sharex = ax3) 

	ax3.plot(tt,ez,'k')
	e_tmp = np.interp(time,tt,ez)
	ax3.plot(time,e_tmp,'ro', linewidth = 4.0)
	ax3.grid()
	if t_max is not None :
		ax3.set_xlim([0.0,t_max])
	#ax3.set_xlabel('Time')
	ax3.set_ylabel('Ez')
	plt.setp(ax3.get_xticklabels(), visible=False)
	# The y-ticks will overlap with "hspace=0", so we'll hide the bottom tick
	ax3.set_yticks(ax3.get_yticks()[1:])

	ax4.plot(tt,az,'k')
	a_tmp = np.interp(time,tt,az)
	ax4.plot(time,a_tmp,'ro', linewidth = 2.0)
	ax4.grid()
	if t_max is not None :
		ax4.set_xlim([0.0,t_max])

	ax4.set_xlabel('Time')
	ax4.set_ylabel('Az')


	colordens = 'jet'#['jet', 64]
	flines    = 8
	axis = np.array([run.GetCoords(axis = 0),run.GetCoords(axis = 1)])

	im1 = ax1.imshow(np.transpose(N_p),
				aspect = 'equal',
	            interpolation = 'nearest',
	            cmap = colordens,
	            origin = 'lower',
	            extent = [axis[0][0],axis[0][-1],
	                      axis[1][0],axis[1][-1]],
	            vmin = vmin,
	            vmax = vmax)
	#plt.colorbar(im1)
	fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04, format =  '%.2f')

	ax1.contour(axis[0],
	            axis[1],
	            np.transpose(dutu),
	            flines,
	            colors = ('k',),
	            origin = 'lower',
	            extent = [axis[0][0], axis[0][-1],
	                      axis[1][0], axis[1][-1]],
	            linestyles = 'solid',
	            linewidths = 2)


	ax1.set_title('proton density\nat t = {:.1f}'.format(time))

	im2 = ax2.imshow(np.transpose(N_e),
				aspect = 'equal',
	            interpolation = 'nearest',
	            cmap = colordens,
	            origin = 'lower',
	            extent = [axis[0][0],axis[0][-1],
	                      axis[1][0],axis[1][-1]],
	            vmin = vmin,
	            vmax = vmax)
	fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04, format =  '%.2f')

	ax2.contour(axis[0],
	            axis[1],
	            np.transpose(dutu),
	            flines,
	            colors = ('k',),
	            origin = 'lower',
	            extent = [axis[0][0], axis[0][-1],
	                      axis[1][0], axis[1][-1]],
	            linestyles = 'solid',
	            linewidths = 2.0)


	ax2.set_title('electron density\nat t = {:.1f}'.format(time))
	plt.tight_layout()
	
	if savepath is not None : 
		savename = savepath + 'dens_t{:.1f}.png'.format(time)
		plt.savefig(savename)

	return

# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/WB_9/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/30deg_T10/'
name = '400x600'


savepath = path+'density/'
try : 
  os.mkdir(savepath)
except OSError:
  pass


tt,ez,az = readEZ(path)

liste =os.listdir(path)
t  = []
Ez = []
beta = []

for folder in liste :

  if folder[0:len(name)] == name :
    
    run  = heckle.Heckle(path + folder, name)
    #
    # .. get the desired data giving time
    #data = run.GetN(time, "a")
    goodtime,grouptime = run.getTimeGroups([0,75])
    vtime =  np.arange(0,np.max(goodtime))

    for time in vtime :
    	if time in goodtime and not (time in t) :
          	
		plotFig(run,time,tt,ez,az, t_max= 70.0, savepath = None)
		t.append(time)
		plt.show()