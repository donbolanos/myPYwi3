#
'''
This script try to fit a hyperbolic tangential function to the Bx component around the current sheet
To further determine a FWHM and fellow the evoution of the lateter.  
'''
#
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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
import potential as pot 

def fct(x,a,b,d):
	return a*np.tanh(b*(x - 60.0)) + d

def readEZ(path,filename ='ez_norm.txt',skiprows = 1): 
	fname = path + filename 
	data_e = np.loadtxt(fname, skiprows = skiprows)
	
	tt = data_e[:,0]
	ez = data_e[:,1]
	az = data_e[:,2]

	return tt, ez, az

def Bx_CS(path,name,time,x): 
	liste =os.listdir(path)
	y  = []
	BX = []
	BZ = []
	JZ = []
	beta = []
	for folder in liste :

		if folder[0:len(name)] == name :
			
			run  = heckle.Heckle(path + folder, name)
			#
			# .. get the desired data giving time
			#data = run.GetN(time, "a")
			goodtime,grouptime = run.getTimeGroups([0,100])
	
			if time in goodtime  :
				y    = run.GetCoords(axis = 1)[250:350]

				BX  = run.GetB(time)[x,250:350,0]
				BZ  = run.GetB(time)[x,250:350,2]
				JZ  = run.GetJ(time)[x,250:350,2]
				
				return y,BX,BZ,JZ

	print '/!\\ TIME NOT FIND /!\\'
	return y,BX,BZ,JZ

#.. use latex fonts
rc('text',usetex = False)
rc('font',size = 16,family ='sans-serif')
rc('axes', labelsize = 'larger')
rc('mathtext',default = 'regular')
#plt.rcParams['figure.constrained_layout.use'] = True

# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/WB_9/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
#path1 = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/45deg_T1/'
#name1 = '400x600'
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/45deg_T10/'
name = 'laser'
#path3 = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/45deg_T100/'
#name3 = 'laser'

savepath = path+'cs/'
#try : 
#  os.mkdir(savepath)
#except OSError:
#  pass

tt,ez,az = readEZ(path)

x_l = 200

vtime = np.arange(0,70.0)

for timing in vtime : 
	y,lb,gf,lj = Bx_CS(path,name,timing,x_l)
	# --- Plot 
	fig = plt.figure(1,(12,6))
	ax1 = plt.subplot2grid((1, 2), (0, 0))
	ax2 = plt.subplot2grid((1, 2), (0, 1))
	
	ax1.set_title('t = {:.1f} and GF_x={:.1e}'.format(timing,gf[len(gf)/2]/np.max(lb)))
	ax1.plot( y, lb, label = 'Bx')
	ax1.plot( y, gf, label = 'GF')
	ax1.plot( y, lj, label = 'Jz')
	ax1.grid()
	ax1.legend(loc='best')

	ax2.plot(tt,ez,'k')
	e_tmp = np.interp(timing,tt,ez)
	ax2.plot(timing,e_tmp,'ro', linewidth = 4.0)
	ax2.grid()
	ax2.set_xlabel('Time')
	ax2.set_xlim([0,60.0])
	plt.tight_layout()

	savename = savepath + 'cs_t{:.1f}.png'.format(timing)
	#plt.savefig(savename)
	plt.show()
