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
	beta = []
	for folder in liste :

		if folder[0:len(name)] == name :
			
			run  = heckle.Heckle(path + folder, name)
			#
			# .. get the desired data giving time
			#data = run.GetN(time, "a")
			goodtime,grouptime = run.getTimeGroups([0,75])
	
			if time in goodtime  :
				y    = run.GetCoords(axis = 1)

				BX  = run.GetJ(time)[x,:,2]
				return y,BX
	print '/!\\ TIME NOT FIND /!\\'
	return y,BX
#.. use latex fonts
rc('text',usetex = False)
rc('font',size = 16,family ='sans-serif')
rc('axes', labelsize = 'larger')
rc('mathtext',default = 'regular')
plt.rcParams['figure.constrained_layout.use'] = True

# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/WB_9/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
path1 = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/45deg_T1/'
name1 = '400x600'
path2 = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/45deg_T10/'
name2 = 'laser'
path3 = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/45deg_T100/'
name3 = 'laser'

x_l = 200
vtime = np.arange(0,50.0)

for timing in vtime : 
	y1,lb1 = Bx_CS(path1,name1,timing,x_l)
	y2,lb2 = Bx_CS(path2,name2,timing,x_l)
	y3,lb3 = Bx_CS(path3,name3,timing,x_l)
	
	# --- Plot 
	plt.title('t = {:.1f}'.format(timing))
	plt.plot( y1[250:350], lb1[250:350], label = 'T1')
	plt.plot( y2[250:350], lb2[250:350], label = 'T10')
	plt.plot( y3[250:350], lb3[250:350], label = 'T100')
	plt.grid()
	plt.legend(loc='best')
	plt.show()
