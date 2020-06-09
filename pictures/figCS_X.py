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

def ly_CS(path,name,time,x): #lineout along y 
	
	liste =os.listdir(path)
	
	y  = []
	BZ = []
	JZ = []
	VY = []
	EZ = []
	
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
				VY  = run.GetV(time,'i')[x,250:350,1]
				
				return y,BX,BZ,JZ,VY

	print '/!\\ TIME NOT FIND /!\\'
	return y,BX,BZ,JZ,VY

def lx_CS(path,name,time,y): #lineout along x 
	liste =os.listdir(path)
	x  = []
	VX = []
	for folder in liste :

		if folder[0:len(name)] == name :
			
			run  = heckle.Heckle(path + folder, name)
			#
			# .. get the desired data giving time
			#data = run.GetN(time, "a")
			goodtime,grouptime = run.getTimeGroups([0,100])
	
			if time in goodtime  :
				x    = run.GetCoords(axis = 0)

				VX  = run.GetV(time,'i')[:,y,0]
				
				return x,VX

	print '/!\\ TIME NOT FIND /!\\'
	return x,VX

#.. use latex fonts
rc('text',usetex = False)
rc('font',size = 18,family ='sans-serif')
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
try : 
  os.mkdir(savepath)
except OSError:
  pass

tt,ez,az = readEZ(path)

x_l = 200
y_l = 300

shiftx = 0.0
shifty = -60.0

vtime = np.arange(0,70.0,0.2)

for timing in vtime : 
	timing = float(int(timing*10.)/10.)
	y,lb,gf,lj,lv = ly_CS(path,name,timing,x_l)
	x,lvx = lx_CS(path,name,timing,y_l)

	# .. shitf in order to center
	y = y + shifty
	x = x + shiftx

	# --- Plot 
	#fig = plt.figure(1,(18,6))
	#ax1 = plt.subplot2grid((1, 3), (0, 0))
	#ax2 = plt.subplot2grid((1, 3), (0, 1))
	#ax3 = plt.subplot2grid((1, 3), (0, 2))
	
	plt.title('$t/\Omega_{ci}^{-1} = $' +'{:.1f}'.format(timing)+' & $zc / \omega_{pi} = 0$')
	plt.plot( y, gf, label = 'GF'  ,linewidth = 2.0)
	plt.plot( y, lj, label = 'Jx'  ,linewidth = 2.0)
	plt.plot( y, lv, label = 'Vi_y',linewidth = 2.0)
	plt.xlim([y[0],y[-1]])
	plt.xlabel("y-axis")
	plt.grid()
	plt.legend(loc='best')
	plt.tight_layout()
	savename = savepath + 'cs_x_t{:.1f}.png'.format(timing)
	plt.savefig(savename)
	plt.close('all')
