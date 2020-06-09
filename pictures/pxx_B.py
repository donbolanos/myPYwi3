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

def DATA_CS(path,name,time,xi=None): 
	liste =os.listdir(path)
	y  = []
	BX = []
	beta = []
#	for folder in liste :
#
#		if folder[0:len(name)] == name :
	
	#		run  = heckle.Heckle(path + folder, name)
	run  = heckle.Heckle(path, name)
	#
	# .. get the desired data giving time
	#data = run.GetN(time, "a")
	goodtime,grouptime = run.getTimeGroups([0,75])

	if time in goodtime  :
		y    = run.GetCoords(axis = 1)
		if xi is None : 
			xi =  int(len(run.GetCoords(axis = 0))/2.0)
		BX  = run.GetB(time)[xi,:,0]
		BY  = run.GetB(time)[xi,:,1]
		BZ  = run.GetB(time)[xi,:,2]
		B   = np.sqrt(BX**2 + BY**2 + BZ**2)
		

		Pyy_e = run.GetPyy(time,"e")[xi,:]
		#Pyy_p = 1./3.*(run.GetPyy(time,"p")[xi,:]+run.GetPxx(time,"p")[xi,:]+ run.GetPzz(time,"p")[xi,:])
		Pyy_p = run.GetPyy(time,"p")[xi,:]
		Pxx_p = run.GetPxx(time,"p")[xi,:]
		Pzz_p = run.GetPzz(time,"p")[xi,:]
		Np    = run.GetN(time,"p")[xi,:]
		Ne    = run.GetN(time,"e")[xi,:]
		Ni    = run.GetN(time,"i")[xi,:]
		return y,B, Pyy_e,Pyy_p,Pxx_p,Pzz_p,Np,Ne,Ni

	print '/!\\ TIME NOT FIND /!\\'
	return y,[],[],[]



#.. use latex fonts
rc('text',usetex = False)
rc('font',size = 16,family ='sans-serif')
rc('axes', labelsize = 'larger')
rc('mathtext',default = 'regular')
plt.rcParams['figure.constrained_layout.use'] = True

# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/WB_9/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T1_fullPart/'
#name = 'laser'
path = '/media/sbolanos/DATA/heckle/'
name = 'heckle'
#path3 = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T100/'
#name3 = 'laser'


savepath = path+'pressure/'
try : 
  os.mkdir(savepath)
except OSError:
  pass

#tt,ez,az = readEZ(path)

x_l = 200
#vtime = np.arange(0,50.)
vtime = [0.0]

for timing in vtime : 
	y,lb,pe,Pyy,Pxx,Pzz,ni,ne,nb = DATA_CS(path,name,timing)
	pb = lb*lb /2
	# --- Plot 
	fig = plt.figure(1,(12,6))
	ax1 = plt.subplot2grid((1, 2), (0, 0))
	ax2 = plt.subplot2grid((1, 2), (0, 1))
	
	ax1.set_title('t = {:.1f}'.format(timing))
	ax1.plot( y, lb , label = 'P_B')
	ax1.plot( y, pe , label = 'Pyy_e')
	ax1.plot( y, Pyy, label = 'Pyy')
	ax1.plot( y, ni , label = 'ni')
	ax1.plot( y, ne , label = 'ne')
	ax1.plot( y, pe / ne , label = 'nb')
	ax1.plot( y, Pyy/ ni , label = 'nb')
	#ax1.plot( y, nb , label = 'nb')
	#ax1.plot( y, pe+lb+Pyy, label = 'total')
	ax1.grid()
	#ax1.set_ylim([0,1.0])
	ax1.legend(loc='best')
	
	#ax2.plot(tt,ez,'k')
	#e_tmp = np.interp(timing,tt,ez)
	#ax2.plot(timing,e_tmp,'ro', linewidth = 4.0)
	#ax2.grid()
	#ax2.set_xlabel('Time')
	#ax2.set_xlim([0,50.0])
	
	#savename = savepath + 'pres_t{:.1f}.png'.format(timing)
	#plt.savefig(savename)
	plt.show()
	