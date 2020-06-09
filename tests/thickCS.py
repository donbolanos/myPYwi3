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

def fwhm(x,yy):
	
	if len(x) != len(yy) : 
		print 'Error x and y don t have the same dimesion'
	
	imax  = np.argmax(yy)
	ymax  = np.max(yy)
	min_i = np.min(yy[:imax])
	min_s = np.min(yy[imax:])
	imin_i = np.argmin(yy[:imax])
	imin_s = np.argmin(yy[imax:]) + imax

	mid = (ymax + (min_s+min_i)/2.0)/2.0
	tmp_fwhm = x[np.max(np.where(yy[imin_i:imin_s]>=mid))] - x[np.min(np.where(yy[imin_i:imin_s]>=mid))]

	return min_i,min_s,ymax,mid,tmp_fwhm

def lo_CS(path,name,time,x): 
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

				#BX  = run.GetB(time)[x,250:350,0]
				#BZ  = run.GetB(time)[x,250:350,2]
				JZ  = run.GetJ(time)[x,250:350,2]
				
				return y,JZ #

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
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T10/'
name = '400x600'
titre = '$\\theta = 0^{\circ} $ & $\\beta = 10$'
#path3 = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/45deg_T100/'
#name3 = 'laser'
vfwhm = []
vjm   = []
vtime = []
xj    = []

tt,ez,az = readEZ(path)

x_l = 200

time = np.arange(0,50.0)

for timing in time :
	try : 
		y,lj = lo_CS(path,name,timing,x_l)
		min_i,min_s,ljmax,mid,tmp_fwhm = fwhm(y,lj)
		vfwhm.append(tmp_fwhm)
		vjm.append(ljmax)
		vtime.append(timing)
		xj.append(lj[int(len(lj)/2.0)])
	except ValueError : 
		break
	# --- Plot 
#	fig = plt.figure(1,(12,6))
#	ax1 = plt.subplot2grid((1, 2), (0, 0))
#	ax2 = plt.subplot2grid((1, 2), (0, 1))
#	
#	ax1.set_title('t = {:.1f} \n fwhm = {:.1f}'.format(timing,tmp_fwhm))
#	#ax1.plot( y, lb, label = 'Bx')
#	#ax1.plot( y, gf, label = 'GF')
#	ax1.plot( y, lj, label = 'Jz')
#	ax1.axhline(y=min_i, linestyle=':',color = 'k', lw=1.5)
#	ax1.axhline(y=min_s, linestyle=':',color = 'k', lw=1.5)
#	ax1.axhline(y=mid  , linestyle=':',color = 'k', lw=1.5)
#	ax1.grid()
#	ax1.legend(loc='best')
#
#	ax2.plot(tt,ez,'k')
#	e_tmp = np.interp(timing,tt,ez)
#	ax2.plot(timing,e_tmp,'ro', linewidth = 4.0)
#	ax2.grid()
#	ax2.set_xlabel('Time')
#	ax2.set_xlim([0,60.0])
#	plt.tight_layout()
#	plt.show()
#	savename = savepath + 'cs_t{:.1f}.png'.format(timing)
#	plt.savefig(savename)

fig = plt.figure(1 , constrained_layout = True)
ax  = plt.subplot(111)

ln1 = ax.plot(vtime,vfwhm,'k', label = 'fwhm')

ax2 = ax.twinx()

ln2 = ax.plot(vtime,vjm,'k:' , label = 'j_max')
ln3 = ax.plot(vtime,xj ,'k--', label = 'j_X')

# The labels in the same box
lns = ln1+ln2+ln3
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc='best')

ax.grid()
ax.set_ylabel('fwhm')
ax.set_xlabel('time')
ax2.set_ylabel('J max')

ax.set_title(titre)

plt.savefig( path + 'cs.png' )
plt.show()
