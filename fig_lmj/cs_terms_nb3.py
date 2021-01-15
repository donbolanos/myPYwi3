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

	print('/!\\ TIME NOT FIND /!\\')
	return y,BX,BZ,JZ

#.. use latex fonts
rc('text',usetex = False)
rc('font',size = 16,family ='sans-serif')
rc('axes', labelsize = 'larger')
rc('mathtext',default = 'regular')
#plt.rcParams['figure.constrained_layout.use'] = True

# .. load the run
path = '/media/sbolanos/BatDRIVE/HECKLE/LMJ/nb2_T1/'
name = 'laser'


savepath = path+'CS/'
try : 
  os.mkdir(savepath)
except OSError:
  pass

tt,ez,az = readEZ(path)

run  = heckle.Heckle(path , name)
#
# .. get the desired data giving time
#data = run.GetN(time, "a")
goodtime,grouptime = run.getTimeGroups([0,40])
vtime =  np.arange(0,np.max(goodtime),0.2)

for time in vtime :
	
	time = int(10.0*time)/10.0

	if time in goodtime  :
		plt.close("all")
		dutu = run.fourierFlux(time)

		Jz  = run.GetJ(time)[...,2] 
		Ez  = run.GetE(time)[...,2]
		Vx  = run.GetV(time,'i')[...,0]
		Vy  = run.GetV(time,'i')[...,1]

		Bx  = run.GetB(time)[...,0]
		By  = run.GetB(time)[...,1]
		Bz  = run.GetB(time)[...,2]
		
		B   = np.sqrt(Bx*Bx + By*By + Bz*Bz)
		N   = run.GetN(time,"e")

	# ------------- GET POSITION -----------
		
	if time < 20.5 :
		#ixmax,iymax = np.unravel_index(np.argmax(Jz[:300,:],axis = None), Jz.shape)
		#alternative
		iymax = 300
		ixmax = np.argmax(Jz[170:225,iymax]) + 170

		xmax = run.GetCoords(axis = 0)[ixmax]
		ymax = run.GetCoords(axis = 1)[iymax]
		print('X point position : X = ',ixmax,' and Y = ',iymax)
	else :
		print('PASS')
		pass


	iyr = iymax
	ixr = ixmax + 15
	
	xr = run.GetCoords(axis = 0)[ixr]
	yr = run.GetCoords(axis = 1)[iyr]
	
	# ------------- LINE-OUT -----------
	
	lb   = Bx[ixmax,:]
	gf   = Bz[ixmax,:]	
	lj   = Jz[ixmax,:]
	vi   = Vy[ixmax,:]

	lb_x = By[:,iymax]
	gf_x = Bz[:,iymax]	
	lj_x = Jz[:,iymax]
	vi_x = Vx[:,iymax]


	lb_r = Bx[ixr  ,:]
	gf_r = Bz[ixr  ,:]	
	lj_r = Jz[ixr  ,:]
	vi_r = Vy[ixr  ,:]

	x    = run.GetCoords(axis = 0)
	y    = run.GetCoords(axis = 1)

	# --- Plot 
	fig = plt.figure(1,(12,6))
	ax1 = plt.subplot2grid((1, 2), (0, 0))
	ax2 = plt.subplot2grid((1, 2), (0, 1))

	fig2 = plt.figure(2,(12,6))
	ax3 = plt.subplot2grid((1, 2), (0, 0))
	ax4 = plt.subplot2grid((1, 2), (0, 1))
	
	ax1.set_title('t = {:.1f} / Y axis'.format(time))
	ax1.plot( y, lb, label = 'Bx')
	ax1.plot( y, gf, label = 'GF')
	ax1.plot( y, lj, label = 'Jz')
	ax1.plot( y, vi, label = 'Vy')
	ax1.grid()
	ax1.set_xlim([40.0,80.0])
	ax1.legend(loc='best')

	ax2.set_title('t = {:.1f} / X axis'.format(time))
	ax2.plot( x, lb_x, label = 'Bx')
	ax2.plot( x, gf_x, label = 'GF')
	ax2.plot( x, lj_x, label = 'Jz')
	ax2.plot( x, vi_x, label = 'Vy')
	ax2.grid()
	ax2.set_xlim([20.0,60.0])
	ax2.legend(loc='best')

	axis = np.array([run.GetCoords(axis = 0),run.GetCoords(axis = 1)])

	im3 = ax3.imshow(np.transpose(Bz),
					 aspect = 'equal',
	            	 interpolation = 'nearest',
	            	 cmap = 'bwr',
	            	 origin = 'lower',
	            	 extent = [axis[0][0],axis[0][-1],
	            	           axis[1][0],axis[1][-1]],
	            	 vmax = 1.0,
	            	 vmin = -1.0)
	            
	#plt.colorbar(im1)
	fig.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04, format =  '%.2f')

	ax3.contour(axis[0],
	            axis[1],
	            np.transpose(dutu),
	            4,
	            colors = ('k',),
	            origin = 'lower',
	            extent = [axis[0][0], axis[0][-1],
	                      axis[1][0], axis[1][-1]],
	            linestyles = 'solid',
	            linewidths = 2)

	ax3.plot(xmax,ymax,'gx')
	ax3.plot(xr,yr,'c+')

	ax3.set_title('Jz\nat t = {:.1f}'.format(time))

	ax4.plot(tt,ez,'k')
	e_tmp = np.interp(time,tt,ez)
	ax4.plot(time,e_tmp,'ro', linewidth = 4.0)
	ax4.grid()
	ax4.set_xlabel('Time')
	ax4.set_xlim([0,40.0])
	plt.tight_layout()

	savename = savepath + 'cs_t{:.1f}.png'.format(time)
	fig.savefig(savename)
	#plt.show()
