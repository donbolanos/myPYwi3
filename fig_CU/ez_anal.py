# coding: utf-8
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
from scipy.signal import find_peaks
#
# .. add the path ---> from Roch
#sys.path.append("/home/roch/code/mypywi/runs")
#sys.path.append("/home/roch/code/mypywi/shapes")
#sys.path.append("/home/roch/code/mypywi/draws")
#
# .. import the modulus
import runs.heckle  as heckle
import fields.fields as tf 
import shapes.field as field
import draws.colp   as colp
import potential as pot 

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
path = '/media/sbolanos/BatDRIVE/HECKLE/TEST/store/las_CU_L18_gf_n10/'
name = 'las'


savepath = path+'CS/'
try : 
  os.mkdir(savepath)
except OSError:
  pass

#tt,ez,az = readEZ(path)
tt = []
ez = []
ez_norm = []

run  = heckle.Heckle(path , name)
#
# .. get the desired data giving time
#data = run.GetN(time, "a")
goodtime,grouptime = run.getTimeGroups([0,200.0])
vtime =  np.arange(1.0,100.0,1.0) # np.arange(0,np.max(goodtime),1.0)

for time in vtime :
	
	#time = int(10.0*time)/10.0

	if time in goodtime  :
		plt.close("all")
		dutu = run.fourierFlux(time)

		Jz  = run.GetJ(time)[...,2] 
		Ez  = run.GetE(time)[...,2]
		Vy  = run.GetV(time,'i')[...,1]

		Bx  = run.GetB(time)[...,0]
		By  = run.GetB(time)[...,1]
		Bz  = run.GetB(time)[...,2]
		
		B   = np.sqrt(Bx*Bx + By*By + Bz*Bz)
		N   = run.GetN(time,"e")

		bmax= np.max(B[:,170:430])
		nmax= np.max(B[:,170:430]) 

#	# ------------- GET POSITION --------- //!\\ TO CHANGE !! #
#------>>>>>   TRASH 		
#	if time < 20.5 :
#
#	#ixmax,iymax = np.unravel_index(np.argmax(Jz[175:325,200:302],axis = None), Ez[175:325,200:302].shape)
#	#
#	#ixmax = ixmax + 175
#	#iymax = iymax + 200
#
#	#alternative
#	iymax = 300
#	ixmax = np.argmax(Jz[170:225,iymax]) + 170
#
#
#	iyr  = iymax
#	ixr  = ixmax + 15
#	
#	xr   = run.GetCoords(axis = 0)[ixr]
#	yr   = run.GetCoords(axis = 1)[iyr]
#	
##             TRASH <<<<< --------

	# ---- EXTRACTInG DATA of interest----- # 
	y    = run.GetCoords(axis = 1)
	x    = run.GetCoords(axis = 0)

	iymax = int(len(y)/2) # FONCTIONNE SSI la configuration est suffisamment symétrique pour avoir la couche de courant au centre
	
	er_h = Ez[  :  ,iymax] # horizontal line-out in the mid plan
	lj_h = Jz[  :  ,iymax] # horizontal line-out in the mid plan
	
	# ---- PEAK DETECTION ----- #
	ipeaks, _ = find_peaks( lj_h, prominence = 0.5, height = 0.4*np.max(lj_h))
	peaks = lj_h[ipeaks]
	xpeaks= x[ipeaks]

	# ---- CHOOSING THE PEAK ---- #
	# Due to the convection of the plasmoid : Ez got reversed in the plasmoid
	# CONSEQUENCE : take the peak where EZ is the flatest
	order = 1
	dx    = x[1] - x[0]
	dezdx = tf.deriv1D(er_h, order, dx)
	
	if len(ipeaks) == 0 : 
		iymax = int(len(y)/2)
		ixmax = int(len(x)/2)
	else : 
		ixmax = ipeaks[np.argmin(np.abs(dezdx[ipeaks]))]

	# --- Get the value ---- #
	ez_tmp = Ez[ixmax,iymax]

	ez.append(ez_tmp)
	ez_norm.append( ez_tmp/(bmax*bmax / np.sqrt(nmax)) )
	tt.append(time)

	print 'at t = ', time, ' -> Ez = ', ez_norm[-1],' AND normalised : ', bmax*bmax / np.sqrt(nmax)
	print 'X point position : X = ',ixmax,' and Y = ',iymax

	# ----- CHECKINGS : LINE-OUT -------- #
	
	iyr  = iymax
	ixr  = ixmax + 15
	
	xr   = run.GetCoords(axis = 0)[ixr]
	yr   = run.GetCoords(axis = 1)[iyr]
	
	lb   = Bx[ixmax,:]
	gf   = Bz[ixmax,:]	
	lj   = Jz[ixmax,:]
	vi   = Vy[ixmax,:]
	er   = Ez[ixmax,:]

	lb_r = Bx[ixr  ,:]
	gf_r = Bz[ixr  ,:]	
	lj_r = Jz[ixr  ,:]
	vi_r = Vy[ixr  ,:]
	er_r = Ez[ixr  ,:]

	lb_h = Bx[  :  ,iymax] # horizontal line-out in the mid plan
	gf_h = Bz[  :  ,iymax] # horizontal line-out in the mid plan	
	lj_h = Jz[  :  ,iymax] # horizontal line-out in the mid plan
	vi_h = Vy[  :  ,iymax] # horizontal line-out in the mid plan
	er_h = Ez[  :  ,iymax] # horizontal line-out in the mid plan


	xmax = run.GetCoords(axis = 0)[ixmax]
	ymax = run.GetCoords(axis = 1)[iymax]
	
	# --- Plot 
	fig  = plt.figure(1,(12,6))
	ax1  = plt.subplot2grid((1, 2), (0, 0))
	ax2  = plt.subplot2grid((1, 2), (0, 1))

	fig2 = plt.figure(2,(12,6))
	ax3  = plt.subplot2grid((1, 3), (0, 0))
	ax4  = plt.subplot2grid((1, 3), (0, 2))
	ax3b = plt.subplot2grid((1, 3), (0, 1))
	
	
	ax1.set_title('t = {:.1f} / X-point'.format(time))
	ax1.plot( y, lb, label = 'Bx')
	ax1.plot( y, gf, label = 'GF')
	ax1.plot( y, lj, label = 'Jz')
	ax1.plot( y, vi, label = 'Vy')
	ax1.plot( y, er, 'k', label = 'Ez')
	ax1.grid()
	ax1.set_xlim([40.0,80.0])
	ax1.legend(loc='best')

	ax2.set_title('t = {:.1f} / horizontal'.format(time))
	
	ax2.plot( x, lb_h, label = 'Bx')
	ax2.plot( x, gf_h, label = 'GF')
	ax2.plot( x, lj_h, label = 'Jz')
	ax2.plot( x, vi_h, label = 'Vy')
	ax2.plot( x, er_h, 'k', label = 'Ez')
	#pointe les pics
	ax2.plot( xpeaks, peaks, 'xg')
	ax2.plot( xpeaks, er_h[ipeaks], 'xg')
	ax2.plot( xmax  , lj_h[ixmax] , 'xr')

	ax2.grid()
	ax2.set_xlim([40.0,60.0])
	ax2.legend(loc='best')

	axis = np.array([run.GetCoords(axis = 0),run.GetCoords(axis = 1)])

	# affichage du champ : ===>> Jz <<=== 
	im3 = ax3.imshow(np.transpose(Jz),
					 aspect = 'equal',
	            	 interpolation = 'nearest',
	            	 cmap = 'bwr',
	            	 origin = 'lower',
	            	 extent = [axis[0][0],axis[0][-1],
	            	           axis[1][0],axis[1][-1]],
	            	 )
	            
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

	ax3.set_xlim([x[175],x[325]])
	ax3.set_ylim([y[200],y[400]])

	ax3.set_title('Jz\nat t = {:.1f}'.format(time))

	# affichage du champ : ===>> Ez <<=== 
	im3b = ax3b.imshow(np.transpose(Ez),
					 aspect = 'equal',
	            	 interpolation = 'nearest',
	            	 cmap = 'bwr',
	            	 origin = 'lower',
	            	 extent = [axis[0][0],axis[0][-1],
	            	           axis[1][0],axis[1][-1]],
	            	 vmax = np.max(Ez),
	            	 vmin = -np.max(Ez),
	            	 )
	            
	#plt.colorbar(im1)
	fig.colorbar(im3b, ax=ax3b, fraction=0.046, pad=0.04, format =  '%.2f')

	ax3b.contour(axis[0],
	            axis[1],
	            np.transpose(dutu),
	            4,
	            colors = ('k',),
	            origin = 'lower',
	            extent = [axis[0][0], axis[0][-1],
	                      axis[1][0], axis[1][-1]],
	            linestyles = 'solid',
	            linewidths = 2)

	ax3b.plot(xmax,ymax,'gx')
	ax3b.plot(xr,yr,'c+')

	ax3b.set_xlim([x[175],x[325]])
	ax3b.set_ylim([y[200],y[400]])

	ax3b.set_title('Ez\nat t = {:.1f}'.format(time))

	# trace MR rate vs time
	ax4.plot(tt,ez,'k')
	ax4.plot(tt,ez_norm,'r')
	e_tmp = np.interp(time,tt,ez)
	#ax4.plot(time,e_tmp,'ro', linewidth = 4.0)
	ax4.grid()
	ax4.set_xlabel('Time')
	ax4.set_xlim([0,150.0])
	ax4.set_ylim([-0.1,1.0])
	plt.tight_layout()

	#savename = savepath + 'cs_t{:.1f}.png'.format(time)
	#fig.savefig(savename)
	plt.show()

# ------ TO SAVE DATA AT THE END ----- #
#sname = path + 'ez_norm.txt'
#with open(sname,'w') as fid :
#    header = '#timing   Ez       Ez[norm] \n'
#    fid.write(header)
#    for i in range(len(ez)):
#        sr = '{:1.4e}   {:1.4e} {:1.4e}\n'.format(tt[i],ez[i],ez_norm[i])
#        fid.write(sr)