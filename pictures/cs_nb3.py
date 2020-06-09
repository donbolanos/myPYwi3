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

path =  '/media/sbolanos/BatDRIVE/HECKLE/LMJ/nb3_T10/'
name = 'nb3'

liste =os.listdir(path)
t  = []
mr = []
mr_n = []
Az = []
beta = []


run  = heckle.Heckle(path , name)
#
# .. get the desired data giving time
#data = run.GetN(time, "a")
goodtime,grouptime = run.getTimeGroups([0,23.2])
vtime =  np.arange(0,np.max(goodtime),0.4)

for time in vtime :
	
	time = int(10.0*time)/10.0

	if time in goodtime and not (time in t) :
		plt.close("all")
		dutu = run.fourierFlux(time)

		Jz  = run.GetJ(time)[...,2] 
		Ez  = run.GetE(time)[...,2]
		
		Bx  = run.GetB(time)[...,0]
		By  = run.GetB(time)[...,1]
		Bz  = run.GetB(time)[...,2]
		
		B   = np.sqrt(Bx*Bx + By*By + Bz*Bz)
		N   = run.GetN(time,"e")

		# ------------- GET DATA -----------
		
		jmax = np.max(Jz[:300,:])
		
		if time < 20.5 :
			#ixmax,iymax = np.unravel_index(np.argmax(Jz[:300,:],axis = None), Jz.shape)
			#alternative
			iymax = 300
			ixmax = np.argmax(Jz[170:225,iymax]) + 170

			xmax = run.GetCoords(axis = 0)[ixmax]
			ymax = run.GetCoords(axis = 1)[iymax]
		else :
			print 'PASS'
			pass
		
		mr_tmp = Ez[ixmax,iymax]
		Az_tmp = dutu[1,1] - dutu[ixmax,iymax]


		Bo   = np.max(B)
		no   = np.max(N)
		beta = np.sqrt(no) / (Bo*Bo) 

		t.append(time)
		mr.append(mr_tmp)
		mr_n.append(mr_tmp*beta)
		Az.append(Az_tmp)

		print 'X point position : X = ',ixmax,' and Y = ',iymax, ' AND BETA = ',beta

		# ------------- GRAPH PART -----------#
		fig = plt.figure(1,(14,6))
		ax1 = plt.subplot2grid((2, 3), (0, 0), rowspan=2)
		ax2 = plt.subplot2grid((2, 3), (0, 1), rowspan=2)
		ax3 = plt.subplot2grid((2, 3), (0, 2))
		ax4 = plt.subplot2grid((2, 3), (1, 2), sharex = ax3) 

		ax3.plot(t,mr,'k')
		ax3.plot(t,mr_n,'g--')
		ax3.plot(time,mr_tmp,'ro', linewidth = 4.0)
		ax3.grid()
		ax3.set_xlim([0.0,50])
		#ax3.set_xlabel('Time')
		ax3.set_ylabel('Ez')
		plt.setp(ax3.get_xticklabels(), visible=False)
		# The y-ticks will overlap with "hspace=0", so we'll hide the bottom tick
		ax3.set_yticks(ax3.get_yticks()[1:])

		ax4.plot(t,Az,'k')
		ax4.plot(time,Az_tmp,'ro', linewidth = 2.0)
		ax4.grid()
		ax4.set_xlim([0.0,50])

		ax4.set_xlabel('Time')
		ax4.set_ylabel('Az')


		colordens = 'bwr'#['jet', 64]
		flines    = 4
		axis = np.array([run.GetCoords(axis = 0),run.GetCoords(axis = 1)])

		im1 = ax1.imshow(np.transpose(Jz),
						 aspect = 'equal',
		            	 interpolation = 'nearest',
		            	 cmap = colordens,
		            	 origin = 'lower',
		            	 extent = [axis[0][0],axis[0][-1],
		            	           axis[1][0],axis[1][-1]],
		            	 vmax = 1.0,
		            	 vmin = -1.0)
		            
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

		ax1.plot(xmax,ymax,'gx')

		ax1.set_title('Jz\nat t = {:.1f}'.format(time))

		im2 = ax2.imshow(np.transpose(Ez),
					aspect = 'equal',
		            interpolation = 'nearest',
		            cmap = colordens,
		            origin = 'lower',
		            extent = [axis[0][0],axis[0][-1],
		                      axis[1][0],axis[1][-1]],
		            vmax = 0.5,
		            vmin = -0.5)

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


		ax2.set_title('Ez\nat t = {:.1f}'.format(time))
		plt.tight_layout()
	


		plt.show()