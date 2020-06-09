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

#.. use latex fonts
rc('text',usetex = False)
rc('font',size = 16,family ='sans-serif')
rc('axes', labelsize = 'larger')
rc('mathtext',default = 'regular')
plt.rcParams['figure.constrained_layout.use'] = True

# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/WB_9/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T1/'
name = 'laser'

liste =os.listdir(path)
t  = []
Ez = []
beta = []
p0 = [1.0, 0.2, 0.0  ]
for folder in liste :

	if folder[0:len(name)] == name :
		
		run  = heckle.Heckle(path + folder, name)
		#
		# .. get the desired data giving time
		#data = run.GetN(time, "a")
		goodtime,grouptime = run.getTimeGroups([0,75])
		vtime =  np.arange(0,np.max(goodtime))

		for time in goodtime :
			if time in goodtime and not (time in t) :
				t.append(time)
				y    = run.GetCoords(axis = 1)

				BX  = run.GetB(time)[...,0]
				JZ  = run.GetJ(time)[...,2]

				lo  = BX[int(BX.shape[0]/2.),:]
				lz  = JZ[int(BX.shape[0]/2.),:]
				i_max = np.argmax(lo[int(1./4.*len(lo)):int(3./4.*len(lo))]) + int(1./4.*len(lo))
				i_min = np.argmin(lo[int(1./4.*len(lo)):int(3./4.*len(lo))]) + int(1./4.*len(lo))
				if i_max< i_min :
					tmp   = i_max
					i_max = i_min
					i_min = tmp
				elif i_max == i_min :
					raise NameError(' BX est constant !!!! ')
				#popt, pcov = curve_fit(fct, y[i_min:i_max], lo[i_min:i_max], p0 = p0)
				#p0 = popt
				##  ----------        GRAPH PART         -----------   ##
				plt.plot(y,lo,label = 'line-out')
				plt.plot(y,lz,label = 'jz')
				#plt.plot(y,fct(y,*popt), label = 'fit: B={:.1f}'.format(popt[1]))
				plt.grid()
				plt.legend(loc = 'best')
				plt.xlabel('y [norm]')
				plt.ylabel('Bx')
				#plt.ylim([-1.0,1.0])
				plt.title('t = {:.1f}'.format(time))
				plt.show()