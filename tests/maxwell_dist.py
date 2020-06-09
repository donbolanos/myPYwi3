# .. import the modulus
import runs.heckle  as heckle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
from scipy.optimize import curve_fit

def maxwellian(x,N0,k):
	return N0*np.exp(-k*x)

# .. load the run
#path = '/media/sbolanos/DATA/BECKLE/'
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/T100_Ta20/400x600_T100_Ta20_t0-55/'

name = '400x600'
run  = heckle.Heckle(path, name)
#
# .. get the desired data giving time
time = 0.0

#data = run.GetN(time, "alpha")
#beta = run.GetB(time)
data = run.GetVSpecies(time,"alpha")

v = data[:,0]**2 + data[:,1]**2 + data[:,2]**2
label = 't= {:d}'.format(int(time))
N,bins, patches = plt.hist(v,bins=1000,histtype='step', log = True, label = label)
# we need to normalize the data to 0..1 for the full
# range of the colormap
#fracs = N.astype(float)/N.max()
#norm = colors.Normalize(fracs.min(), fracs.max())
#
#for thisfrac, thispatch in zip(fracs, patches):
#    color = cm.viridis(norm(thisfrac))
#    thispatch.set_facecolor(color)
xdata = (bins[1:]+bins[0:-1])/2
p0=[2e5,0.02]
popt, pcov = curve_fit(maxwellian, xdata[150:600], N[150:600],p0=p0)
plt.plot(xdata,maxwellian(xdata,*popt))

plt.legend(loc = 'best')
plt.show()