# .. import the modulus
import runs.heckle  as heckle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np

# .. load the run
#path = '/media/sbolanos/DATA/BECKLE/'
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/T1_Ta0p2/400x600_T1_Ta0p2/'

name = '400x600'
run  = heckle.Heckle(path, name)

##  ------ >> CONDITION DESIRED << ------- ##
xmin = 20.0
xmax = 40.0
ymin = 0.0
ymax = 5.0

##  ------   GET DATA    --------  ##
time  = 20.0 

dataI = run.GetIndexSpecies(time, "p") 
dataR = run.GetRSpecies(time    , "p")
dataV = run.GetVSpecies(time    , "p")

vx = []
vy = []
vz = []

for i in range(len(dataR)) : 
	if dataR[i,0] < xmax and dataR[i,0] > xmin and dataR[i,1] < ymax and dataR[i,1] > ymin and dataV[i,1]<0:
		
		vx.append(dataV[i,0])
		vy.append(dataV[i,1])
		vz.append(dataV[i,2])

plt.hist(vx,bins=100,histtype='step', label = 'vx')
plt.hist(vy,bins=100,histtype='step', label = 'vy')
plt.hist(vz,bins=100,histtype='step', label = 'vz')
# we need to normalize the data to 0..1 for the full
# range of the colormap
#fracs = N.astype(float)/N.max()
#norm = colors.Normalize(fracs.min(), fracs.max())
#
#for thisfrac, thispatch in zip(fracs, patches):
#    color = cm.viridis(norm(thisfrac))
#    thispatch.set_facecolor(color)
plt.legend(loc = 'best')
plt.show()