
# .. import the modulus
import runs.heckle  as heckle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np

# .. load the run
#path = '/media/sbolanos/DATA/BECKLE/'
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/T100_Ta20/400x600_T100_Ta20_t0-55/'

name = '400x600'
run  = heckle.Heckle(path, name)
#
# .. get the desired data giving time
goodtime,grouptime = run.getTime([0,100],'species.h5')
e = []
e_back = []
mag = []
mu0 = 4*np.pi*1e-7
#data = run.GetN(time, "alpha")
for time in goodtime : 
	
	B = run.GetB(time)

	data = run.GetVSpecies(time,"p")
	dataa= run.GetVSpecies(time,"b")

	masse= run.GetMass("p")
	massA= run.GetMass("b")

	v      = 0.5 * masse * np.sum(data[:,0]**2 + data[:,1]**2 + data[:,2]**2)
	va     = 0.5 * massA * np.sum(dataa[:,0]**2 + dataa[:,1]**2 + dataa[:,2]**2)
	tmp_b  = 0.5 * np.sum(B[...,0]**2 + B[...,1]**2 + B[...,2]**2)*(0.15*0.15*0.15)/ mu0
	
	e.append(v)
	e_back.append(va)

	mag.append(tmp_b)

T, Ep, Ea, Emag = zip(*sorted(zip(goodtime,e,e_back,mag)))
#T,  Emag = zip(*sorted(zip(goodtime,mag)))
Etot = np.array(Ep)+np.array(Emag)+np.array(Ea)

plt.plot( T, Ep, label = 'proton')
plt.plot( T, Ea, label = 'back ')
plt.plot( T, Emag,label= 'mag')
plt.plot( T, Etot, label = 'E Tot ')
plt.xlabel('TIME')
plt.ylabel('sum(v2)')
plt.grid()
plt.legend(loc = 'best')
plt.show()