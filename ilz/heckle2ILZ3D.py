#
#
import sys
import os 
import numpy as np
import matplotlib.pyplot as plt
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
#
#
# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/750x500_ths/750x500_0-13/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
#path = '/media/sbolanos/BatDRIVE/HECKLE/750x500_ths/T100_Ta20/750x500_T100_Ta20_t74-82/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
path = '/media/sbolanos/BatDRIVE/HECKLE/ANDREE/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
dirf = '3D_single/'
name = ''
run  = heckle.Heckle(path+dirf, name)
#
# .. get the desired data giving time
time = 0.5
bz   = run.GetB(time)[...,2]
by   = run.GetB(time)[...,1]
bx   = run.GetB(time)[...,0]

ez   = run.GetE(time)[...,2]
ey   = run.GetE(time)[...,1]
ex   = run.GetE(time)[...,0]

x    = run.GetCoords(axis=0)
y    = run.GetCoords(axis=1)
z    = run.GetCoords(axis=2)

#      >>>>>>>      INITIALISATION PLASMA   <<<<<<<<<
n0   = 1.0e25 # electronic density			<<<<<<<<<
b0   = 200   # Tesla						<<<<<<<<<
A    = 1	  #								<<<<<<<<<
Z    = 1	  #								<<<<<<<<<
##  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# constant :

mu   = 4*np.pi*1.0e-7
kb   = 1.38e-23
mp   = 1.67e-27
qe   = 1.6e-19
mi   = A*mp

va = b0/np.sqrt(mu*n0*mi/Z)       # alfven velocity
dp = np.sqrt(mi/(mu*n0*Z*qe*qe))  # inertial length
tc = mi/(Z*qe*b0)                 # gyroperiod

Bz   = bz*b0   
By   = by*b0   
Bx   = bx*b0   
 
Ez   = ez*va*b0   
Ey   = ey*va*b0   
Ex   = ex*va*b0   
 
X    = x*dp    
Y    = y*dp
Z    = z*dp
t 	 = time*tc

#CREATE THE FOLDER TO SAVE THE FILE
path_dir = path + dirf + 'ILZ_PRO/'
try : 
	os.mkdir(path_dir)
except OSError :
	pass 

#path_save='/media/sbolanos/BatDRIVE/HECKLE/750x500_ths/750x500_ILZ/n_{:1.0e}_t_{:1.1e}/'.format(n0,t*1.0e9)
path_save=path_dir+'n_{:1.0e}_B_{:1.1f}/'.format(n0,b0)

try : 
	os.mkdir(path_save)
except OSError :
	pass 

filename = name + '_t_' + '{:1.1f}.txt'.format(time)
fname    = path_save + filename
with open(fname,'w+') as f :
	
	#header :#x		y		Bx		By		Bz		Ex		Ey		Ez
	header = '#x		y		z		Bx		By		Bz		Ex		Ey		Ez\n'
	f.write(header)

	for k, valZ in enumerate(Z):
		for j, valY in enumerate(Y):
			for i, valX in enumerate(X):
				tmp = '{:1.4e}	{:1.4e}	{:1.4e}	{:1.4e}	{:1.4e}	{:1.4e}	{:1.4e}	{:1.4e}	{:1.4e}\n'.format(valX,valY,valZ,Bx[i,j,k],By[i,j,k],Bz[i,j,k],Ex[i,j,k],Ey[i,j,k],Ez[i,j,k]) 
				f.write(tmp)

#... data for ppILZ
print 'X : n = ', len(X),' | min = ',np.min(X)*1e3,' mm | max = ', np.max(X)*1e3,' mm'
print 'Y : n = ', len(Y),' | min = ',np.min(Y)*1e3,' mm | max = ', np.max(Y)*1e3,' mm'
print 'Z : n = ', len(Z),' | min = ',np.min(Z)*1e3,' mm | max = ', np.max(Z)*1e3,' mm'


#... path integrated
dz = Z[1] - Z[0]

Bxdz = np.sum(Bx, axis = 2) * dz
Bydz = np.sum(By, axis = 2) * dz
Exdz = np.sum(Ex, axis = 2) * dz
Eydz = np.sum(Ey, axis = 2) * dz

#... visualisation : PLOTTING 
extent = np.array([X[0],X[-1],Y[0],Y[-1]])*1e3 #to express in mm

fig, [ax1, ax2] = plt.subplots(1, 2, figsize = (10,5),sharey=True)

ax1.set_xlabel('X [mm]')
ax1.set_ylabel('Y [mm]')
ax1.set_title('Bxdz')
im1 = ax1.imshow(np.transpose(Bxdz), extent = extent,
			   cmap = 'bwr',
			   origin = 'lower')
fig.colorbar(im1, ax=ax1)

ax2.set_xlabel('X [mm]')
ax2.set_ylabel('Y [mm]')
ax2.set_title('Bydz')
im2 = ax2.imshow(np.transpose(Bydz), extent = extent,
			   cmap = 'bwr',
			   origin = 'lower')
fig.colorbar(im2, ax=ax2)

fig1, [ax11, ax21] = plt.subplots(1, 2, figsize = (10,5),sharey=True)

ax11.set_xlabel('X [mm]')
ax11.set_ylabel('Y [mm]')
ax11.set_title('Exdz')
im11 = ax11.imshow(np.transpose(Exdz), extent =extent,
			   cmap = 'bwr',
			   origin = 'lower')
fig1.colorbar(im11, ax=ax11)


ax21.set_xlabel('X [mm]')
ax21.set_ylabel('Y [mm]')
ax21.set_title('Eydz')
im21 = ax21.imshow(np.transpose(Eydz), extent =extent,
			   cmap = 'bwr',
			   origin = 'lower')
fig1.colorbar(im21, ax=ax21)

plt.show()