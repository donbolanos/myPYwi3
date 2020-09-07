# -*- coding: utf-8 -*-
#
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
#

# .. import the modulus
import runs.heckle  as heckle
import shapes.field as field
import draws.colp   as colp
#
#
# .. load the run
path = '/media/sbolanos/BatDRIVE/HECKLE/laser/las_OR_plan/'

name =  ''
run  = heckle.Heckle(path, name)
#
time_max = 100.0
goodtime,grouptime = run.getTimeGroups([0.0,time_max])
times = sorted(goodtime) 
#... paramètre pour l'analyse
dp = 0.2 #resolution
#position centre bulle #1 
xp = int(50./dp) 
yp = int(50./dp) 

#position milieu boîte
xc = int(50./dp)
yc = int(75./dp)

#init flux
FL = [] # B-flux on left side 
FR = [] # B-flux on right side
FT = [] # B-flux at the top
FD = [] # B-flux at the buttom

BL = [] # B-flux on left side 
BR = [] # B-flux on right side
BT = [] # B-flux at the top
BD = [] # B-flux at the buttom

for time in times : 
	
	# B : magnetic
	#... gauche du centre
	dataBx_l = run.GetB(time)[0:xp,yp,0]
	dataBy_l = run.GetB(time)[0:xp,yp,1]
	dataBz_l = run.GetB(time)[0:xp,yp,2]

	dataB_l  = np.sqrt(dataBx_l**2 + dataBy_l**2 + dataBz_l**2)

	fL_tmp   = np.sum(dataB_l)*dp

	FL.append(fL_tmp)
	BL.append(np.sum(dataBy_l)*dp)

	#... droite du centre
	dataBx_r = run.GetB(time)[xp:-1,yp,0]
	dataBy_r = run.GetB(time)[xp:-1,yp,1]
	dataBz_r = run.GetB(time)[xp:-1,yp,2]

	dataB_r  = np.sqrt(dataBx_r**2 + dataBy_r**2 + dataBz_r**2)

	fR_tmp   = np.sum(dataB_r)*dp

	FR.append(fR_tmp)
	BR.append(np.abs(np.sum(dataBy_r)*dp))
	
	#... bas du centre
	dataBx_d = run.GetB(time)[xp,0:yp,0]
	dataBy_d = run.GetB(time)[xp,0:yp,1]
	dataBz_d = run.GetB(time)[xp,0:yp,2]

	dataB_d  = np.sqrt(dataBx_d**2 + dataBy_d**2 + dataBz_d**2)

	fD_tmp   = np.sum(dataB_d)*dp

	FD.append(fD_tmp)
	BD.append(np.abs(np.sum(dataBx_d)*dp))
	
	#... haut du centre
	dataBx_t = run.GetB(time)[xp,yp:yc,0]
	dataBy_t = run.GetB(time)[xp,yp:yc,1]
	dataBz_t = run.GetB(time)[xp,yp:yc,2]

	dataB_t  = np.sqrt(dataBx_t**2 + dataBy_t**2 + dataBz_t**2)

	fT_tmp   = np.sum(dataB_t)*dp

	FT.append(fT_tmp)
	BT.append(np.sum(dataBx_t)*dp)

	#... check the B-profil
	list_time = [0.0,60.0]
	if time in list_time :
		plt.figure()
		#Bpara profil 
		plt.subplot(131)
		plt.plot(run.GetB(time)[xp,0:yc,0],label = 'bx')
		plt.plot(run.GetB(time)[:,yp,1]   ,label = 'by')
		plt.xlabel('x ou y [cell]')
		plt.ylabel('B-field')
		plt.legend(loc ='best')
		plt.title('B para')
		
		#Bperp profil 
		plt.subplot(132)
		plt.plot(run.GetB(time)[xp,0:yc,1],label = 'by')
		plt.plot(run.GetB(time)[:,yp,0]   ,label = 'bx')
		plt.xlabel('x ou y [cell]')
		plt.ylabel('B-field')
		plt.legend(loc ='best')
		plt.title('B perp')
		
		#Bz profil 
		plt.subplot(133)
		plt.plot(run.GetB(time)[xp,0:yc,2],label = 'y-axis')
		plt.plot(run.GetB(time)[:,yp,2]   ,label = 'x-axis')
		plt.xlabel('x ou y [cell]')
		plt.ylabel('B-field')
		plt.legend(loc ='best')
		plt.title('B Z')
		
"""
Btot ne représente que la somme du champ magnétique le long du rayon d'une des bulle. 
Le flux magnétique est donné par F = int B.dS. En 2D comme ici, dS est par conséquent 1D,
il représente donc une longueur. Dans les simulations HECKLE, le dS est donc le pas de résolution.
"""
plt.figure()
plt.plot(times, FL,'b' ,label ='FL' ) 
plt.plot(times, FR,'r' ,label ='FR' ) 
plt.plot(times, FT,'g' ,label ='FT' ) 
plt.plot(times, FD,'m' ,label ='FD' )
plt.plot(times, BL,'b.-' ,label ='BL' ) 
plt.plot(times, BR,'r--' ,label ='BR' ) 
plt.plot(times, BT,'g:' ,label ='BT' ) 
plt.plot(times, BD,'m:' ,label ='BD' )
plt.legend(loc = 'best') 
plt.xlabel('time')
plt.ylabel('B-flux')
plt.grid()
plt.show() 