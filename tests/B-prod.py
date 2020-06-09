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
path = '/media/sbolanos/BatDRIVE/HECKLE/TEST/store/las_CU_L16/'

name =  ''
run  = heckle.Heckle(path, name)
#
time_max = 11.0
goodtime,grouptime = run.getTimeGroups([0.0,time_max])

#paramètre pour l'analyse
Y = 170
X = 250 

x = np.arange(0.0,X*0.2 ,0.2)
Btot =[]

for time in goodtime : 
	
	# B : magnetic
	dataBx = run.GetB(time)[0:X,Y,0]
	dataBy = run.GetB(time)[0:X,Y,1]
	dataBz = run.GetB(time)[0:X,Y,2]

	dataB  = np.sqrt(dataBx**2 + dataBy**2 + dataBz**2)

	B_tmp  = np.sum(dataB)

	Btot.append(B_tmp)
"""
Btot ne représente que la somme du champ magnétique le long du rayon d'une des bulle. 
Le flux magnétique est donné par F = int B.dS. En 2D comme ici, dS est par conséquent 1D,
il représente donc une longueur. Dans les simulations HECKLE, le dS est donc le pas de résolution.
"""
plt.plot(goodtime,Btot)
plt.xlabel('time')
plt.ylabel('sum B')
plt.grid()
plt.show() 