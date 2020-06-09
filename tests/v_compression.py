# -*- coding: utf-8 -*-
#
'''
This script extract from the h5 file the data on the velocity 
In order to see any velocity which can affect the compressin of the magnetic field lines
'''
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
import potential as pot 
#
# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/WB_9/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/45deg_T10/'
name = 'laser'

liste =os.listdir(path)

t  = []
Ez = []
Vx_e = []
Vx_p = []
Vx_a = []
Vx_i = []

ysup = 5

for folder in liste :

  if folder[0:len(name)] == name :
    
        run  = heckle.Heckle(path + folder, name)
        #
        # .. get the desired data giving time
        #data = run.GetN(time, "a")
        x    = run.GetCoords(axis=0)
        xmin = len(x)/2 - ysup
        xmoy = len(x)/2 
        xmax = len(x)/2 + ysup 
        goodtime,grouptime = run.getTimeGroups([0,75])

        for time in np.arange(0,50) :
          if (not time in t) and (time in goodtime)  :
            t.append(time)
            # V : velocity 
            dataVx = run.GetV(time,'i')[xmoy,:,0]
            dataVy = run.GetV(time,'i')[xmoy,:,1]
            dataVz = run.GetV(time,'i')[xmoy,:,2]
            #Vx_e.append(dataVx_e[0,ymoy])#[ymin:ymax,:])
            
            # V : velocity 
            dataVx_p = run.GetV(time,'p')[...,0]
            #dataVy_i = run.GetV(time,'p')[...,1]
            #dataVz_i = run.GetV(time,'p')[...,2]
            #Vx_p.append(dataVx_p[0,ymoy])#[ymin:ymax,:])
            
            # V : velocity 
            dataVx_i = run.GetV(time,'i')[...,0]
            #dataVy_i = run.GetV(time,'p')[...,1]
            #dataVz_i = run.GetV(time,'p')[...,2]
            #Vx_i.append(dataVx_i[0,ymoy])#[ymin:ymax,:])
            
            # V : velocity 
            dataVx_a = run.GetV(time,"alpha")[...,0]
            #dataVy_a = run.GetV(time,"alpha")[...,1]
            #dataVz_a = run.GetV(time,"alpha")[...,2]
            #Vx_a.append(dataVx_a[0,ymoy])#[ymin:ymax,:])
            #

#t,Vx_e,Vx_p,Vx_i,Vx_a = zip(*sorted(zip(t,Vx_e,Vx_p,Vx_i,Vx_a)))

            plt.plot(dataVx,label = 'X')
            plt.plot(dataVy,label = 'Y')
            plt.plot(dataVz,label = 'Z')
            plt.grid()
            plt.xlabel('x')
            plt.ylabel('velocity (norm)')
            plt.legend(loc = 'best')
            plt.title('t = {:.1f}'.format(time))
            plt.show()