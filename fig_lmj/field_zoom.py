# -*- coding: utf-8 -*-
"""
Comment on the script
...
"""

import sys
import numpy as np
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
flatten = lambda l: [item for sublist in l for item in sublist]
# .. load the run
path = '/media/sbolanos/BatDRIVE/HECKLE/LMJ/nb3_T10/'
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T10/400x600_vA_T10_t32-48/'
#path = '/media/sbolanos/DATA/BECKLE/'
name = ''
run  = heckle.Heckle(path, name)
#
# .. get the desired data giving time
#time = 53.0
time = 10.0
#data = run.GetN(time, "e")
#data = run.GetE(time)[...,2]
#data1 = run.GetB(time)[...,1]
data = run.GetB(time)[...,2]
#dataM = np.sqrt(data*data +data1*data1+data2*data2)
dutu = run.fourierFlux(time)
#data = run.GetJ(time)[...,1]
#
# .. set the needed parameter to plot the 2d field
domain    = None
shifts    = [-60, -60]
bounds    = [-0.6,0.6] #[0.1,0.35]
colormap  = ['bwr',64]#['Greys', 64]
labelcb   = '$n_e/n_0$'
flines    = 5
ticks     = None #[10,5]
labels    = ['x [mm]','y [mm]']
subticks  = [5,5]
figsize   = [7,4]
filetype  = None#'pdf'
filename  = name+'_t_'+str(time)
vlines    = None#[0.0,-4.0]
#
zoomx = [100,400]
zoomy = [150,450]

normal = {'b0': 1.0e-8, 'n0': 2} # for hedp : b0 in megaGauss & n0 in cm-3
#
# .. load the data for image
im0 = field.Field(run = run,
                  data = data,
                  domain = domain,
                  shifts = shifts)
#
# .. load the data for contours
im1 = field.Field(run = run,
                  data = dutu,
                  domain = domain,
                  shifts = shifts)

#.. passage d'unité adimensionnée à physique 
# constants :
A  = 1.
Z  = 1.
mu = 4*np.pi*1.0e-7
kb = 1.38e-23
mp = 1.67e-27 
qe = 1.6e-19
mi = A*mp 
b0 = 200  #express in T to denormalized
n0 = 0.5e26 #express in cm-3 to denormalized
  
va = b0/np.sqrt(mu*n0*mi/Z)
dp = np.sqrt(mi/(mu*n0*Z*qe*qe))
tc = mi/(Z*qe*b0)

im0.axis[0] = im0.axis[0]*dp*1e3 #express in mm 
im0.axis[1] = im0.axis[1]*dp*1e3 #express in mm 

im1.axis[0] = im1.axis[0]*dp*1e3 #express in mm 
im1.axis[1] = im1.axis[1]*dp*1e3 #express in mm 

#SKIP THE NORMALISED

#vectim0 = flatten(abs(im0.data))
maxim0= np.max(np.max(np.abs(im0.data)))
#bounds=[-maxim0,maxim0]
# .. draw the plot
plo = colp.Colp(coloraxis = (im0.axis[0][zoomx[0]:zoomx[1]],im0.axis[1][zoomy[0]:zoomy[1]]),
                colordata = im0.data[zoomx[0]:zoomx[1],zoomy[0]:zoomy[1]],
                bounds = bounds,
                colormap = colormap,
                labelcb = labelcb,
                contouraxis = (im1.axis[0][zoomx[0]:zoomx[1]],im1.axis[1][zoomy[0]:zoomy[1]]),
                contourdata = im1.data[zoomx[0]:zoomx[1],zoomy[0]:zoomy[1]],
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filename,
                vlines = None)
#
