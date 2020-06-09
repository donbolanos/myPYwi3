#
'''
This script write in txt files the potential map in order to keep all in a unique folder
'''
#
import sys
import os
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
# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/new45/new45_27-40/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/T100_Ta20/400x600_T100_Ta02_t0-55/'
name = '400x600'

savepath = path+'Az/'
liste =os.listdir(path)

t  = []
Ez = []

run  = heckle.Heckle(path , name)

time = 0.0
dutu = run.fourierFlux(time)
#data = run.GetJ(time)[...,1]
#
# .. set the needed parameter to plot the 2d field
domain    = None
shifts    = [0, 0]
bounds    = None
colormap  = ['jet', 64]
flines    = 10
ticks     = [20, 20]
subticks  = [2, 2]
figsize   = [6, 9]
filetype  = None #'pdf'
filename  = name+'_t_'+str(time)
#
normal = {'b0': 1.0e-8, 'n0': 2} # for hedp : b0 in megaGauss & n0 in cm-3
#
## .. load the data for image
#im0 = field.Field(run = run,
#                  data = data,
#                  domain = domain,
#                  shifts = shifts)
##
## .. load the data for contours
im1 = field.Field(run = run,
                  data = dutu,
                  domain = domain,
                  shifts = shifts)

plo = colp.Colp(coloraxis = im1.axis,
                colordata = im1.data,
                bounds = bounds,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = im1.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filename)

#fname = savepath + 'Az_t_' + str(time) +'.txt'
#with open(fname,'w') as fid :
#  header = "#  x             y            Az\n"
#  fid.write(header)
#  for j,x in enumerate(im1.axis[0]):
#    for i,y in enumerate(im1.axis[1]) : 
#      sr ='{:1.4e}  {:1.4e} {:1.4e}\n'.format(x,y,im1.data[j,i])
#      fid.write(sr)
#
