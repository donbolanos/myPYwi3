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
path = '/media/sbolanos/BatDRIVE/HECKLE/PRO/T10/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
name = '400x600'
savepath = path+'bz_zoom_video/'

try : 
  os.mkdir(savepath)
except OSError:
  pass

liste =os.listdir(path)
os.chdir(savepath)

# .. set the needed parameter to plot the 2d field
domain    = None
shifts    = [0, 0]
bounds    = [-0.5,0.5]
colormap  = ['bwr', 64]
flines    = 10
ticks     = [20, 20]
subticks  = [2, 2]
figsize   = None #[6,9]
filetype  = 'png'
#
zoomx = [100,300]
zoomy = [250,350]

for folder in liste :

  if folder[0:len(name)] == name :

    run  = heckle.Heckle(path + folder, name)
    #
    # .. get the desired data giving time
    #data = run.GetN(time, "a")
    goodtime,grouptime = run.getTimeGroups([0,35])
    for time in goodtime :
      
      #
      # .. set the needed parameter to plot the 2d field
      filename  = '{:03d}'.format(int(time*5))
      title     = '$ t = {:.1f}$'.format(time)

      #dataX = run.GetB(time)[...,0]
      #dataY = run.GetB(time)[...,1]
      #dataZ = run.GetB(time)[...,2]
      data  = run.GetB(time)[...,2]
      dutu  = run.fourierFlux(time)
      #
      #normal = {'b0': 1.0e-8, 'n0': 2} # for hedp : b0 in megaGauss & n0 in cm-3
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
      
      #vectim0 = flatten(abs(im0.data))
      maxim0= np.max(np.max(np.abs(im0.data)))
      #bounds=[-maxim0,maxim0]
      # .. draw the plot
      plo = colp.Colp(coloraxis = (im0.axis[0][zoomx[0]:zoomx[1]],im0.axis[0][zoomy[0]:zoomy[1]]),
                      colordata = im0.data[zoomx[0]:zoomx[1],zoomy[0]:zoomy[1]],
                      bounds = bounds,
                      colormap = colormap,
                      contouraxis = (im1.axis[0][zoomx[0]:zoomx[1]],im1.axis[0][zoomy[0]:zoomy[1]]),
                      contourdata = im1.data[zoomx[0]:zoomx[1],zoomy[0]:zoomy[1]],
                      flines = flines,
                      arrowaxis = None,
                      arrowdata = None,
                      title = title,
                      labels = im0.labels,
                      ticks = ticks,
                      subticks = subticks,
                      figsize = figsize,
                      filetype = filetype,
                      filename = filename)
      #

os.system('bash ~/Documents/code/mypywi/videos/mp4.sh')
