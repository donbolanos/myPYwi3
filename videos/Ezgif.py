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
path = '/media/sbolanos/BatDRIVE/HECKLE/laser/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
name = 'las_OR_plan'
savepath = path+name+'/bz_video/'

try : 
  os.mkdir(savepath)
except OSError:
  pass

liste =os.listdir(path)
os.chdir(savepath)
t  = []
Ez = []

# .. set the needed parameter to plot the 2d field
domain    = None
shifts    = [0, 0]
bounds    = [-0.5,0.5]
colormap  = ['bwr', 64]
flines    = - np.arange(0.,1000.,2.0)[::-1]
ticks     = [20, 20]
subticks  = [2, 2]
figsize   = None #[6,9]
filetype  = 'png'

for folder in liste :

  if folder[0:len(name)] == name :

    run  = heckle.Heckle(path + folder, name)
    #
    # .. get the desired data giving time
    #data = run.GetN(time, "a")
    goodtime,grouptime = run.getTimeGroups([0,100.0])
    for time in goodtime :
      #dutu = run.fourierFlux(time)
      #data = run.GetJ(time)[...,1]
      #
      # .. set the needed parameter to plot the 2d field
      filename  = '{:03d}'.format(int(time*5))
      title     = '$ t = {:.1f}$'.format(time)

      #dataX = run.GetB(time)[...,0]
      #dataY = run.GetB(time)[...,1]
      #dataZ = run.GetB(time)[...,2]
      data  = run.GetB(time)[...,2]
      dutu  = run.fourierFlux(time)
      dutu  = dutu - np.max(dutu) + 0.05 # <<<=== APPLY A "GAUGE" 
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
      plo = colp.Colp(coloraxis = im0.axis,
                      colordata = im0.data,
                      bounds = bounds,
                      colormap = colormap,
                      contouraxis = im1.axis,
                      contourdata = im1.data,
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
