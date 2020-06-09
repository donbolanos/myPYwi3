#
'''
This script write in txt files the potential map in order to keep all in a unique folder
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
#
# .. load the run
path = '/media/sbolanos/BatDRIVE/HECKLE/TEST/store/gf0_T10_B0.92/'
name = 'plan'

savepath = path+'Az/'
try : 
  os.mkdir(savepath)
except OSError:
  pass

liste =os.listdir(path)

t  = []
Ez = []

#for folder in liste :
#
#  if folder[0:len(name)] == name :
#
#run  = heckle.Heckle(path + folder, name)
run =  heckle.Heckle(path, name)
#
# .. get the desired data giving time
#data = run.GetN(time, "a")
goodtime,grouptime = run.getTimeGroups([0.0,35.0])
#goodtime = np.arange(0.0,50.0,10.0)
for time in goodtime :
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
  
  fname = savepath + 'Az_t_' + str(time) +'.txt'
  with open(fname,'w') as fid :
    header = "#  x             y            Az\n"
    fid.write(header)
    for j,x in enumerate(im1.axis[0]):
      for i,y in enumerate(im1.axis[1]) : 
        sr ='{:1.4e}  {:1.4e} {:1.4e}\n'.format(x,y,im1.data[j,i])
        fid.write(sr)
  
    
  # --- PLOT PART --- #
  #plt.figure()
  #plt.imshow(im1.data)
  #plt.colorbar()
  #plt.title('t = ' + str(time)+'\n max = '+ str(int(np.max(dutu)*10.0)/10.0)+ ' and min = '+ str(int(np.min(dutu)*10.0)/10.0)+ ' and diff = '+ str(int((np.max(dutu) - np.min(dutu))*10.0)/10.0))

#plt.show()