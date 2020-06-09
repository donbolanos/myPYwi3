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
import potential as pot 
#
# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/WB_9/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T1/'
name = 'laser'

liste =os.listdir(path)

t  = []
Ez = []

for folder in liste :

  if folder[0:len(name)] == name :
    
    run  = heckle.Heckle(path + folder, name)
    #
    # .. get the desired data giving time
    #data = run.GetN(time, "a")
    goodtime,grouptime = run.getTimeGroups([0,100])
    for time in goodtime :
      if not time in t :
        #dutu = run.fourierFlux(time)
        data = run.GetE(time)[...,2]
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
        ## .. load the data for contours
        im1 = field.Field(run = run,
                          data = data,
                          domain = domain,
                          shifts = shifts)
      
        Ez.append(im1.data[200,300])
        t.append(time)

T, ez = zip(*sorted(zip(t,Ez)))

# RECONNECTION RATE CALCULATED FROM THE POTENTIAL
righttime = np.arange(0.0,75,0.2)
pathAz = path + 'Az/'
temps,EE, Az = pot.getEzTemp(pathAz,t,401,601,[200,300],[1,300],returnPot = True)


plt.plot(T,ez, 'k',label = 'Ez')
plt.plot(temps,EE, label = 'dA/dt')

plt.xlabel('Time (norm)',fontsize=26)
plt.ylabel('Z Electric Field (norm)',fontsize=26)
plt.xlim([0,75.])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#plt.ylim([0,0.3])
plt.legend(loc = 'best')
plt.grid()

#   ---------   SAVE PART  --------  # 
savename = path+'ez_tmp'
plt.savefig(savename)

plt.show()

sname = savename + '.txt'
with open(sname,'w') as fid :
    header = '#timing[norm] Ez[norm] Az[norm] timing[norm] dA/dt[norm]\n'
    fid.write(header)
    for i in range(len(EE)):
        sr = '{:1.4e}   {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n'.format(T[i],ez[i],Az[i],temps[i],EE[i])
        fid.write(sr)
