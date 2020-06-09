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
path = '/media/sbolanos/BatDRIVE/HECKLE/TEST/store/gf0_T10_B0.92/'
name = 'plan'

liste =os.listdir(path)

t  = []
Ez = []

#for folder in liste :
#
#  if folder[0:len(name)] == name : # SI plusieurs dossier
#run  = heckle.Heckle(path + folder, name)

run  = heckle.Heckle(path , name)
#
# .. get the desired data giving time
#data = run.GetN(time, "a")
goodtime,grouptime = run.getTimeGroups([0,30.])

for time in goodtime :
  if not time in t :
    #dutu = run.fourierFlux(time)
    data = run.GetE(time)[...,2]
    Bx   = run.GetB(time)[...,0]
    By   = run.GetB(time)[...,1]
    Bz   = run.GetB(time)[...,2]
    B    = np.sqrt(Bx*Bx + By*By + Bz*Bz)
    N    = run.GetN(time,"p") 
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
    
    Bo = np.max(B)
    no = np.max(N)
    enorm = im1.data[200,300] / (Bo*Bo / np.sqrt(no))
    Ez.append(enorm)
    t.append(time)

T, ez  = zip(*sorted(zip(t,Ez)))
pathAz = path + 'Az/'
az     = pot.getAzTemp(pathAz,T,401,601,[200,300],[1,300])
# RECONNECTION RATE CALCULATED FROM THE POTENTIAL
#righttime = np.arange(0.0,75,0.2)
#pathAz = path + 'Az/'
#temps,EE, Az = pot.getEzTemp(pathAz,t,401,601,[200,300],[1,300],returnPot = True)

fig, ax1 = plt.subplots(tight_layout = True)

ax2  = ax1.twinx()
lns1 = ax1.plot(T,ez, 'k' , label = 'Ez')
lns2 = ax2.plot(T,az, 'k:', label = 'A' )

ax1.set_xlabel('Time (norm)',fontsize=16)
ax1.set_ylabel('Z Electric Field (norm)',fontsize=16)
ax2.set_ylabel('Reconnected B Flux',fontsize=16)
plt.xlim([0,50.])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#plt.ylim([0,0.3])
# The labels in the same box
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc='best')
#ax1.legend(loc = 0)
#ax2.legend(loc = 0)
ax1.grid()

#   ---------   SAVE PART  --------  # 
savename = path+'ez_norm'
plt.savefig(savename)

plt.show()

sname = savename + '.txt'
with open(sname,'w') as fid :
    header = '#timing[norm] Ez[norm] Az[norm] \n'
    fid.write(header)
    for i in range(len(ez)):
        sr = '{:1.4e}   {:1.4e} {:1.4e}\n'.format(T[i],ez[i],az[i])
        fid.write(sr)
