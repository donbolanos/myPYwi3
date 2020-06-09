#
'''
Calcul l'evolution temporel du parametre Beta 
'''
#
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
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
#.. use latex fonts
rc('text',usetex = False)
rc('font',size = 16,family ='sans-serif')
rc('axes', labelsize = 'larger')
rc('mathtext',default = 'regular')
plt.rcParams['figure.constrained_layout.use'] = True

# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/WB_9/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T5_n0.2/'
name = 'laser'
Tth  = 5.0
liste =os.listdir(path)

t  = []
Ez = []
beta = []
beta_e = []
n_m = []
b_m = []
for folder in liste :

  if folder[0:len(name)] == name :
    
        run  = heckle.Heckle(path + folder, name)
        #
        # .. get the desired data giving time
        #data = run.GetN(time, "a")
        goodtime,grouptime = run.getTimeGroups([0,75])

        for time in goodtime :
          if not time in t :
          	
            Bx   = run.GetB(time)[...,0]
            By   = run.GetB(time)[...,1]
            Bz   = run.GetB(time)[...,2]
            B    = np.sqrt(Bx*Bx + By*By + Bz*Bz)
            N    = run.GetN(time,"p")
            P    = 1.0/3.0 * (run.GetPxx(time,"p") + run.GetPyy(time,"p") + run.GetPzz(time,"p"))
            B0   = np.max(B)
            n0   = np.max(N)
            p_tmp= np.max(P)
            tmp_b = n0 / ( B0 * B0 / 2.0 )
            beta_e.append(p_tmp / (B0*B0/2.0))
            beta.append(tmp_b)
            t.append(time)
            n_m.append(n0)
            b_m.append(B0)

t,beta,beta_e,n_m,b_m = zip(*sorted(zip(t,beta,beta_e,n_m,b_m)))
beta          = np.array(beta) *Tth
## --------------- GRAPH PART -------------  ##

plt.plot(t,beta, label='mix')
plt.plot(t,beta_e,label  ='beta')
plt.plot(t,n_m, label='n_p')
plt.plot(t,b_m,label='B0')
plt.xlabel('Time')
plt.ylabel('$\\beta$')
#plt.ylim([0.3,5.0])
plt.legend()
plt.grid()
plt.show()