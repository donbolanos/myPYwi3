# -*- coding: utf-8 -*-
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
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T10/'
name = '400x600'

liste =os.listdir(path)

t  = []
mean_n = []
mean_x = []
mean_y = []
mean_z = []
mean_P = []
mean_T = []

for folder in liste :

    if folder[0:len(name)] == name :

    	run = heckle.Heckle(path + folder, name)
    	   
        x   = run.GetCoords(axis = 0)
        nx  = len(x)
        y   = run.GetCoords(axis = 1)
        ny  = len(y)
        
        middle_y = ny/2
        middle_x = nx/2
        #
        # .. get the desired data giving time
        #data = run.GetN(time, "a")
        goodtime,grouptime = run.getTimeGroups([0,75])

        for time in goodtime :
          if not time in t :
            pxx    = run.GetPxx(time,"p") 
            pyy    = run.GetPyy(time,"p")
            pzz    = run.GetPzz(time,"p")
            N      = run.GetN(time,"p") 
            
            t.append(time)
            
            tmp = N[:,middle_y]
            lineout_n = N[:,middle_y]
            lineout_x = pxx[:,middle_y]
            lineout_y = pyy[:,middle_y]
            lineout_z = pzz[:,middle_y]
            lineout_P = (lineout_x + lineout_z + lineout_y ) / 3.
            lineout_T = lineout_P/lineout_n
            lineout_E = lineout_n**(5/2)/np.sqrt(lineout_P)
            if time in np.arange(70.):
                #plt.plot(lineout, label = '{:1.1f}'.format(time) )
                plt.plot(x,lineout_n,label = 'n : {:1.1f}'.format(time),linewidth = 3 )
                #plt.plot(x,lineout_x,label = 'pxx:{:1.1f}'.format(time) )
                #plt.plot(x,lineout_y,label = 'pyy:{:1.1f}'.format(time) )
                #plt.plot(x,lineout_z,label = 'pzz:{:1.1f}'.format(time) )
                plt.plot(x,lineout_T/100.0,label = 'T : {:1.1f}'.format(time) )
                plt.plot(x,lineout_P/10.0,label = 'P : {:1.1f}'.format(time) )

                plt.xlabel('x (norm)',fontsize=18)
                plt.ylabel('Density (norm)',fontsize=18)
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                #plt.ylim([0,3.0])
                plt.legend(loc = 'best')
                plt.grid()
                plt.show()

            mean_n.append(np.mean(lineout_n))
            mean_x.append(np.mean(lineout_x))
            mean_y.append(np.mean(lineout_y))
            mean_z.append(np.mean(lineout_z))
            mean_P.append(np.mean(lineout_P))
            mean_T.append(np.mean(lineout_T))


t,mean_n, mean_x, mean_y, mean_z, mean_P,mean_T = zip(*sorted(zip(t,mean_n, mean_x, mean_y, mean_z, mean_P,mean_T)))

xx,temps = np.meshgrid(x,t) 

#epsilon = 1e-6 # pourne pas diviser par 0
#eMAT = lineout_n**(5/2)/np.sqrt(lineout_P + epsilon) # matrix proportionnelle a l'émission bremstrahlung
#   ---------   SAVE PART  --------  # 
	
# -------------  PLOT SUM-UP ----------- #
plt.figure()

dens = mean_n #np.mean(lineout_n,axis = 1)
pres = mean_P #np.mean(lineout_P,axis = 1)
tmpt = mean_T #np.mean(lineout_T,axis = 1)
#emis = np.mean(eMAT     ,axis = 1)

plt.plot(t,dens,'k'  , label = 'N')
plt.plot(t,pres,'k:' , label = 'P')
plt.plot(t,tmpt,'k--' , label = 'T')
plt.plot(t,emis*1e2,'k--', label = '$N^{5/2}/P^{1/2}$ (x10)')

plt.grid()
plt.xlabel('t (norm)',fontsize=18)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#plt.ylim([0,np.max(dens)])
plt.ylim([0,0.05])

plt.legend(loc = 'best')

savename = path +'mean_emission_zoom'

plt.savefig(savename)
plt.show()

#   ---------   SAVE PART  --------  # 

nname = path+'density_i.txt'
pname = path+'pressure_i.txt'
ename = path+'emission.txt'
xname = path+'x.txt'
tname = path+'time.txt'

#np.savetxt(nname ,lineout_n , fmt='%1.4e')
#np.savetxt(pname ,lineout_P , fmt='%1.4e')
#np.savetxt(ename ,eMAT      , fmt='%1.4e')
#np.savetxt(xname ,xx        , fmt='%1.4e')
#np.savetxt(tname ,temps     , fmt='%1.4e')


# tenter avec un save txt
#with open(sname,'w') as fid :
#    header = '#timing[norm] Ez[norm] Az[norm] timing[norm] dA/dt[norm]\n'
#    fid.write(header)
#    for i in range(len(EE)):
#        sr = '{:1.4e}   {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n'.format(T[i],ez[i],Az[i],temps[i],EE[i])
#        fid.write(sr)
