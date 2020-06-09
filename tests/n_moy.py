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
lineout = []
lineout_e = []

for folder in liste :

    if folder[0:len(name)] == name :

    	run  = heckle.Heckle(path + folder, name)
    	#
        # .. get the desired data giving time
        #data = run.GetN(time, "a")
        goodtime,grouptime = run.getTimeGroups([0,100])

        for time in goodtime :
          if not time in t :
            Ni    = run.GetN(time,"p") 
            Ne    = run.GetN(time,"e")
            nx,ny = np.shape(Ni)
            middle_y = ny/2
            middle_x = nx/2
            t.append(time)
            if len(lineout) == 0 :
            	lineout = [Ni[:,middle_y]]
            	lineout_e = [Ne[:,middle_y]]
            else  :
                lineout   = np.concatenate( (lineout, [Ni[:,middle_y]]), axis = 0 ) 
                lineout_e = np.concatenate( (lineout_e, [Ne[:,middle_y]]), axis = 0 ) 

x = run.GetCoords()

xx,T = meshgrid(x,t) 
#   ---------   SAVE PART  --------  # 

for i in range(len(t))[::50] :
	#plt.plot(lineout[i,:], label = '{:1.1f}'.format(t[i]) )
	plt.plot(x,lineout_e[i,:],label = '{:1.1f}'.format(t[i]) )

plt.xlabel('x (norm)',fontsize=26)
plt.ylabel('Density (norm)',fontsize=26)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#plt.ylim([0,0.3])
plt.legend(loc = 'best')
plt.grid()


#plt.savefig(savename)
plt.show()

#   ---------   SAVE PART  --------  # 

ename = path+'density_e.txt'
pname = path+'density_p.txt'
xname = path+'x.txt'
tname = path+'time.txt'

np.savetxt(ename ,lineout_e , fmt='%1.4e')
np.savetxt(pname ,lineout   , fmt='%1.4e')
np.savetxt(xname ,xx        , fmt='%1.4e')
np.savetxt(tname ,T         , fmt='%1.4e')


# tenter avec un save txt
#with open(sname,'w') as fid :
#    header = '#timing[norm] Ez[norm] Az[norm] timing[norm] dA/dt[norm]\n'
#    fid.write(header)
#    for i in range(len(EE)):
#        sr = '{:1.4e}   {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n'.format(T[i],ez[i],Az[i],temps[i],EE[i])
#        fid.write(sr)
