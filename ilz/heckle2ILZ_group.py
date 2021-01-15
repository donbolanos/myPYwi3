#
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
#...initialisation plasma
n0   = 5.0e25 # electronic density         <<<<<<<<<
b0   = 200    # Tesla                       <<<<<<<<<
A    = 1
Z    = 1

# constant :
mu   = 4*np.pi*1.0e-7
kb   = 1.38e-23
mp   = 1.67e-27
qe   = 1.6e-19
mi   = A * mp

va   = b0/np.sqrt(mu*n0*mi/Z)  # alfven velocity
dp   = np.sqrt(mi/(mu*n0*Z*qe*qe))  # inertial length
tc   = mi/(Z*qe*b0)            # gyroperiod

# .. load the run
path = '/media/sbolanos/KINGSTON/nb2_T1/'
name = 'nlasers'
foname = 'lmj'

#create the folder 
path_dir = path + 'ILZ/'
try : 
    os.mkdir(path_dir)
except OSError :
    pass 

path_save=path_dir+'n_{:1.0e}_B_{:1.1f}/'.format(n0,b0)

try : 
    os.mkdir(path_save)
except OSError :
    pass 


liste = os.listdir(path)

tlist = [3.,10.,20.,30.,40.]#np.arange(0.0,50.0,10.) 

#for folder in liste :
#
#  if folder[0:len(name)] == name :
  
#       run  = heckle.Heckle(path + folder, '')
run  = heckle.Heckle(path, '')

goodtime,grouptime = run.getTimeGroups([0,201.])
#
# .. get the desired data giving time
for time in goodtime :
  if time in tlist :
    bz   = run.GetB(time)[...,2]
    by   = run.GetB(time)[...,1]
    bx   = run.GetB(time)[...,0]

    ez   = run.GetE(time)[...,2]
    ey   = run.GetE(time)[...,1]
    ex   = run.GetE(time)[...,0]

    x    = run.GetCoords(axis=0)
    y    = run.GetCoords(axis=1)

    Bz   = bz*b0   
    By   = by*b0   
    Bx   = bx*b0   
     
    Ez   = ez*va*b0   
    Ey   = ey*va*b0   
    Ex   = ex*va*b0   
     
    X    = x*dp    
    Y    = y*dp
    t 	 = time*tc

    filename = foname + '_t_' + str(time) + '.txt'
    fname    = path_save + filename
    with open(fname,'w+') as f :
    
        #header :
        header = '#x	y	Bx	By	Bz	Ex	Ey	Ez\n'
        f.write(header)

        for i, valX in enumerate(X):
            for j, valY in enumerate(Y):
                tmp = '{:1.4e}	{:1.4e}	{:1.4e}	{:1.4e}	{:1.4e}	{:1.4e}	{:1.4e}	{:1.4e}\n'.format(valX,valY,Bx[i,j],By[i,j],Bz[i,j],Ex[i,j],Ey[i,j],Ez[i,j]) 
                f.write(tmp)
