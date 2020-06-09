#
#
import sys
import numpy as np
from mayavi import mlab

# .. import the modulus
import runs.heckle  as heckle
import shapes.field as field
import draws.colp   as colp
#

# .. load the run
path =  '/media/sbolanos/BatDRIVE/HECKLE/ANDREE/3D_d/'
name = ''
run  = heckle.Heckle(path, name)

#
# .. get the desired data giving time
time = 0.2

Bx = run.GetB(time)[...,0]
By = run.GetB(time)[...,1]
Bz = run.GetB(time)[...,2]

dataNe = run.GetN(time,'e')

obj = mlab.flow(Bx,By,Bz)
mlab.outline()
mlab.show()
