#
#
import sys
import numpy as np
#
# .. add the path
#sys.path.append("/home/roch/code/mypywi/runs")
#sys.path.append("/home/roch/code/mypywi/fourier")
#sys.path.append("/home/roch/code/mypywi/draws")
#
# .. import the modulus
from runs.heckle import Heckle
from fourier.fourier2d import Fourier2D
from draws.colp import Colp
#
#
# .. load the run
path = '/media/sbolanos/DATA/HOPPER/400x600_ths/'
name = 'zobi'
run  = Heckle(path, name)
#
# .. set the needed parameter to plot the 2d field
#plane     = "xy"
#domain    = None
field     = 'modulus'
display   = 'full'
#cut       = 0.0
time      = [0, 20]
#shifts    = [-4, -2]
bounds    = [0, 1000]
colormap  = ['jet', 64]
#flines    = None
ticks     = None
subticks  = None
figsize   = [9, 6]
filetype  = None
#
#
# get the list of times & timegroups
mytimes, mygroups = run.getTimeGroups(time)
#
#ds = float(run.dl) -> run.dl is already a float64
dt = float(mytimes[1]-mytimes[0])
#
numoftimes = mytimes.__len__()
numofcells = np.int_(run.ncells)+1
#
data = np.empty( np.insert(numofcells,0,numoftimes), dtype = complex, order = 'C')
#
# then fill the data for input fft
for it in range(mytimes.__len__()) :
    data[it, :].real = run.GetB(mytimes[it])[..., 2]
    data[it, :].imag = 0.0
#
axeX = run.GetCoords(axis=0)
axeY = run.GetCoords(axis=1)
axis = [axeX,axeY]

f    = Fourier2D(data, axis = axis, field = field, display = display)
#
#
# .. draw the plot
plo = Colp(colordata = [f.data[0]/float(run.dl[0]), f.data[1]/float(mytimes[1]-mytimes[0]), f.data[2]],
           bounds = bounds,
           colormap = colormap,
           contourdata = None,
           flines = None,
           arrowdata = None,
           labels = ['$k$', '$\omega$'],
           ticks = ticks,
           subticks = subticks,
           figsize = figsize,
           filetype = filetype)
#
