#
#
import sys
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
#
# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x400_ohs/1hs_17-30/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
path = '/media/sbolanos/DATA/HECKLE/'
name = '400x400'
run  = heckle.Heckle(path, name)
#
# .. get the desired data giving time
time = 0.0
#data = run.GetN(time, "a")
dataX = run.GetB(time)[...,0]
dataY = run.GetB(time)[...,1]
dataZ = run.GetB(time)[...,2]
data  = np.sqrt(dataX**2 + dataY**2 + dataZ**2)
dutu  = run.fourierFlux(time)
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
figsize   = [6, 6]
filetype  = None#'pdf'
filename  = path+'B_'+name+'_t_'+str(time)
#
normal = {'b0': 1.0e-8, 'n0': 2} # for hedp : b0 in megaGauss & n0 in cm-3
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
                labels = im0.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filename)
#
