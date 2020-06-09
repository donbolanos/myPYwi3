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
import draws.constellation   as cst

# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/750x500_ths/750x500_2half_6-12/'
path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_ths/T1_Ta0p2/400x600_T1_Ta0p2/'
#path = '/media/sbolanos/DATA/BECKLE/'
name = '400x600'
run  = heckle.Heckle(path, name)
#
# .. get the desired data giving time
time = 10.0

dutu = run.fourierFlux(time)
data = run.GetRSpecies(time,'p')
#
# .. set the needed parameter to plot the 2d field
domain    = None
shifts    = [0, 0]
bounds    = None
colormap  = ['bwr', 64]
flines    = 10
ticks     = [60, 30]
subticks  = [2, 2]
figsize   = [6,9]
filetype  = None#'pdf'
filename  = name+'_t_'+str(time)
#
normal = {'b0': 1.0e-8, 'n0': 2} # for hedp : b0 in megaGauss & n0 in cm-3
#

# .. load the data for contours
im1 = field.Field(run = run,
                  data = dutu,
                  domain = domain,
                  shifts = shifts)

# .. draw the plot
plo = cst.Constellation(coloraxis = None,
                        colordata = data,
                        bounds = bounds,
                        colormap = colormap,
                        contouraxis = im1.axis,
                        contourdata = im1.data,
                        flines = flines,
                        arrowaxis = None,
                        arrowdata = None,
                        labels = im1.labels,
                        ticks = ticks,
                        subticks = subticks,
                        figsize = figsize,
                        filetype = filetype,
                        filename = filename)
#
