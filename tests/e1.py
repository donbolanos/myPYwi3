#
#
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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
path = '/media/sbolanos/TOURO/HECKLE/400x600_ths/0deg_11-20/'#'/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
name = '400x600'
run  = heckle.Heckle(path, name)
#
# .. get the desired data giving time
time = 20.0
data = run.GetV(time, "p")
dutu = run.fourierFlux(time)
#data = run.GetJ(time)[...,1]
#
# .. set the needed parameter to plot the 2d field
domain    = None
shifts    = [0, 0]
bounds    = None
colormap  = ['jet', 64]
flines    = 10
arrownorm = 1
ticks     = [20, 20]
subticks  = [2, 2]
figsize   = None
filetype  = None #'pdf'
filename  = name + '_t_' + str(time)
#
normal = {'b0': 1.0e-8, 'n0': 2} # for hedp : b0 in megaGauss & n0 in cm-3
#
# .. load the data for image
im0 = field.Field(run = run,
                  data = data[...,2],
                  domain = domain,
                  shifts = shifts)
#
axisX, axisY = np.meshgrid(im0.axis[0],im0.axis[1])
## .. load the data for contours
#im1 = field.Field(run = run,
#                  data = data[...,0],
#                  domain = domain,
#                  shifts = shifts)
#
## .. load the data for contours
#im2 = field.Field(run = run,
#                  data = data[...,1],
#                  domain = domain,
#                  shifts = shifts)

# .. draw the plot
#plo = colp.Colp(coloraxis = im0.axis,
#                colordata = im0.data,
#                bounds = bounds,
#                colormap = colormap,
#                contouraxis = None,
#                contourdata = None,
#                flines = flines,
#                arrowaxis = [axisX[::10,::10],axisY[::10,::10]],
#                arrowdata = [data[::10,::10,0],data[::10,::10,1]],
#                arrownorm = arrownorm,
#                labels = im0.labels,
#                ticks = ticks,
#                subticks = subticks,
#                figsize = figsize,
#                filetype = filetype,
#                filename = filename)

fig, ax = plt.subplots(num = 0, dpi = 100)

# .. draw the color image
colormap = cm.get_cmap('jet', 64)
col = np.hypot(data[...,0],data[...,1])
#im = ax.imshow(np.transpose(im0.data),
#               aspect = 'auto',
#               interpolation = 'nearest',
#               cmap = colormap,
#               origin = 'lower',
#               extent = [im0.axis[0][0], im0.axis[0][-1],
#                         im0.axis[1][0], im0.axis[1][-1]])
ar = ax.quiver(axisX[::5,::5], axisY[::5,::5],np.transpose(data[::5,::5,0]),np.transpose(data[::5,::5,1]),np.transpose(col[::5,::5]))

plt.show()