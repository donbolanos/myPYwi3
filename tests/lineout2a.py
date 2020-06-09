#
#
import sys
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
#

rc('text', usetex = True)
rc('font', size=16)
rc('axes', labelsize='larger')
rc('mathtext', default='regular')

flatten = lambda l: [item for sublist in l for item in sublist]
# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/45deg_T10/laser_046r/'
#path = '/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
path = '/media/sbolanos/DATA/init_hkl/'
name = ''
run  = heckle.Heckle(path, name)
#
# .. get the desired data giving time
time = 0.0
data = run.GetN(time, "e")
#data = run.GetE(time)[...,2]
dataB = np.sqrt(run.GetB(time)[...,0]**2 + run.GetB(time)[...,1]**2 + run.GetB(time)[...,2]**2)
#data2 = run.GetB(time)[...,2]
#dataM = np.sqrt(data*data +data1*data1+data2*data2)
dutu = run.fourierFlux(time)
#data = run.GetJ(time)[...,1]
#
# .. set the needed parameter to plot the 2d field
domain    = None
shifts    = [0, 0]
bounds    = None #[-0.5,0.5]
colormap  = ['bwr', 64]
flines    = 5
ticks     = None#[60, 30]
labels    = ['$zc / \omega_{pi}$','$yc / \omega_{pi}$']
subticks  = [5,5]
figsize   = [6,7]
filetype  = None#'pdf'
filename  = name+'_t_'+str(time)
vlines = [0.0,-4.0]
#

y_lineout = 250

x = np.linspace(-75.,75,num=751)

plt.plot(x,data[y_lineout,:] ,label = '$n_e/n_0$')
plt.plot(x,dataB[y_lineout,:],label = '$|B|/B_0$')
plt.grid()
plt.xlabel(labels[1])
plt.legend(loc="best")
plt.tight_layout()
plt.show()

normal = {'b0': 1.0e-8, 'n0': 2} # for hedp : b0 in megaGauss & n0 in cm-3
#
# .. load the data for image
im0 = field.Field(run = run,
                  data = dataB,
                  domain = domain,
                  shifts = shifts)
#
# .. load the data for contours
im1 = field.Field(run = run,
                  data = dutu,
                  domain = domain,
                  shifts = shifts)

#vectim0 = flatten(abs(im0.data))
maxim0= np.max(np.max(np.abs(im0.data)))
#bounds=[-maxim0,maxim0]
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
                labels = labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filename,
                vlines = None)
#
