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
#path1 = '/media/sbolanos/BatDRIVE/HECKLE/PRO/T1/laser_008/'
#path1 = '/media/sbolanos/DATA/heckle/'
path2 = '/media/sbolanos/BatDRIVE/HECKLE/PRO/T10/400x600_T10/'
#path3 = '/media/sbolanos/BatDRIVE/HECKLE/PRO/T1/400x600_T1/'
path3 = '/media/sbolanos/BatDRIVE/HECKLE/TEST/gf30_T10/gf30_T10_v2/'
name  = ''
#run1  = heckle.Heckle(path1, name)
run2  = heckle.Heckle(path2, name)
run3  = heckle.Heckle(path3, name)
#
# .. get the desired data giving time
time = 0.0

#data = run.GetN(time, "e")
#data = run.GetE(time)[...,2]
#dataB1 = run1.GetB(time)[...,0]
dataB2 = run2.GetB(time)[...,0]
dataB3 = run3.GetB(time)[...,0]

#dataN1 = run1.GetN(time,"e")
dataN2 = run2.GetN(time,"p")
dataN3 = run3.GetN(time,"p")

dataV2 = run2.GetV(time,"i")[...,1]
dataV3 = run3.GetV(time,"i")[...,1]

dataV2x = run2.GetV(time,"i")[...,0]
dataV3x = run3.GetV(time,"i")[...,0]

Y = run2.GetCoords(axis=1)
Y = Y - np.max(Y)/2.0
y = np.arange(len(dataN2[200,:]))*0.2 - 0.5*0.2 * len(dataN2[200,:])

X = run2.GetCoords(axis=0)
X = X - np.max(X)/2.0
x = np.arange(len(dataN2[:,300]))*0.2 - 0.5*0.2 * len(dataN2[:,300])



#dataM = np.sqrt(data*data +data1*data1+data2*data2)
#dutu = run.fourierFlux(time)
#data = run.GetJ(time)[...,1]
#
# .. set the needed parameter to plot the 2d field
domain    = None
shifts    = [-40, -60]
bounds    = None #[-0.5,0.5]
colormap  = ['jet', 64]
flines    = 5
ticks     = None#[60, 30]
labels    = ['$zc / \omega_{pi}$','$yc / \omega_{pi}$']
subticks  = [5,5]
figsize   = [6,7]
filetype  = None#'pdf'
filename  = name+'_t_'+str(time)
vlines = [0.0,-4.0]
#

#y_lineout = 210

#x = np.linspace(-40.,40.,num=401)
#.. fig 1
plt.figure()
#plt.plot(dataB1[200,:] ,label = 'no 1')
plt.plot(y,dataB2[200,:] ,label = 'chap. 2')
plt.plot(y,dataB3[200,:] ,label = 'chap. 3')
plt.grid()
plt.xlabel(labels[1])
plt.ylabel('$B_x$')
plt.legend(loc="best")
plt.xlim([-15,15])
plt.tight_layout()

#.. fig 2
plt.figure()
#plt.plot(dataN1[200,:] ,label = 'no 1')
plt.plot(y,dataN2[200,:] ,label = 'chap. 2')
plt.plot(y,dataN3[200,:] ,label = 'chap. 3')
plt.grid()
plt.title('density')
plt.xlabel(labels[1])
plt.legend(loc="best")
plt.tight_layout()


#.. fig 3
plt.figure()
#plt.plot(dataN1[200,:] ,label = 'no 1')
plt.plot(y,dataV2[200,:] ,label = 'chap. 2')
plt.plot(y,dataV3[200,:] ,label = 'chap. 3')
plt.grid()
plt.title('vitess')
plt.xlabel(labels[1])
plt.legend(loc="best")
plt.tight_layout()


#.. fig 3
plt.figure()
#plt.plot(dataN1[200,:] ,label = 'no 1')
plt.plot(x,dataV2x[:,300] ,label = 'chap. 2')
plt.plot(x,dataV3x[:,300] ,label = 'chap. 3')
plt.grid()
plt.title('vitess outflow')
plt.xlabel(labels[1])
plt.legend(loc="best")
plt.tight_layout()

#.. fig 3
#
#plt.figure()
#plt.plot(Y,dataN1[200,:] ,label = '$n_e$')
#plt.plot(Y,np.abs(dataB1[200,:] ),label = '$\\vert B \\vert$')
#plt.xlim([-40.0,40.0])
#plt.grid()
#plt.ylabel('$norm.$')
#plt.xlabel(labels[1])
#plt.legend(loc="best")
#plt.tight_layout()
#

plt.show()

normal = {'b0': 1.0e-8, 'n0': 2} # for hedp : b0 in megaGauss & n0 in cm-3
#
# .. load the data for image
#im0 = field.Field(run = run,
#                  data = data,
#                  domain = domain,
#                  shifts = shifts)
##
## .. load the data for contours
#im1 = field.Field(run = run,
#                  data = dutu,
#                  domain = domain,
#                  shifts = shifts)
#
##vectim0 = flatten(abs(im0.data))
#maxim0= np.max(np.max(np.abs(im0.data)))
##bounds=[-maxim0,maxim0]
## .. draw the plot
#plo = colp.Colp(coloraxis = im0.axis,
#                colordata = im0.data,
#                bounds = bounds,
#                colormap = colormap,
#                contouraxis = im1.axis,
#                contourdata = im1.data,
#                flines = flines,
#                arrowaxis = None,
#                arrowdata = None,
#                labels = labels,
#                ticks = ticks,
#                subticks = subticks,
#                figsize = figsize,
#                filetype = filetype,
#                filename = filename,
#                vlines = None)
#
