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
flatten = lambda l: [item for sublist in l for item in sublist]
# .. load the run
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/45deg_T10/laser_046r/'
#path = '/media/sbolanos/TOURO/HECKLE/270x540_ths/t_8-38/'
#path = '/media/sbolanos/BatDRIVE/HECKLE/ANDREE/102/'
path =  '/media/sbolanos/BatDRIVE/HECKLE/TEST/store/CU_L16/'
#path = '/media/sbolanos/BatDRIVE/HECKLE/PRO/anti60_T10/400x600_pro60_T10/'
name = ''
run  = heckle.Heckle(path, name)
#
# .. get the desired data giving time
time = 60.0
#data = run.GetN(time, "e")
#data1 = run.GetB(time)[...,2]
#data = np.sqrt(run.GetB(time)[...,0]**2 + run.GetB(time)[...,1]**2 + run.GetB(time)[...,2]**2)
#dataM = np.sqrt(data*data +data1*data1+data2*data2)
dutu = run.fourierFlux(time)
#data = run.GetJ(time)[...,2]
data = run.GetN(time,'p')

#
# .. set the needed parameter to plot the 2d field
domain    = None
shifts    = [-40, -60]
bounds    = None#[0.0,1.0]#[-0.4,0.4]
colormap  = ['viridis', 64]#['bwr', 64]
flines    = 8
ticks     = None#[60, 30]
labels    = ['$zc / \omega_{pi}$','$yc / \omega_{pi}$']
subticks  = [5,5]
figsize   = [6,7]
filetype  = None#'pdf'
filename  = name+'_t_'+str(time)
vlines    = None #[30*0.2]
#
zoomx = [100,300]
zoomy = [250,350]

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
                vlines = vlines)
#
