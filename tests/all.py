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
#
# .. load the run
#path = '/media/sbolanos/DATA/init_hkl/'
#path = '/media/sbolanos/BatDRIVE/HECKLE/TEST/store/las_CU_L18_gf_n10/'
path = '/media/sbolanos/BatDRIVE/HECKLE/laser/las_Al_plan/'
#path = '/media/sbolanos/BatDRIVE/HECKLE/400x600_te/T10/400x600_vA_T10_t0-18/'
#path = '/media/sbolanos/BatDRIVE/HECKLE/PRO/T10/400x600_T10/'
#path = '/media/sbolanos/BatDRIVE/HECKLE/TEST/store/CU_D48_L16_r1/'
#path = '/media/sbolanos/BatDRIVE/HECKLE/PRO/T10/400x600_T10/'
name =  'py3'
run  = heckle.Heckle(path, name)
#
# .. get the desired data giving time
time = 0.0

dataN = run.GetN(time,"p")
dataNe= run.GetN(time,"e")
dataNb= run.GetN(time,"b")
# B : magnetic
dataBx = run.GetB(time)[...,0]
dataBy = run.GetB(time)[...,1]
dataBz = run.GetB(time)[...,2]
# E : electric
dataEx = run.GetE(time)[...,0]
dataEy = run.GetE(time)[...,1]
dataEz = run.GetE(time)[...,2]
# J : Current 
dataJx = run.GetJ(time)[...,0]
dataJy = run.GetJ(time)[...,1]
dataJz = run.GetJ(time)[...,2]
# V : velocity 
dataVx_e = run.GetV(time,'e')[...,0]
dataVy_e = run.GetV(time,'e')[...,1]
dataVz_e = run.GetV(time,'e')[...,2]
# V : velocity 
dataVx_i = run.GetV(time,'p')[...,0]
dataVy_i = run.GetV(time,'p')[...,1]
dataVz_i = run.GetV(time,'p')[...,2]
# V : velocity 
dataVx_a = run.GetV(time,"alpha")[...,0]
dataVy_a = run.GetV(time,"alpha")[...,1]
dataVz_a = run.GetV(time,"alpha")[...,2]

dutu = run.fourierFlux(time)

#
# .. set the needed parameter to plot the 2d field
domain    = None
shifts    = [0, 0]
bounds    = None
colormap  = ['bwr', 64]
colordens = ['jet', 64]
flines    = - np.arange(0.,150.,2.0)[::-1]
ticks     = [20, 20]
subticks  = [2, 2]
figsize   = [9,6]
filetype  = 'png'
filename   = name+'_t_'+str(time)
filenameN  = 'density_'+name+'_t_'+str(time)
filenameNe = 'densityelec_'+name+'_t_'+str(time)
filenameNb = 'densityback_'+name+'_t_'+str(time)
filenameBx = 'bx_'+name+'_t_'+str(time)
filenameBy = 'by_'+name+'_t_'+str(time)
filenameBz = 'bz_'+name+'_t_'+str(time)
filenameEx = 'ex_'+name+'_t_'+str(time)
filenameEy = 'ey_'+name+'_t_'+str(time)
filenameEz = 'ez_'+name+'_t_'+str(time)
filenameJx = 'jx_'+name+'_t_'+str(time)
filenameJy = 'jy_'+name+'_t_'+str(time)
filenameJz = 'jz_'+name+'_t_'+str(time)
filenameVx_e = 'vex_'+name+'_t_'+str(time)
filenameVy_e = 'vey_'+name+'_t_'+str(time)
filenameVz_e = 'vez_'+name+'_t_'+str(time)
filenameVx_i = 'vix_'+name+'_t_'+str(time)
filenameVy_i = 'viy_'+name+'_t_'+str(time)
filenameVz_i = 'viz_'+name+'_t_'+str(time)
#
normal = {'b0': 1.0e-8, 'n0': 2} # for hedp : b0 in megaGauss & n0 in cm-3
#
if filetype is not None:
  os.mkdir(path+filename)
  os.chdir(path+filename)

# .. load the data for image
imN = field.Field(run = run,
                  data = dataN,
                  domain = domain,
                  shifts = shifts)

#
imNe = field.Field(run = run,
                  data = dataNe,
                  domain = domain,
                  shifts = shifts)

#
imNb = field.Field(run = run,
                  data = dataNb,
                  domain = domain,
                  shifts = shifts)

#
imBx = field.Field(run = run,
                  data = dataBx,
                  domain = domain,
                  shifts = shifts)

maxBx= np.max(np.max(np.abs(imBx.data)))
boundsBx = [-maxBx,maxBx]
#
imBy = field.Field(run = run,
                  data = dataBy,
                  domain = domain,
                  shifts = shifts)

maxBy= np.max(np.max(np.abs(imBy.data)))
boundsBy = [-maxBy,maxBy]
#
imBz = field.Field(run = run,
                  data = dataBz,
                  domain = domain,
                  shifts = shifts)

maxBz= np.max(np.max(np.abs(imBz.data)))
boundsBz = [-maxBz,maxBz]
#
imEx = field.Field(run = run,
                  data = dataEx,
                  domain = domain,
                  shifts = shifts)

maxEx= np.max(np.max(np.abs(imEx.data)))
boundsEx = [-maxEx,maxEx]
#
imEy = field.Field(run = run,
                  data = dataEy,
                  domain = domain,
                  shifts = shifts)

maxEy= np.max(np.max(np.abs(imEy.data)))
boundsEy = [-maxEy,maxEy]
#
imEz = field.Field(run = run,
                  data = dataEz,
                  domain = domain,
                  shifts = shifts)

maxEz= np.max(np.max(np.abs(imEz.data)))
boundsEz = [-maxEz,maxEz]
#
imJx = field.Field(run = run,
                  data = dataJx,
                  domain = domain,
                  shifts = shifts)

maxJx= np.max(np.max(np.abs(imJx.data)))
boundsJx = [-maxJx,maxJx]
#
imJy = field.Field(run = run,
                  data = dataJy,
                  domain = domain,
                  shifts = shifts)

maxJy= np.max(np.max(np.abs(imJy.data)))
boundsJy = [-maxJy,maxJy]
#
imJz = field.Field(run = run,
                  data = dataJz,
                  domain = domain,
                  shifts = shifts)

maxJz= np.max(np.max(np.abs(imJz.data)))
boundsJz = [-maxJz,maxJz]
#
imVx_e = field.Field(run = run,
                  data = dataVx_e,
                  domain = domain,
                  shifts = shifts)

maxVx_e= np.max(np.max(np.abs(imVx_e.data)))
boundsVx_e = [-maxVx_e,maxVx_e]
#
imVy_e = field.Field(run = run,
                  data = dataVy_e,
                  domain = domain,
                  shifts = shifts)

maxVy_e= np.max(np.max(np.abs(imVy_e.data)))
boundsVy_e = [-maxVy_e,maxVy_e]
#
imVz_e = field.Field(run = run,
                  data = dataVz_e,
                  domain = domain,
                  shifts = shifts)

maxVz_e= np.max(np.max(np.abs(imVz_e.data)))
boundsVz_e = [-maxVz_e,maxVz_e]
#
imVx_i = field.Field(run = run,
                  data = dataVx_i,
                  domain = domain,
                  shifts = shifts)

maxVx_i= np.max(np.max(np.abs(imVx_i.data)))
boundsVx_i = [-maxVx_i,maxVx_i]
#
imVy_i = field.Field(run = run,
                  data = dataVy_i,
                  domain = domain,
                  shifts = shifts)

maxVy_i= np.max(np.max(np.abs(imVy_i.data)))
boundsVy_i = [-maxVy_i,maxVy_i]
#
imVz_i = field.Field(run = run,
                  data = dataVz_i,
                  domain = domain,
                  shifts = shifts)

maxVz_i= np.max(np.max(np.abs(imVz_i.data)))
boundsVz_i = [-maxVz_i,maxVz_i]

#
#
# .. load the data for contours
im1 = field.Field(run = run,
                  data = dutu,
                  domain = domain,
                  shifts = shifts)

# .. draw the plot
plo = colp.Colp(coloraxis = imN.axis,
                colordata = imN.data,
                bounds = bounds,
                colormap = colordens,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameN)

ploNe = colp.Colp(coloraxis = imNe.axis,
                colordata = imNe.data,
                bounds = bounds,
                colormap = colordens,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameNe)

ploNb = colp.Colp(coloraxis = imNb.axis,
                colordata = imNb.data,
                bounds = bounds,
                colormap = colordens,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameNb)

ploBx = colp.Colp(coloraxis = imBx.axis,
                colordata = imBx.data,
                bounds = boundsBx,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameBx)

ploBy = colp.Colp(coloraxis = imBy.axis,
                colordata = imBy.data,
                bounds = boundsBy,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameBy)

ploBz = colp.Colp(coloraxis = imBz.axis,
                colordata = imBz.data,
                bounds = boundsBz,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameBz)

ploEx = colp.Colp(coloraxis = imEx.axis,
                colordata = imEx.data,
                bounds = boundsEx,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameEx)

ploEy = colp.Colp(coloraxis = imEy.axis,
                colordata = imEy.data,
                bounds = boundsEy,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameEy)

ploEz = colp.Colp(coloraxis = imEz.axis,
                colordata = imEz.data,
                bounds = boundsEz,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameEz)

ploJx = colp.Colp(coloraxis = imJx.axis,
                colordata = imJx.data,
                bounds = boundsJx,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameJx)

ploJy = colp.Colp(coloraxis = imJy.axis,
                colordata = imJy.data,
                bounds = boundsJy,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameJy)

ploJz = colp.Colp(coloraxis = imJz.axis,
                colordata = imJz.data,
                bounds = boundsJz,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameJz)
ploVx_e = colp.Colp(coloraxis = imVx_e.axis,
                colordata = imVx_e.data,
                bounds = boundsVx_e,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameVx_e)

ploVy_e = colp.Colp(coloraxis = imVy_e.axis,
                colordata = imVy_e.data,
                bounds = boundsVy_e,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameVy_e)

ploVz_e = colp.Colp(coloraxis = imVz_e.axis,
                colordata = imVz_e.data,
                bounds = boundsVz_e,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameVz_e)
#
ploVx_i = colp.Colp(coloraxis = imVx_i.axis,
                colordata = imVx_i.data,
                bounds = boundsVx_i,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameVx_i)

ploVy_i = colp.Colp(coloraxis = imVy_i.axis,
                colordata = imVy_i.data,
                bounds = boundsVy_i,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameVy_i)

ploVz_i = colp.Colp(coloraxis = imVz_i.axis,
                colordata = imVz_i.data,
                bounds = boundsVz_i,
                colormap = colormap,
                contouraxis = im1.axis,
                contourdata = im1.data,
                flines = flines,
                arrowaxis = None,
                arrowdata = None,
                labels = imN.labels,
                ticks = ticks,
                subticks = subticks,
                figsize = figsize,
                filetype = filetype,
                filename = filenameVz_i)
##
