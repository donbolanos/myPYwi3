#
#! /usr/bin/env python
#
#
#import matplotlib as mpl
#mpl.use('Agg')
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import meshgrid
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc
import operator
#
#
class i1a:

    """

    """

    def __init__(self,
                 handle = None,
                 colordata = None,
                 contourdata = None,
                 flines = 10,
                 ticks = None,
                 subticks = None,
                 strticks = None,
                 title = None,
                 cmap = None,
                 figsize = [4.65, 3.51],
                 filetype = None):

        # ..
        self.handle = handle

        # ..
        self.colordata = colordata

        # ..
        self.contourdata = contourdata

        # ..
        self.flines = flines

        # .. set the position of ticks
        self.ticks = ticks

        # .. set the # of subticks
        self.subticks = subticks

        # .. set strticks
        self.strticks = strticks

        # .. set the title
        self.title = title

        # ..
        if cmap == None: self.cmap = cm.get_cmap('jet', 8)
        else : self.cmap = cmap

        # ..
        self.figsize = figsize

        # ..
        if filetype is not None : self.filetype = filetype.lower()
        else : self.filetype = None

        # ..
        self.Draw()


    def Draw(self):

        """
        draw the plot
        """

        # .. use latex fonts
        rc('text', usetex = True)
        #rc('text.latex', preamble = r'\usepackage{cmbright}')

        # .. just 1 plot
        fig, ax = plt.subplots(num = 0, figsize = self.figsize, dpi = 100)

        # .. build the mesh
        xp, yp = meshgrid(self.handle.xdata, self.handle.ydata)

        # .. draw the color image
        if self.colordata is not None:
          im = ax.imshow(np.transpose(self.colordata),
                         aspect = 'auto',
                         interpolation = 'nearest',
                         cmap = self.cmap,
                         origin = 'lower',
                         extent = reduce(operator.add, self.handle.extent),
                         vmin = self.handle.bounds[0],
                         vmax = self.handle.bounds[1])

        # .. draw the isocontour
        if self.contourdata is not None:
          co = ax.contour(xp, yp, np.transpose(self.contourdata), self.flines,
                          colors = ('k',),
                          origin = 'lower',
                          extent = reduce(operator.add, self.handle.extent),
                          linestyles =  'solid',
                          linewidths = 1)

        # .. set the x ticks & subticks
        if self.ticks != None : ax.xaxis.set_major_locator(MultipleLocator(self.ticks[0]))
        if self.strticks != None : ax.xaxis.set_major_formatter(FormatStrFormatter(self.strticks[0]))
        if self.subticks != None : ax.xaxis.set_minor_locator(MultipleLocator(self.subticks[0]))

        # .. set the x ticks & subticks
        if self.ticks != None : ax.yaxis.set_major_locator(MultipleLocator(self.ticks[1]))
        if self.strticks != None : ax.yaxis.set_major_formatter(FormatStrFormatter(self.strticks[1]))
        if self.subticks != None : ax.yaxis.set_minor_locator(MultipleLocator(self.subticks[1]))

        # .. set the limits of the plot
        ax.set_xlim([self.handle.extent[0][0], self.handle.extent[0][1]])
        ax.set_ylim([self.handle.extent[1][0], self.handle.extent[1][1]])

        # .. write the title
        if self.title is not None : ax.text(0.98*self.handle.extent[0][0]+0.02*self.handle.extent[0][1], -0.04*self.handle.extent[1][0]+1.04*self.handle.extent[1][1], self.title)

        # .. write the labels
        ax.set_xlabel(self.handle.xlabel)
        ax.set_ylabel(self.handle.ylabel)

        # .. 3 ticks for the colorbar
        ticks = np.linspace(self.handle.bounds[0], self.handle.bounds[1], num = 3)

        # .. set the colorbar
        cbar = fig.colorbar(im, ticks = ticks, pad = 0.03, aspect = 40)

        # .. draw the plot
        if self.filetype is not None : fig.savefig(self.handle.filename + '.'+self.filetype.lower(), dpi = 200, bbox_inches = 'tight')
        else : plt.show()

        # .. clear the figure
        fig.clf()

