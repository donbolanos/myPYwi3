#
#! /usr/bin/env python
#
#
#import matplotlib as mpl
#mpl.use('Agg')
#
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc
#
#
class l1a :

    """

    """

    def __init__(self,
                 handle = None,
                 axis,
                 data,
                 shade = None,
                 diff = None,
                 step = None,
                 linestyle = None,
                 linecolor = None,
                 linewidth = None,
                 ticks = None,
                 subticks = None,
                 strticks = None,
                 marker = None,
                 legend = None,
                 title = None,
                 figsize = [5, 4],
                 filetype = None):

        # .. set the handle
        self.handle = handle

        # .. set the x data
        self.axis = axis

        # .. set the y data
        self.data = data

        # .. set the position of shaded rectangle(s)
        self.shade = shade

        # .. set the curved needed for diff
        self.diff = diff

        # .. set the position of ticks
        self.ticks = ticks

        # .. set the # of subticks
        self.subticks = subticks

        # .. set strticks
        self.strticks = strticks

        # .. set the linestyle
        if linestyle != None : self.linestyle = linestyle
        else :
          if self.axis.ndim == 1 :
            self.linestyle = "-"
          else :
            npl = self.axis.shape[0]
            self.linestyle = npl*["-"]

        # .. set the linecolor
        if linecolor != None : self.linecolor = linecolor
        else :
          if self.axis.ndim == 1 :
            self.linecolor = "0.0"
          else :
            npl = self.axis.shape[0]
            self.linecolor = npl*["0.0"]

        # .. set the linewidth
        if linewidth != None : self.linewidth = linewidth
        else :
          if self.axis.ndim == 1 :
            self.linewidth = 1.0
          else :
            npl = self.axis.shape[0]
            self.linewidth = npl*[1.0]

        # .. set the linewidth
        if marker != None : self.marker = marker
        else :
          if self.axis.ndim == 1 :
            self.marker = [None]
          else :
            npl = self.axis.shape[0]
            self.marker = npl*[None]

        # .. set the legend
        self.legend = legend

        # .. set the title
        self.title = title

        # .. set the figsize
        self.figsize = figsize

        # .. which filetype
        if filetype is not None : self.filetype = filetype.lower()
        else : self.filetype = None

        # .. draw the plot...
        self.Draw()


    def Draw(self):

        """
        """

        # .. use latex fonts
        rc('text', usetex = True)

        # .. just 1 plot
        fig, ax = plt.subplots(num = 0, figsize = self.figsize, dpi = 100)

        # .. plot the data
        if self.axis.ndim == 1 :
          ax.plot(self.axis, self.data, linestyle = self.linestyle, color = self.linecolor, linewidth = self.linewidth)
          # drawstyle = "steps-mid"
        else :
          npl = self.axis.shape[0]
          for i in range(npl) :
            ax.plot(self.axis[i], self.data[i], linestyle = self.linestyle[i], color = self.linecolor[i], marker = self.marker[i], linewidth = self.linewidth[i])

        # fill between 2 curves
        if self.diff != None :
          ax.fill_between(self.axis[0], self.data[self.diff[0]], self.data[self.diff[1]], where = self.data[self.diff[1]] >= self.data[self.diff[0]], facecolor = '0.4', interpolate = True)
          ax.fill_between(self.axis[0], self.data[self.diff[0]], self.data[self.diff[1]], where = self.data[self.diff[1]] <= self.data[self.diff[0]], facecolor = '0.8', interpolate = True)

        # .. draw shaded rectangle(s)
        if self.shade != None :
          for i in range(self.shade.__len__()) : plt.axvspan(self.shade[i][0], self.shade[i][1], facecolor = "0.8", alpha = 1.0, linewidth = 0.0)

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

        # .. set a legend if wanted
        if self.legend != None : ax.legend(self.legend)

        # .. write the title
        if self.title is not None : ax.text(0.98*self.handle.extent[0][0]+0.02*self.handle.extent[0][1], -0.04*self.handle.extent[1][0]+1.04*self.handle.extent[1][1], self.title)

        # .. write the labels
        ax.set_xlabel(self.handle.xlabel)
        ax.set_ylabel(self.handle.ylabel)

        # .. draw the plot
        if self.filetype is not None : fig.savefig(self.handle.filename + '.'+self.filetype.lower(), dpi=200, bbox_inches = 'tight')
        else : plt.show()

        # .. clear the figure
        fig.clf()
