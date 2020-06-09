#
#
#
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc
#
#
class Linp :

    """

    """

    def __init__(self,
                 axis,
                 data,
                 bounds = None,
                 step = None,
                 labels = None,
                 ticks = None,
                 subticks = None,
                 linestyle = None,
                 linecolor = None,
                 linewidth = None,
                 markers = None,
                 legend = None,
                 figsize = [5, 4],
                 filetype = None):

        # .. set the x data
        self.axis = axis

        # .. set the y data
        self.data = data

        self.labels = labels
        # .. set the position of shaded rectangle(s)
        #self.shade = shade

        # .. set the curves needed for diff
        #self.diff = diff

        # .. set the position of ticks
        self.ticks = ticks

        # .. set the # of subticks
        self.subticks = subticks

        # .. set the linestyle
        if linestyle != None :
           self.linestyle = linestyle

        else :
            self.linestyle = self.axis.__len__()*["-"]

        # .. set the linecolor
        if linecolor != None :
           self.linecolor = linecolor

        else :
            self.linecolor = self.axis.__len__()*["0.0"]

        # .. set the linewidth
        if linewidth != None :
           self.linewidth = linewidth

        else :
            self.linewidth = self.axis.__len__()*[1.0]

        # .. set the markers
        if markers != None :
           self.markers = markers

        else :
            self.markers = self.axis.__len__()*[None]


        # .. set the drawstyle
        if step != None :
           self.step = step

        else :
            self.step = self.axis.__len__()*['default']
        print self.step


        # .. set the legend
        self.legend = legend

        # .. set the figsize
        self.figsize = figsize

        # .. which filetype
        if filetype is not None :
           self.filetype = filetype.lower()
        else :
           self.filetype = None

        # .. draw the plot...
        self.Draw(bounds)



    def Draw(self,
             bounds):

        """
        """

        # .. use latex fonts
        rc('text', usetex = True)

        # .. just 1 plot
        fig, ax = plt.subplots(num = 0, figsize = self.figsize, dpi = 100)

        # .. plot the data
        #if self.axis.__len__() == 1 :
        #   ax.plot(self.axis[0],
        #           self.data,
        #           linestyle = self.linestyle,
        #           color = self.linecolor,
        #           linewidth = self.linewidth)
        # drawstyle = "steps-mid"
        #else :
        npl = self.axis.__len__()
        for i in range(npl) :
           ax.plot(self.axis[i],
                   self.data[i],
                   drawstyle = self.step[i],
                   linestyle = self.linestyle[i],
                   color = self.linecolor[i],
                   linewidth = self.linewidth[i],
                   marker = self.markers[i])

        # fill between 2 curves
        #if self.diff != None :
        #  ax.fill_between(self.axis[0], self.data[self.diff[0]], self.data[self.diff[1]], where = self.data[self.diff[1]] >= self.data[self.diff[0]], facecolor = '0.4', interpolate = True)
        #  ax.fill_between(self.axis[0], self.data[self.diff[0]], self.data[self.diff[1]], where = self.data[self.diff[1]] <= self.data[self.diff[0]], facecolor = '0.8', interpolate = True)

        # .. draw shaded rectangle(s)
        #if self.shade != None :
        #  for i in range(self.shade.__len__()) : plt.axvspan(self.shade[i][0], self.shade[i][1], facecolor = "0.8", alpha = 1.0, linewidth = 0.0)

        # .. set the ticks
        if self.ticks != None :
           ax.xaxis.set_major_locator(MultipleLocator(self.ticks[0]))
           ax.yaxis.set_major_locator(MultipleLocator(self.ticks[1]))

        # .. set the subticks
        if self.subticks != None :
           ax.xaxis.set_minor_locator(MultipleLocator(self.subticks[0]))
           ax.yaxis.set_minor_locator(MultipleLocator(self.subticks[1]))

        # .. set the limits of the plot
        ax.set_xlim([np.min(self.axis[:][0]), np.max(self.axis[:][-1])])
        ax.set_ylim(bounds)

        # .. set a legend if wanted
        if self.legend != None :
           ax.legend(self.legend)

        # .. write the labels
        ax.set_xlabel(self.labels[0])
        ax.set_ylabel(self.labels[1])

        # .. draw the plot
        if self.filetype is not None :
           fig.savefig('linp.'+self.filetype.lower(),
           dpi = 200,
           bbox_inches = 'tight')
        else :
           plt.show()

        # .. clear the figure
        fig.clf()
