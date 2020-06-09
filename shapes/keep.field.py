
import numpy as np


class Shape():

   """

   """

   def __init__(self,
                run,
                data,
                domain = None,
                shift = None,
                plane = None,
                direction = None,
                cut = None,
                stride = None,
                label = None):

      # ..... private members
      self._run = run
      self._domain = domain
      self._cut = cut

      # ..... public members
      self.xdata = None
      self.ydata = None
      self.zdata = None

      # .. for the 1d case, direction & cut have to be defined
      self._setDirection(direction)

      # .. for the 2d case, plane & cut have to be defined
      self._setPlane(plane)

      # .. for the 3d case, no slice needed !

      # .. set the shift for the axis
      self._setShift(shift)

      # .. set extent : cell index in a list
      self._setExtent()

      # .. set domain : in physical units, including shift
      self._setDomain()

      # .. get a 1d field with the right shape
      if self.direction is not None :
         self._shape1D(data)
         self._make1DLabel(label)

      # .. get a 2d field with the right shape
      if self.plane is not None :
         self._shape2D(data)
         self._make2DLabel(label)

      # .. get a 3d field with the right shape :


   def _setPlane(self,
                 plane) :

      if plane is not None :

         if self._run.ndim == 1 :
            raise NotImplementedError("what the fuck you want a 2D plot from a 1D run !")

         elif self._run.ndim == 2 :
            self.plane = "xy"
            self._cut = None

         elif self._run.ndim == 3 :
            self.plane = plane.lower()
            if (self.plane in ["xy", "yx"]) :
               self._cut = np.rint(cut/self._run.dl[2])
            elif (self.plane in ["xz", "zx"]) :
               self._cut = np.rint(cut/self._run.dl[1])
            elif (self.plane in ["yz", "zy"]) :
               self._cut = np.rint(cut/self._run.dl[0])
            else :
               raise ValueError( "what the fuck is this plane ?" )

         else :
            raise ValueError( "what the fuck is the dim of your run ?" )

      else :
         self.plane = None


   def _setDirection(self,
                     direction) :

      if direction is not None :

         if self._run.ndim == 1 :
            self.direction = "x"

         elif self._run.ndim == 2 :

            self.direction = direction.lower()

            if self.direction == "x" :
               self._cut = np.rint(cut/self._run.dl[1])

            elif self.direction == "y" :
               self._cut = np.rint(cut/self._run.dl[0])

            else :
               raise ValueError( "what the fuck is the direction you want to deal with ?" )

         elif self._run.ndim == 3 :

            self.direction = direction.lower()

            if self.direction == "x" :
               self._cut = [np.rint(cut[0]/self._run.dl[1]),
                            np.rint(cut[1]/self._run.dl[2])]

            elif self.direction == "y" :
               self._cut = [np.rint(cut[0]/self._run.dl[0]),
                            np.rint(cut[1]/self._run.dl[2])]

            elif self.direction == "z" :
               self._cut = [np.rint(cut[0]/self._run.dl[0]),
                            np.rint(cut[1]/self._run.dl[1])]
            else :
               raise ValueError( "what the fuck is the direction you want to deal with ?" )

         else :
            raise ValueError( "what the fuck is the dim of your run ?" )

      else :
         self.direction = None


   def _setShift(self,
                 shift) :

      if shift is None :

         if self._run.ndim == 1 :
            self._shift = [0]

         elif self._run.ndim == 2 :
            self._shift = [0, 0]

         elif rself._un.ndim == 3 :
            self._shift = [0, 0, 0]

      else :
         self._shift = shift


   def _setExtent(self) :

      if self._domain == None :

         if self._run.ndim == 1 :
            self._extent = [[0, self._run.ncells[0]]]

         elif self._run.ndim == 2 :
            self._extent = [[0, self._run.ncells[0]],\
                            [0, self._run.ncells[1]]]

         elif self._run.ndim == 3 :
            self._extent = [[0, self._run.ncells[0]],\
                            [0, self._run.ncells[1]],\
                            [0, self._run.ncells[2]]]

         else :
            raise ValueError( "what the fuck is the dim of your run ?" )

      else :

         if self._run.ndim == 1 :
            xlims = self._run.coord2Index(self._domain[0], axis = 0)
            self._extent = [xlims]

         elif self._run.ndim == 2 :
            xlims = self._run.coord2Index(self._domain[0], axis = 0)
            ylims = self._run.coord2Index(self._domain[1], axis = 1)
            self._extent = [xlims, ylims]

         elif self._run.ndim == 3 :
            xlims = self._run.coord2Index(self._domain[0], axis = 0)
            ylims = self._run.coord2Index(self._domain[1], axis = 1)
            zlims = self._run.coord2Index(self._domain[2], axis = 2)
            self._extent = [xlims, ylims, zlims]

         else :
            raise ValueError( "what the fuck is the dim of your run ?" )


   def _setDomain(self) :

      if self._run.ndim == 1 :
         xlims = self._run.index2Coord(self._extent[0], axis = 0)
         self._domain = [xlims+self._shift[0]]

      elif self._run.ndim == 2 :
         xlims = self._run.index2Coord(self._extent[0], axis = 0)
         ylims = self._run.index2Coord(self._extent[1], axis = 1)
         self._domain = [xlims+self._shift[0],\
                         ylims+self._shift[1]]

      elif self._run.ndim == 3 :
         xlims = self._run.index2Coord(self._extent[0], axis = 0)
         ylims = self._run.index2Coord(self._extent[1], axis = 1)
         zlims = self._run.index2Coord(self._extent[2], axis = 2)
         self._domain = [xlims+self._shift[0],\
                         ylims+self._shift[1],\
                         zlims+self._shift[2]]

      else :
         raise ValueError( "what the fuck is the dim of your run ?" )


   def _shape2D(self, rawdata):

      """
      input : rawdata is a numpy array of rank 2 (for a 2d run)
              or 3 (for a 3d run)

      set value : xdata is a 1d array containing the x axis values (size nx)
                  ydata is a 1d array containing the y axis values (size ny)
                  zdata is a 2d array with the rawdata values (size nx by ny)
      """

      self.xdata = np.linspace(self._domain[0][0],\
                               self._domain[0][1],\
                               self._extent[0][1]-self._extent[0][0]+1)
      self.ydata = np.linspace(self._domain[1][0],\
                               self._domain[1][1],\
                               self._extent[1][1]-self._extent[1][0]+1)

      if self._run.ndim == 1 :
         raise ValueError( "what the fuck you want to reshape data a 1d run ?" )

      elif self._run.ndim == 2 :
           self.zdata = rawdata

      elif self._run.ndim == 3 :
         if self._plane == "xy" :
            self.zdata = rawdata[:, :, self._cut]

         elif self._plane == "xz" :
            self.zdata = rawdata[:, self._cut, :]

         elif self._plane == "yz" :
            self.zdata = rawdata[self._cut, :, :]

         else :
            raise ValueError( "what the fuck is this plane ?" )

      else :
         raise ValueError( "what the fuck is the dim of your run ?" )


   def _shape1D(self, rawdata) :

      """
      """

      self.xdata = np.linspace(self._domain[0][0],\
                               self._domain[0][1],\
                               self._extent[0][1]-self._extent[0][0]+1)

      if self._run.ndim == 1 :
         self.ydata = rawdata

      elif self._run.ndim == 2 :
         if self.direction == "x" :
            self.ydata = rawdata[:, self._cut[0]]

         elif self.direction == "y" :
            self.ydata = rawdata[self._cut[0], :]

         else :
            raise ValueError( "what the fuck is this direction ?" )

      elif self._run.ndim == 3 :
         if self.direction == "x" :
            self.ydata = rawdata[:, self._cut[0], self._cut[1]]

         elif self.direction == "y" :
            self.ydata = rawdata[self._cut[0], :, self._cut[1]]

         elif self.direction == "z" :
            self.ydata = rawdata[self._cut[0], self._cut[1], :]

         else :
            raise ValueError( "what the fuck is this direction ?" )


   def _make2DLabel(self, label) :

      """
      should be virtual & defined in field, fourier... !!!!!!!!!!!!!!!!!!!!!!!!!
      """

      if self._plane == "xy" :
         xstr = 'x'
         ystr = 'y'

      if self._plane == "xz" :
         xstr = 'x'
         ystr = 'z'

      if self._plane == "yz" :
         xstr = 'y'
         ystr = 'z'

      if label is None or label[0] is None :
         xlabel = r'$'+xtrs+'\,\Omega_p/c$'

      else :
         xlabel = label[0]

      if label is None or label[1] is None :
         xlabel = r'$'+ytrs+'\,\Omega_p/c$'

      else :
         ylabel = label[1]

      self.label = [xlabel, ylabel]


   def _make1DLabel(self, label) :

      """
      should be virtual & defined in field, fourier... !!!!!!!!!!!!!!!!!!!!!!!!!
      """

      if label is None or label[0] is None :
         if self.direction == "x" :
             xlabel = r'$x \, \Omega_p/c$'

         if self.direction == "y" :
             xlabel = r'$y \, \Omega_p/c$'

         if self.direction == "z" :
             xlabel = r'$z \, \Omega_p/c$'

      else :
         xlabel = label[0]

      if label is None or label[1] is None :
          ylabel = r'$ $'

      else :
         ylabel = label[1]

      self.label = [xlabel, ylabel]


