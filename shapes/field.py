
import numpy as np
from . import shape


class Field(shape.Shape):

   """

   """

   def __init__(self,
                run,
                data,
                domain = None,
                shifts = None,
                stride = None,
                section = None,
                labels = None):

      # set the "limit" list
      limit = []
      for i in range(run.ndim) :
         limit.append([0, run.domsize[i]])

      super(Field, self).__init__(data, limit, domain, stride, section)

      # ..... private members
      self._run = run

      # ..... public members

      # .. set the shifts for the axis
      self._setShifts(shifts)

      # .. set labels
      self._setLabels(section, labels)


   def _setShifts(self,
                  shifts) :


      if shifts == None :
         shifts = self.axis.__len__()*[0]

      else :
         for i, val in enumerate(shifts) :
            if val == None :
               shifts[i] = 0.0
            else :
               pass

      for i in range(self._rank) :
         self.axis[i] += shifts[i]


   def _setLabels(self,
                  section,
                  labels) :

      """
      should be virtual & defined in field, fourier... !!!!!!!!!!!!!!!!!!!!!!!!!
      """

      if labels != None :
         self.labels = labels

      else :
         self.labels = []
         s = ['x', 'y', 'z']

         if section == None :
            for i in range(self.data.ndim) :
               self.labels.append('$'+s[i]+' c / \omega_p$')

         else :
            for i in range(section.__len__()) :
               if section[i] == None :
                  self.labels.append('$'+s[i]+' c / \omega_p$')

               else :
                  pass

