"""
.. module:: heckle
    :synopsis: class for data produced by the code Heckle
.. moduleauthor:: Nicolas Aunai <nicolas.aunai@lpp.polytechnique.fr>
"""




from . import run
import h5py
import os
import numpy as np

if __name__ == '__main__':
    main()


class Heckle(run.Run):

    """
    test
    """

    #----------------------------------------------------------
    #----------------------------------------------------------
    def __init__(self, path, runID = 'r', fieldname = 'fields.h5'):

        """
        Exemple  :

        Creation : 2015-02-27 16:49:57.016701

        """
        super(Heckle, self).__init__(path, runID=runID)

        # ----- heckle private members
        self._fieldfilename   = 'fields.h5'
        self._speciefilename  = 'species.h5'
        self._timefilename    = 'time.h5'
        self._restartfilename = 'restarts.h5'
        self._proton   = ['protons', 'proton', 'p']
        self._electron = ['electrons', 'electron', 'e']
        self._ion      = ['ions', 'ion', 'i']
        self._alpha    = ['alpha','background','b']
        # we often have a second proton population : could we call it "alphas" ?
        # (even if backgrounds are generally not alpha particles)

        # ----- heckle public members
        self.modelname  = 'Hybrid'
        self.codename   = 'Heckle'
        # Attribute of fields file
        self.ncells     = np.array(self._getFieldsAttrib('nbrOfCells'))
        self.domsize    = np.array(self._getFieldsAttrib('domainSize'))
        self.dl         = np.array(self._getFieldsAttrib('meshSize'))
        self.hprsty     = self._getFieldsAttrib('hyperResistivity')
        self.rsty       = self._getFieldsAttrib('resistivity')
        self.dumpfields = self._getFieldsAttrib('dumpFields')
        # Attribute of species file
        try : 
            self.charge     = np.array(self._getSpeciesAttrib('charge'))
            self.mass       = np.array(self._getSpeciesAttrib('mass'))
            self.weight     = np.array(self._getSpeciesAttrib('weight'))
            self.nparticles = np.array(self._getSpeciesAttrib('nbrOfParticles'))
            self.dumpspecies = self._getSpeciesAttrib('dumpSpecies')
        except : 
            pass
        # Attribute of species file
        # 1D case
        if self.ncells[1] == 1 :
            self.ncells  = np.array([self.ncells[0]])
            self.domsize = np.array([self.domsize[0]])
            self.dl      = np.array([self.dl[0]])

        # 2D case
        elif self.ncells[2] == 1 :
            self.ncells  = self.ncells[:-1]
            self.domsize = self.domsize[:-1]
            self.dl      = self.dl[:-1]

        self.ndim         = self.ncells.size
        self.sfield_shape = tuple(self.ncells+1)
        self.vfield_shape = tuple(self.ncells+1) + (3,)
        self.tfield_shape = tuple(self.ncells+1) + (3,3,)
        self.dtype        = np.float32 # what is the meaning ? only for fields.h5 ?
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def _getFieldsAttrib(self, attributeName):

        """
        Creation : 2015-02-27 19:12:14.062301
        """

        f = h5py.File(os.path.join(self.path, self._fieldfilename),'r')
        attribute = f.attrs[attributeName]
        f.close()

        return attribute

    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def _getSpeciesAttrib(self, attributeName):

        """
        Creation : 2018-05-28 
        """

        f = h5py.File(os.path.join(self.path, self._speciefilename),'r')
        attribute = f.attrs[attributeName]
        f.close()

        return attribute

    #==========================================================





    def getTimeGroups(self, time) :

        """
        get time (a single value or a list with format [tmin, tmax])
        and return a list of strings, containing all the time groupnames
        associated in "fields.h5"

        return a list with the times, and a list with the associated groups
        """

        f = h5py.File(os.path.join(self.path, self._fieldfilename),'r')

        # build the list of all time recorded in fields.h5
        groups = f.keys()
        timesfromfile = np.empty([0])

        # first parsed string is "time", second ":" & third the time we want
        for grp in groups :
           timesfromfile = np.append(timesfromfile , float(grp.strip().split()[2]))

        if time.__len__() == 1 :
            # good times, bad times, you know i've had my share
            goodtimes = timesfromfile[np.argmin(np.fabs(np.array(time)-timesfromfile))]

            timegroups = 'time : %f' % (goodtimes)

        elif time.__len__() == 2 :
            goodtimes = timesfromfile[(timesfromfile >= time[0]) & (timesfromfile <= time[1])]

            timegroups = ['time : %f' % (tim) for tim in goodtimes]

        else :
            raise ValueError( "time has to be a list of 1 or 2 float values")

        return goodtimes, timegroups

    def getTime(self, time, file) :

        """
        get time (a single value or a list with format [tmin, tmax])
        and return a list of strings, containing all the time groupnames
        associated in "file"

        return a list with the times, and a list with the associated groups
        """

        f = h5py.File(os.path.join(self.path, file),'r')

        # build the list of all time recorded in fields.h5
        groups = f.keys()
        timesfromfile = np.empty([0])

        # first parsed string is "time", second ":" & third the time we want
        for grp in groups :
           timesfromfile = np.append(timesfromfile , float(grp.strip().split()[2]))

        if time.__len__() == 1 :
            # good times, bad times, you know i've had my share
            goodtimes = timesfromfile[np.argmin(np.fabs(np.array(time)-timesfromfile))]

            timegroups = 'time : %f' % (goodtimes)

        elif time.__len__() == 2 :
            goodtimes = timesfromfile[(timesfromfile >= time[0]) & (timesfromfile <= time[1])]

            timegroups = ['time : %f' % (tim) for tim in goodtimes]

        else :
            raise ValueError( "time has to be a list of 1 or 2 float values")

        return goodtimes, timegroups


    #----------------------------------------------------------
    #----------------------------------------------------------
    def _readfieldsdataset(self, time, fieldname):

        """
        Creation : 2015-02-27 19:12:14.062301

           all time groups are under the "root dir" : "/"
           all dataset are in a given time group
        """

        mytimes, mygroups = self.getTimeGroups([time])

        f = h5py.File(os.path.join(self.path, self._fieldfilename),'r')

        data = f[mygroups+'/'+fieldname][()]

        f.close()

        # run 1D
        if self.ndim == 1 :
           mydata = data[:,0,0]

        # run 2D
        elif self.ndim == 2 :
           mydata = data[:,:,0]

        # run 3D
        elif self.ndim == 3 :
           mydata = data

        else :
            raise ValueError("what the fuck is this dim value : '%d' ?" % self.ndim)

        return mydata
    #==========================================================


    #----------------------------------------------------------
    #----------------------------------------------------------
    def _readspeciesdataset(self, time, species, fieldname):

        """
        Creation : 2015-02-27 19:12:14.062301

           all time groups are under the "root dir" : "/"
           all dataset are in a given time group
        """

        mytimes, mygroups = self.getTime([time],self._speciefilename)

        f = h5py.File(os.path.join(self.path, self._speciefilename),'r')

        data = f[mygroups+'/'+species+'/'+fieldname][()]

        f.close()
        
        mydata = data
        
        return mydata
    #==========================================================



    #----------------------------------------------------------
    #----------------------------------------------------------
    def _fillvec(self, v1, v2, v3):

        shape =  v1.shape + (3,)
        V     = np.zeros(shape = shape, dtype = v1.dtype)
        V[...,0] = v1
        V[...,1] = v2
        V[...,2] = v3

        return V

    # private method _fillvec could be changed by _fillall
    # this one is not restricted to vector : can fill tensor of whatever rank

    def _fillall(self, mylist) :

        MyListLen = mylist.__len__()
        shape = mylist[0].shape+(MyListLen,)
        Tab = np.zeros(shape = shape, dtype = mylist[0].dtype)

        for i in range(MyListLen) :
            Tab[...,i] = mylist[i]

        return Tab
    #==========================================================




    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetB(self, time, origin_grid = False):

        """
        Returns the magnetic field vector at desired time.

        Args:
            :param time: the time at which B is returned

        Kwargs:
            :param origin_grid (bool): [default:False] If True the\
                    magnetic components are defined on the Yee Grid.

        Returns:
            :returns B: numpy.ndarray of shape (nx,3), (nx,ny,3),(nx,ny,nz,3) \
                     for 1D, 2D and 3D runs, respectively
        """

        Bx = self._readfieldsdataset(time, 'Bx')
        By = self._readfieldsdataset(time, 'By')
        Bz = self._readfieldsdataset(time, 'Bz')

        return self._fillvec(Bx, By, Bz)
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetJ(self, time):

        """
        Creation : 2015-02-27 19:12:14.062301
        """

        Jx = self._readfieldsdataset(time, 'Jx')
        Jy = self._readfieldsdataset(time, 'Jy')
        Jz = self._readfieldsdataset(time, 'Jz')

        return self._fillvec(Jx, Jy, Jz)
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetE(self, time, origin_grid = False):

        """
        Creation : 2015-02-27 19:12:14.062301
        """

        Ex = self._readfieldsdataset(time, 'Ex')
        Ey = self._readfieldsdataset(time, 'Ey')
        Ez = self._readfieldsdataset(time, 'Ez')

        return self._fillvec(Ex, Ey, Ez)
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetCoords(self, extent = None, axis = 0):

        """
        Creation : 2015-02-27 19:12:14.062301

           returns an array with the coordinates of grid points along
           a given axis, for a given extent
        """

        dl = self.dl[axis]
        if extent is None:
            return dl * np.arange(self.ncells[axis]+1)
        else:
            npts = (extent[1]-extent[0])/dl + 1
            return extent[0] + dl*np.arange(npts)
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetV(self, time, species):

        """
        Returns the Bulk velocity of a specific species
        Exemple  : Ve = run.GetV(14.2, species='electrons')
        Creation : 2015-02-27 16:15:40.856938
        """

        if species.lower() in self._electron:
            Vxn = 'Vx[0]'
            Vyn = 'Vy[0]'
            Vzn = 'Vz[0]'

        elif species.lower() in self._proton:
            Vxn = 'Vx[1]'
            Vyn = 'Vy[1]'
            Vzn = 'Vz[1]'
        
        elif species.lower() in self._alpha:
            Vxn = 'Vx[2]'
            Vyn = 'Vy[2]'
            Vzn = 'Vz[2]'

        elif species.lower() in self._ion:
            Vxn = 'Vix'
            Vyn = 'Viy'
            Vzn = 'Viz'

        else:
            raise ValueError("Unknown species '%s'" % species)

        Vx = self._readfieldsdataset(time, Vxn)
        Vy = self._readfieldsdataset(time, Vyn)
        Vz = self._readfieldsdataset(time, Vzn)

        return self._fillvec(Vx, Vy, Vz)
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetN(self, time, species):

        """
        Returns the particle density of a specific species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938

        """

        proton   = ["protons", "proton", "p"]
        electron = ["electrons", "electron", "e"]
        ion      = ["ions", "ion", "i"]

        if species.lower() in self._electron:
            name = 'n[0]'

        elif species.lower() in self._proton:
            name = 'n[1]'

        elif species.lower() in self._ion:
            name = 'n[0]'

        elif species.lower() in self._alpha:
            name = 'n[2]'

        else:
            raise ValueError("Unknown species '%s'" % species)

        data = self._readfieldsdataset(time, name)

        return data
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetVSpecies(self, time, species):

        """
        Returns the velocity of a specific species
        Creation : 2018-05-28 
        """

        if species.lower() in self._proton:
            espece = 'specie 1'
        elif species.lower() in self._alpha:
            espece = 'specie 2'
        else:
            raise ValueError("Unknown species '%s'" % species)

        Vx = self._readspeciesdataset(time, espece, 'v[0]')
        Vy = self._readspeciesdataset(time, espece, 'v[1]')
        Vz = self._readspeciesdataset(time, espece, 'v[2]')

        return self._fillvec(Vx, Vy, Vz)
    #==========================================================



    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetRSpecies(self, time, species):

        """
        Returns the coordinates of a specific species
        Creation : 2018-05-28 
        """

        if species.lower() in self._proton:
            espece = 'specie 1'
        elif species.lower() in self._alpha:
            espece = 'specie 2'
        else:
            raise ValueError("Unknown species '%s'" % species)

        Rx = self._readspeciesdataset(time, espece, 'r[0]')
        Ry = self._readspeciesdataset(time, espece, 'r[1]')
        Rz = self._readspeciesdataset(time, espece, 'r[2]')

        return self._fillvec(Rx, Ry, Rz)
    #==========================================================




    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetIndexSpecies(self, time, species):

        """
        Returns the index of all the particles of a specific species
        Creation : 2018-05-28 
        """

        if species.lower() in self._proton:
            espece = 'specie 1'
        elif species.lower() in self._alpha:
            espece = 'specie 2'
        else:
            raise ValueError("Unknown species '%s'" % species)

        index = self._readspeciesdataset(time, espece, 'index')
        
        return index
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPxx(self, time, species):

        """
        Returns the component of the pressure tensor of a species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        if species.lower() in self._electron:
            name = 'Pxx[0]'

        elif species.lower() in self._proton:
            name = 'Pxx[1]'

        elif species.lower() in self._ion:
            name = 'Pxx[2]' # Add on the 16th of february to test 
            # full Pressure tensor not additive and not written by heckle
            #raise valueError("No pressure tensor for species %s" % species)

        else:
            raise ValueError("Unknown species '%s'" % species)

        data = self._readfieldsdataset(time, name)

        return data
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPxy(self, time, species):

        """
        Returns the component of the pressure tensor of a species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        if species.lower() in self._electron:
            name = 'Pxy[0]'

        elif species.lower() in self._proton:
            name = 'Pxy[1]'

        elif species.lower() in self._ion:
            # full Pressure tensor not additive and not written by heckle
            raise valueError("No pressure tensor for species %s" % species)

        else:
            raise ValueError("Unknown species '%s'" % species)

        data = self._readfieldsdataset(time, name)

        return data
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPxz(self, time, species):

        """
        Returns the component of the pressure tensor of a species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        if species.lower() in self._electron:
            name = 'Pxz[0]'

        elif species.lower() in self._proton:
            name = 'Pxz[1]'

        elif species.lower() in self._ion:
            # full Pressure tensor not additive and not written by heckle
            raise valueError("No pressure tensor for species %s" % species)

        else:
            raise ValueError("Unknown species '%s'" % species)

        data = self._readfieldsdataset(time, name)

        return data
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPyy(self, time, species):

        """
        Returns the component of the pressure tensor of a species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        if species.lower() in self._electron:
            name = 'Pyy[0]'

        elif species.lower() in self._proton:
            name = 'Pyy[1]'

        elif species.lower() in self._ion:
            name = 'Pxx[2]' # Add on the 16th of february to test 
            # full Pressure tensor not additive and not written by heckle
            #raise valueError("No pressure tensor for species %s" % species)

        else:
            raise ValueError("Unknown species '%s'" % species)

        data = self._readfieldsdataset(time, name)

        return data
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPyz(self, time, species):

        """
        Returns the component of the pressure tensor of a species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        if species.lower() in self._electron:
            name = 'Pyz[0]'

        elif species.lower() in self._proton:
            name = 'Pyz[1]'

        elif species.lower() in self._ion:
            # full Pressure tensor not additive and not written by heckle
            raise valueError("No pressure tensor for species %s" % species)

        else:
            raise ValueError("Unknown species '%s'" % species)

        data = self._readfieldsdataset(time, name)

        return data
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPzz(self, time, species):

        """
        Returns the component of the pressure tensor of a species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        if species.lower() in self._electron:
            name = 'Pzz[0]'

        elif species.lower() in self._proton:
            name = 'Pzz[1]'

        elif species.lower() in self._ion:
            # full Pressure tensor not additive and not written by heckle
            raise valueError("No pressure tensor for species %s" % species)

        else:
            raise ValueError("Unknown species '%s'" % species)

        data = self._readfieldsdataset(time, name)

        return data
    #==========================================================



    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetDriverXX(self, time):
        return self._readfieldsdataset(time, 'Drxx')

    def GetDriverXY(self, time):
        return self._readfieldsdataset(time, 'Drxy')

    def GetDriverXZ(self, time):
        return self._readfieldsdataset(time, 'Drxz')

    def GetDriverYY(self, time):
        return self._readfieldsdataset(time, 'Dryy')

    def GetDriverYZ(self, time):
        return self._readfieldsdataset(time, 'Dryz')

    def GetDriverZZ(self, time):
        return self._readfieldsdataset(time, 'Drzz')
    #==========================================================

    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetMass(self, species):

        """
        Returns the mass of particles of a particular species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938

        """
        # these values should be read from attribute of species.h5
        if species.lower() in self._electron:
            index = 0
        elif species.lower() in self._proton:
            index = 1
        elif species.lower() in self._alpha:
            index = 2
        else:
            raise ValueError("Unknown species '%s'" % species)
        return self.mass[index]
    #==========================================================




    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetCharge(self, species):

        """
        Returns the electric charge of particles of a particular species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938

        """
        # these values should be read from attribute of species.h5
        if species.lower() in self._electron:
            return -1.
        elif species.lower() in self._proton:
            return 1.
        else:
            raise ValueError("Unknown species '%s'" % species)
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetHyperResistivity(self):

        """
        @todo: Brief Docstring for GetHyperResistivity
        Longer description here

        @return: @todo

        Exemple  :

        Creation : 2015-03-03

        """
        return self.hprsty
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetResistivity(self):

        """
        @todo: Brief Docstring for GetHyperResistivity
        Longer description here

        @return: @todo

        Exemple  :

        Creation : 2015-03-03

        """

        return self.rsty
    #==========================================================









    #----------------------------------------------------------
    #----------------------------------------------------------
    def indices2coord(self, indices, quantity=None):

        """
        Converts an array of shape (2, N) in physical coordinates
        normalized to the ion inertial length.

        Args :
            quantity : Not used in the Heckle code since all
            quantities are defined on the same grid.

        Exemple  :

        Creation: 2015-10-08
        """

        # could more general ? indices could be a np array of whatever rank
        # (1, 2 or 3) and whatever size... even if right now, it only seems
        # to be used in recrate.py (and picgsfc)
        coords = np.ndarray(indices.shape)
        ndim   = indices.shape[0]

        for c in np.arange(ndim):
            coords[c,:] = self.dl[c]*indices[c,:]

        return coords
    #==========================================================




    def index2Coord(self,
                    index,
                    axis) :

        """
        index : 1d numpy array containing a set of indices

        axis : 0, 1 or 2 is the direction along which indices are given

        output : a 1d numpy array of same shape as index withj the coordinates
                 in physical units
        """

        index = np.array(index)
        coord = np.empty(index.shape, dtype = float)

        if axis == 0 or axis == 1 or axis == 2 :
           coord = index*self.dl[axis]
        else :
           raise ValueError("axis value is mandatory & can only be 0, 1 or 2")

        return coord


    def coord2Index(self,
                    coord,
                    axis) :

        """
        coord : 1d numpy array containing a set of coordinates

        axis : 0, 1 or 2 is the direction along which indices are given

        output : a 1d numpy array of same shape as coord with the associateds
                 indices on G1 grid
        """

        coord = np.array(coord)
        index = np.empty(coord.shape, dtype = int)

        if axis == 0 or axis == 1 or axis == 2 :
           index = int(coord/self.dl[axis])
        else :
           raise ValueError("axis value is mandatory & can only be 0, 1 or 2")

        return index


    def spawnAxis(self,
                  index,
                  coord,
                  axis) :
        """
        """

        # spawn axis from a given indices
        if index is not None :
           if index.__len__() == 2 :
              myaxis = np.linspace(index[0], index[1], index[1]-index[0]+1)
           else :
              raise ValueError("index has to be of length 2")

        # spawn axis from given coordinates
        if coord is not None :
           if coord.__len__() == 2 :
              if axis == 0 or axis == 1 or axis == 2 :
                 minbound = int(coord[0]/self.dl[axis])*self.dl[axis]
                 maxbound = int(coord[1]/self.dl[axis])*self.dl[axis]
                 myaxis = np.arange(minbound, maxbound, self.dl[axis])
              else :
                 raise ValueError("axis has to be specified")
           else :
              raise ValueError("coord has to be of length 2")

        return myaxis



    def fourierFlux(self, time) :

       import pyfftw
       from pylab import meshgrid


       numOfThreads = 8

       dim = self.ncells+1

       # arrays have to be aligned... for pyfftw to work
       bx = pyfftw.empty_aligned(dim, dtype = 'complex128')
       by = pyfftw.empty_aligned(dim, dtype = 'complex128')
       BX = pyfftw.empty_aligned(dim, dtype = 'complex128')
       BY = pyfftw.empty_aligned(dim, dtype = 'complex128')
       AZ = pyfftw.empty_aligned(dim, dtype = 'complex128')
       az = pyfftw.empty_aligned(dim, dtype = 'complex128')

       # 2d fftw plans
       fftForx = pyfftw.FFTW(bx, BX, axes=(0, 1), direction='FFTW_FORWARD', flags=('FFTW_MEASURE', ))
       fftFory = pyfftw.FFTW(by, BY, axes=(0, 1), direction='FFTW_FORWARD', flags=('FFTW_MEASURE', ))
       fftBack = pyfftw.FFTW(AZ, az, axes=(0, 1), direction='FFTW_BACKWARD', flags=('FFTW_MEASURE', ))

       # fill the input arrays
       bx.real[...] = self.GetB(time)[..., 0]
       by.real[...] = self.GetB(time)[..., 1]
       bx.imag = np.zeros(dim, dtype = float)
       by.imag = np.zeros(dim, dtype = float)

       # update needed before execute plan
       fftForx.update_arrays(bx, BX)
       fftFory.update_arrays(by, BY)

       pyfftw.FFTW.execute(fftForx)
       pyfftw.FFTW.execute(fftFory)

       # this intend to build the k values, first half is ascending
       # positives values, last one being the descending negative values
       axis = 0
       N = dim[axis]
       k = np.linspace(0, N-1, N)
       kplus = k[:int(N/2)+1]
       kminus =-k[int((N-1)/2):0:-1]
       kx = np.concatenate((kplus, kminus))*2*np.pi/self.domsize[axis]

       axis = 1
       N = dim[axis]
       k = np.linspace(0, N-1, N)
       kplus = k[:int(N/2)+1]
       kminus =-k[int((N-1)/2):0:-1]
       ky = np.concatenate((kplus, kminus))*2*np.pi/self.domsize[axis]

       # create 2d values to avoid nested "for" loops
       wy, wx = meshgrid(ky, kx)

       k2 = np.square(wx)+np.square(wy)
       k2[0][0] = 1.0

       AZ = np.divide(1j*(wx*BY-wy*BX), k2)

       # remove dc component while k=0 mode is not defined
       AZ[0][0] = 0.0

       fftBack.update_arrays(AZ, az)

       pyfftw.FFTW.execute(fftBack)

       _flux = az.real/(self.ncells[0]*self.ncells[1])

       return _flux
