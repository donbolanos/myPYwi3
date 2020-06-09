

import numpy as np


# class interface :
# ------------------

#   - modelname     - gives the name of the model (pic, hybrid, etc.)
#   - display()     - display some informations about the run
#   - GetCoords()   - returns coordinate along a direction with default extent
#   - GetNbDiag()   - returns the number of dumps for fields/part. data
#   - GetMass()     - mass of 'species'
#   - GetCharge()   - charge of 'species'
#   - GetRunTime()  - returns total time of the run (even not finished)
#   - GetTimeStep() - dump step for fields, particles
#   - whereAmI()    - where is the run stored

#   - indices2coord()
#       - based on indices of the quantity, return the coordinates
#   - coord2indices()


#   - GetVe()
#   - GetVi()
#   - GetV(species)
#   - GetB(, original_grid=False)
#   - GetE(, "")
#   - GetJ()
#   - GetNe()
#   - GetNi()
#   - GetN(species)
#   - GetFlux(warning if 3D)
#   - GetP(species)
#   - GetPxx(species)
#   - GetPxy(species)
#   - GetPxz(species)
#   - GetPyy(species)
#   - GetPyz(species)
#   - GetPzz(species)
#   - GetThermalEnergy(species) (utherm)
#   - VgradV(species)
#   - divP()
#   - GetParticles(time, species, location=None)


# model dependant
# ------------------
#   - GetDebye()


# should be private :
# ------------------

#   - iter2time()
#   - time2iter()
#       - all should be accessed by time

#   - GetField(fieldname..., time)
#   - GetFieldData()
#   - divT(T)



# obsoletes ?

#   - GetFIeldFiles()
#   - GetFileTime()
#       - not valid anymore for 1 file HDF5 groups
#   - GetNbFieldFiles()

#   - GetOutOfPlaneDir() ?
#       - should always put the ignorable coordinate in Z (3rd coord)




class Run (object):

    """
    The :class:`Run` is the base class for all run objects implemented
    in pywi. It is not intended to be instanciated but defines the interface
    that all the run objects must implement. Some routines are already
    implemented in :class:`Run` (such as GetVp()) since they will be
    the same for all runs.


    In this class, should be methods that concern only basic/raw quantities
    that are available in data files or specific run properties. For instance,
    a routine returning the pressure tensor is typically one that will be
    implemented here, while one returning the reconnection rate or the degree
    of nongyrotropy will not, those are higher-level products that can be
    calculated the same way for all runs based on basics quantities

    """



    #----------------------------------------------------------
    #----------------------------------------------------------
    def __init__(self, path, runID = 'r'):
        """Exemple  :

        Creation : 2015-02-27 16:04:58.868849

        """
        self.modelname  = 'Unknown'
        self.codename   = 'Unknown code'
        self.path       = path
        self.runID      = runID
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def display(self):
        """
        Shows some information on screen about the run, like its path, 
        its ID, and the name of the physical model.
        """
        print 'Model name       : %s '  % (self.modelname)
        print 'Run ID           : %s'   % (self.runID)
        print 'Path             : %s'   % (self.path)
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def whereAmI(self):
        """ Exemple  : path = run.whereAmI()
        """
        return self.path
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetVe(self, time):
        """
        Returns the electron Bulk velocity at desired time

        Parameters:

            :time: is the time in inverse ion cyclotron frequency

        :Returns: numpy.ndarray, axis -1 is the component of the vector

        Exemple:

           >>> Ve = run.GetVe(0.1)

        """
        return self.GetV(time, species='electrons')
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetVi(self, time):
        """
        Returns the Ion Bulk velocity at desired time. Ions means
        that it is the average of ALL ion species

        Parameters:

            :time: is the time in inverse ion cyclotron frequency

        :Returns: numpy.ndarray, axis -1 is the component of the vector

        Exemple:

           >>> Ve = run.GetVe(0.1)

        """
        return self.GetV(time, species='ions')
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetVp(self, time):
        """
        Returns the proton Bulk velocity

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        return self.GetV(time, species='protons')
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetV(self, time, species):
        """
        Returns the Bulk velocity of a specific species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetB(self, time, original_grid = False):
        """
        Returns the magnetic field at a specific time
        original_grid = True returns the field on the
        grid that is used in the simulation. This is not
        necessarily implemented for all models

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================




    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetE(self, time, original_grid = False):
        """
        Returns the electric field at a specific time
        original_grid = True returns the field on the
        grid that is used in the simulation. This is not
        necessarily implemented for all models

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================




    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetJ(self, time):
        """
        Returns the electric current density vector

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetNe(self, time):
        """
        Returns the electron particle density

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        return self.GetN(time, species='electrons')
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetNi(self, time):
        """
        Returns the ion particle density, this include all
        ion species.

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        return self.GetN(time, species='ions')
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetNp(self, time):
        """
        Returns the proton particle density
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938

        """
        return self.GetN(time, species='protons')
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetN(self, time, species):
        """
        Returns the particle density of a specific species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetFLux(self, time):
        """
        Returns the magnetic flux function
        this raises an exception if the run is not 2D

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetP(self, time, species):
        """
        Returns the pressure tensor of a specific species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        Pxx = self.GetPxx(time, species)
        Pxy = self.GetPxy(time, species)
        Pxz = self.GetPxz(time, species)
        Pyy = self.GetPyy(time, species)
        Pyz = self.GetPyz(time, species)
        Pzz = self.GetPzz(time, species)

        P  = np.ndarray(dtype=self.dtype, shape=self.tfield_shape)
        P[...,0,0] = Pxx
        P[...,1,1] = Pyy
        P[...,2,2] = Pzz

        P[...,0,1] = Pxy
        P[...,0,2] = Pxz
        P[...,1,2] = Pyz

        P[...,1,0] = Pxy
        P[...,2,0] = Pxz
        P[...,2,1] = Pyz

        return P

    #==========================================================


    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetD(self, time):
        Dxx = self.GetDriverXX(time)
        Dxy = self.GetDriverXY(time)
        Dxz = self.GetDriverXZ(time)
        Dyy = self.GetDriverYY(time)
        Dyz = self.GetDriverYZ(time)
        Dzz = self.GetDriverZZ(time)

        D  = np.ndarray(dtype=self.dtype, shape=self.tfield_shape)
        D[...,0,0] = Dxx
        D[...,1,1] = Dyy
        D[...,2,2] = Dzz

        D[...,0,1] = Dxy
        D[...,0,2] = Dxz
        D[...,1,2] = Dyz

        D[...,1,0] = Dxy
        D[...,2,0] = Dxz
        D[...,2,1] = Dyz

        return D

    #==========================================================

    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPe(self, time):
        """
        Returns the electron pressure tensor
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """
        return self.GetP(time, species='electrons')
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPp(self, time):
        """
        Returns the proton pressure tensor
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """
        return self.GetP(time, species='protons')
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPxx(self, time, species):
        """
        Returns the component of the pressure tensor of a species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPxy(self, time, species):
        """
        Returns the component of the pressure tensor of a species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================








    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPxz(self, time, species):
        """
        Returns the component of the pressure tensor of a species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPyy(self, time, species):
        """
        Returns the component of the pressure tensor of a species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPyz(self, time, species):
        """
        Returns the component of the pressure tensor of a species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPzz(self, time, species):
        """
        Returns the component of the pressure tensor of a species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetThermalEnergy(self, time, species):
        """
        Returns the thermal energy for a specific species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """
        P = self.GetP(time, species=species)
        return 0.5 * (P[0,0,:,:,:] + P[1,1,:,:,:] + P[2,2,:,:,:])
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def VgradV(self, time, species):
        """
        Returns the (V dot Grad )(v)

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def divP(self, time, species):
        """
        Returns the divergence of the pressure tensor of a
        specific species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def __divT(self, time, T):
        """
        Returns the divergence of the tensor T

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetParticles(self, time, species, location=None):
        """
        Returns a group of particles from a particular species
        at some location for some time
        If no location is not specified, it returns the whole
        box

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetMass(self, species):
        """
        Returns the mass of particles of a particular species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetCharge(self, species):
        """
        Returns the electric charge of particles of a particular species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetFlux(self, time):

        #import matplotlib.pyplot as plt

        B = self.GetB(time)

        if len(B.shape) == 4:
            raise ValueError("magnetic Flux() only defined for 2D")

        shape = (B.shape[0], B.shape[1])
        flux  = np.zeros(shape=shape, dtype=B.dtype)

        dl1 = self.dl[0]
        dl2 = self.dl[1]

        n1 = shape[0]
        n2 = shape[1]

        b1 = B[...,0]
        b2 = B[...,1]


        flux[1:,0] = flux[0,0] + np.cumsum(b2[:-1,0]*dl1)
        flux[0,1:] = flux[0,0] - np.cumsum(b1[0,:-1]*dl2)


        flux = self._fast(flux, b1, b2, dl1, dl2)

        return flux
    #==========================================================







    # SEBERG (#scipy) method to calculate the 2D integral flux from B
    #http://stackoverflow.com/questions/11854522/nested-loop-with-array-indexing-in-numpy
    def _fast(self, flux, b1, b2, dl1=1., dl2=1.):

        from scipy.ndimage import convolve

        _flux = np.zeros((flux.shape[0]+1, flux.shape[1]+1), dtype=flux.dtype)
        temp_b1 = np.zeros((b1.shape[0]+1, b1.shape[1]+1), dtype=b1.dtype)
        temp_b2 = np.zeros((b2.shape[0]+1, b2.shape[1]+1), dtype=b2.dtype)

        _flux[:-1,:-1] = flux
        convolve(_flux[:-1,:-1], [[0, 0.5], [0.5, 0]], _flux[1:,1:])

        temp_b2[1:,1:-1] = b2[:,1:]*dl1
        temp_b1[1:-1,1:] = b1[1:,:]*dl2

        conv_b = np.array([[0.0, 0.5], [0.5, 0.5]])
        convolve(temp_b2[:-1,:-1], [[0.5, 0.5], [0.5, 0.]], temp_b2[1:,1:])
        convolve(temp_b1[:-1,:-1], [[-0.5, 0.5], [0.5, 0.]], temp_b1[1:,1:])

        _flux += temp_b2
        _flux += temp_b1

        return _flux[:-1,:-1]

