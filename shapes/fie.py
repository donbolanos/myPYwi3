#
#! /usr/bin/env python
#
#
import numpy as np
#
#
class fie() :

    """
    class associated to a "field" data
    """

    def __init__(self,
                 run = None,
                 field = None,
                 plane = None,
                 direction = None,
                 cut = None,
                 time = 0.0,
                 ylabel = None,
                 extent = None,
                 shifts = None,
                 bounds = None):

        """
        """

        # .. the run
        self.run = run

        # .. the data
        self.xdata = None
        self.ydata = None
        self.zdata = None

        # .. name of the field to display
        self.field = field

        # .. the index for the 2d cut
        if plane is not None :
          self.plane = plane.lower()
          if   (self.plane in ["xy", "yx"]) : self.cut = np.rint(cut/self.run.dl[2])
          elif (self.plane in ["xz", "zx"]) : self.cut = np.rint(cut/self.run.dl[1])
          elif (self.plane in ["yz", "zy"]) : self.cut = np.rint(cut/self.run.dl[0])
          else: print "wrong plane..."
        else : self.plane = None

        # .. the indexes for the 1d cut
        if direction is not None :
          self.direction = direction.lower()
          if   self.direction == "x" : self.idx = [np.rint(cut[0]/self.run.dl[1]), np.rint(cut[1]/self.run.dl[2])]
          elif self.direction == "y" : self.idx = [np.rint(cut[0]/self.run.dl[0]), np.rint(cut[1]/self.run.dl[2])]
          elif self.direction == "z" : self.idx = [np.rint(cut[0]/self.run.dl[0]), np.rint(cut[1]/self.run.dl[1])]
          else: print "wrong direction..."
        else : self.direction = None

        # .. set the extent of the domain to plot
        self.extent = extent

        # .. set the shifts for the domain to plot
        if shifts == None : self.shifts = [0, 0]
        else : self.shifts = shifts

        # .. set the time of the fig
        self.time = run.SetTime(time)

        # ..
        self.ylabel = ylabel

        # ..
        self.bounds = bounds

        # .. set the name of the saved file
        self.filename = 'fie'
        #self.filename = 'f' + '%d' %run.Time2Iter(self.time, 'f')

        # .. the labels of the fig
        self.MakeLabels()

        # .. get a 1d field with the right shape
        if self.direction is not None: self.Shape1D()

        # .. get a 2d field with the right shape
        if self.plane is not None: self.Shape2D()


    def Shape2D(self):

        """
        reshape "field" data in "plane" sliced at "cut" on "domain" @ "time"
        """

        # .. fill "rawdata" @ the given "time"
        if   self.field == "db"   : rawdata = self.run.GetDivB(time = self.time)[0]
        elif self.field == "bt"   : rawdata = self.run.GetModulusB(time = self.time)[0]
        elif self.field == "pm"   : rawdata = self.run.GetMagneticPressure(time = self.time)[0]
        elif self.field == "n0"   : rawdata = self.run.GetElectronDensity(time = self.time)[0]
        elif self.field == "n1"   : rawdata = self.run.GetProtonDensity(time = self.time)[0]
        elif self.field == "n2"   : rawdata = self.run.GetAlphaDensity(time = self.time)[0]
        elif self.field == "v0"   : rawdata = self.run.GetModulusElectronVelocity(time = self.time)[0]
        elif self.field == "v1"   : rawdata = self.run.GetModulusProtonVelocity(time = self.time)[0]
        elif self.field == "v2"   : rawdata = self.run.GetModulusAlphaVelocity(time = self.time)[0]
        elif self.field == "vi"   : rawdata = self.run.GetModulusIonVelocity(time = self.time)[0]
        elif self.field == "p0"   : rawdata = self.run.GetElectronKineticPressure(time = self.time)[0]
        elif self.field == "p1"   : rawdata = self.run.GetProtonKineticPressure(time = self.time)[0]
        elif self.field == "p2"   : rawdata = self.run.GetAlphaKineticPressure(time = self.time)[0]
        elif self.field == "pi"   : rawdata = self.run.GetIonKineticPressure(time = self.time)[0]
        elif self.field == "pk"   : rawdata = self.run.GetTotalKineticPressure(time = self.time)[0]
        elif self.field == "t0"   : rawdata = self.run.GetElectronTemperature(time = self.time)[0]
        elif self.field == "t1"   : rawdata = self.run.GetProtonTemperature(time = self.time)[0]
        elif self.field == "t2"   : rawdata = self.run.GetAlphaTemperature(time = self.time)[0]
        elif self.field == "ti"   : rawdata = self.run.GetIonTemperature(time = self.time)[0]
        elif self.field == "r1"   : rawdata = self.run.GetProtonRamPressure(time = self.time)[0]
        elif self.field == "r2"   : rawdata = self.run.GetAlphaRamPressure(time = self.time)[0]
        elif self.field == "ri"   : rawdata = self.run.GetIonRamPressure(time = self.time)[0]
        elif self.field == "jde"  : rawdata = self.run.GetJdotE(time = self.time)[0]
        elif self.field == "psi"  : rawdata = self.run.GetPsiAngle(time = self.time)[0]
        elif self.field == "j0"   : rawdata = self.run.GetModulusElectronCurrent(time = self.time)[0]
        elif self.field == "j1"   : rawdata = self.run.GetModulusprotonCurrent(time = self.time)[0]
        elif self.field == "j2"   : rawdata = self.run.GetModulusAlphaCurrent(time = self.time)[0]
        elif self.field == "jt"   : rawdata = self.run.GetModulusTotalCurrent(time = self.time)[0]
        elif self.field == "bx"   : rawdata = self.run.GetB(time = self.time)[0]
        elif self.field == "by"   : rawdata = self.run.GetB(time = self.time)[1]
        elif self.field == "bz"   : rawdata = self.run.GetB(time = self.time)[2]
        elif self.field == "bxy"  : rawdata = self.run.GetB(time = self.time)[3]
        elif self.field == "ax"   : rawdata = self.run.GetA(time = self.time)[0]
        elif self.field == "ay"   : rawdata = self.run.GetA(time = self.time)[1]
        elif self.field == "az"   : rawdata = self.run.GetA(time = self.time)[2]
        elif self.field == "ex"   : rawdata = self.run.GetE(time = self.time)[...,0]
        elif self.field == "ey"   : rawdata = self.run.GetE(time = self.time)[...,1]
        elif self.field == "ez"   : rawdata = self.run.GetE(time = self.time)[...,2]
        elif self.field == "exy"  : rawdata = self.run.GetE(time = self.time)[...,3]
        elif self.field == "exy"  : rawdata = self.run.GetE(time = self.time)[3]
        elif self.field == "jx"   : rawdata = self.run.GetJ(time = self.time)[0]
        elif self.field == "jy"   : rawdata = self.run.GetJ(time = self.time)[1]
        elif self.field == "jz"   : rawdata = self.run.GetJ(time = self.time)[2]
        elif self.field == "jxy"  : rawdata = self.run.GetJ(time = self.time)[3]
        elif self.field == "vix"  : rawdata = self.run.GetIonVelocity(time = self.time)[0]
        elif self.field == "viy"  : rawdata = self.run.GetIonVelocity(time = self.time)[1]
        elif self.field == "viz"  : rawdata = self.run.GetIonVelocity(time = self.time)[2]
        elif self.field == "v0x"  : rawdata = self.run.GetElectronVelocity(time = self.time)[0]
        elif self.field == "v0y"  : rawdata = self.run.GetElectronVelocity(time = self.time)[1]
        elif self.field == "v0z"  : rawdata = self.run.GetElectronVelocity(time = self.time)[2]
        elif self.field == "v0xy" : rawdata = self.run.GetElectronVelocity(time = self.time)[3]
        elif self.field == "v1x"  : rawdata = self.run.GetProtonVelocity(time = self.time)[0]
        elif self.field == "v1y"  : rawdata = self.run.GetProtonVelocity(time = self.time)[1]
        elif self.field == "v1z"  : rawdata = self.run.GetProtonVelocity(time = self.time)[2]
        elif self.field == "v2x"  : rawdata = self.run.GetAlphaVelocity(time = self.time)[0]
        elif self.field == "v2y"  : rawdata = self.run.GetAlphaVelocity(time = self.time)[1]
        elif self.field == "v2z"  : rawdata = self.run.GetAlphaVelocity(time = self.time)[2]
        elif self.field == "j0x"  : rawdata = self.run.GetElectronCurrent(time = self.time)[0]
        elif self.field == "j0y"  : rawdata = self.run.GetElectronCurrent(time = self.time)[1]
        elif self.field == "j0z"  : rawdata = self.run.GetElectronCurrent(time = self.time)[2]
        elif self.field == "j1x"  : rawdata = self.run.GetProtonCurrent(time = self.time)[0]
        elif self.field == "j1y"  : rawdata = self.run.GetProtonCurrent(time = self.time)[1]
        elif self.field == "j1z"  : rawdata = self.run.GetProtonCurrent(time = self.time)[2]
        elif self.field == "j2x"  : rawdata = self.run.GetAlphaCurrent(time = self.time)[0]
        elif self.field == "j2y"  : rawdata = self.run.GetAlphaCurrent(time = self.time)[1]
        elif self.field == "j2z"  : rawdata = self.run.GetAlphaCurrent(time = self.time)[2]
        elif self.field == "idx"  : rawdata = self.run.GetVcrossB(time = self.time)[0]
        elif self.field == "idy"  : rawdata = self.run.GetVcrossB(time = self.time)[1]
        elif self.field == "idz"  : rawdata = self.run.GetVcrossB(time = self.time)[2]
        elif self.field == "idxy" : rawdata = self.run.GetVcrossB(time = self.time)[3]
        elif self.field == "hax"  : rawdata = self.run.GetHall(time = self.time)[0]
        elif self.field == "hay"  : rawdata = self.run.GetHall(time = self.time)[1]
        elif self.field == "haz"  : rawdata = self.run.GetHall(time = self.time)[2]
        elif self.field == "haxy" : rawdata = self.run.GetHall(time = self.time)[3]
        elif self.field == "gp0x" : rawdata = self.run.GetDivElectronPressure(time = self.time)[0]
        elif self.field == "gp0y" : rawdata = self.run.GetDivElectronPressure(time = self.time)[1]
        elif self.field == "gp0z" : rawdata = self.run.GetDivElectronPressure(time = self.time)[2]
        elif self.field == "gp0xy": rawdata = self.run.GetDivElectronPressure(time = self.time)[3]
        elif self.field == "gp1x" : rawdata = self.run.GetDivProtonPressure(time = self.time)[0]
        elif self.field == "gp1y" : rawdata = self.run.GetDivProtonPressure(time = self.time)[1]
        elif self.field == "gp1z" : rawdata = self.run.GetDivProtonPressure(time = self.time)[2]
        elif self.field == "gp2x" : rawdata = self.run.GetDivAlphaPressure(time = self.time)[0]
        elif self.field == "gp2y" : rawdata = self.run.GetDivAlphaPressure(time = self.time)[1]
        elif self.field == "gp2z" : rawdata = self.run.GetDivAlphaPressure(time = self.time)[2]
        elif self.field == "gpix" : rawdata = self.run.GetDivIonPressure(time = self.time)[0]
        elif self.field == "gpiy" : rawdata = self.run.GetDivIonPressure(time = self.time)[1]
        elif self.field == "gpiz" : rawdata = self.run.GetDivIonPressure(time = self.time)[2]
        elif self.field == "rex"  : rawdata = self.run.GetResistivity(time = self.time)[0]
        elif self.field == "rey"  : rawdata = self.run.GetResistivity(time = self.time)[1]
        elif self.field == "rez"  : rawdata = self.run.GetResistivity(time = self.time)[2]
        elif self.field == "hvx"  : rawdata = self.run.GetHyperViscosity(time = self.time)[0]
        elif self.field == "hvy"  : rawdata = self.run.GetHyperViscosity(time = self.time)[1]
        elif self.field == "hvz"  : rawdata = self.run.GetHyperViscosity(time = self.time)[2]
        elif self.field == "p0xx" : rawdata = self.run.GetFullElectronPressure(time = self.time)[0]
        elif self.field == "p0xy" : rawdata = self.run.GetFullElectronPressure(time = self.time)[1]
        elif self.field == "p0xz" : rawdata = self.run.GetFullElectronPressure(time = self.time)[2]
        elif self.field == "p0yy" : rawdata = self.run.GetFullElectronPressure(time = self.time)[3]
        elif self.field == "p0yz" : rawdata = self.run.GetFullElectronPressure(time = self.time)[4]
        elif self.field == "p0zz" : rawdata = self.run.GetFullElectronPressure(time = self.time)[5]
        elif self.field == "p1xx" : rawdata = self.run.GetFullprotonPressure(time = self.time)[0]
        elif self.field == "p1xy" : rawdata = self.run.GetFullprotonPressure(time = self.time)[1]
        elif self.field == "p1xz" : rawdata = self.run.GetFullprotonPressure(time = self.time)[2]
        elif self.field == "p1yy" : rawdata = self.run.GetFullprotonPressure(time = self.time)[3]
        elif self.field == "p1yz" : rawdata = self.run.GetFullprotonPressure(time = self.time)[4]
        elif self.field == "p1zz" : rawdata = self.run.GetFullprotonPressure(time = self.time)[5]
        elif self.field == "p2xx" : rawdata = self.run.GetFullIonPressure(time = self.time)[0]
        elif self.field == "p2xy" : rawdata = self.run.GetFullIonPressure(time = self.time)[1]
        elif self.field == "p2xz" : rawdata = self.run.GetFullIonPressure(time = self.time)[2]
        elif self.field == "p2yy" : rawdata = self.run.GetFullIonPressure(time = self.time)[3]
        elif self.field == "p2yz" : rawdata = self.run.GetFullIonPressure(time = self.time)[4]
        elif self.field == "p2zz" : rawdata = self.run.GetFullIonPressure(time = self.time)[5]
        elif self.field == "pixx" : rawdata = self.run.GetFullIonPressure(time = self.time)[0]
        elif self.field == "pixy" : rawdata = self.run.GetFullIonPressure(time = self.time)[1]
        elif self.field == "pixz" : rawdata = self.run.GetFullIonPressure(time = self.time)[2]
        elif self.field == "piyy" : rawdata = self.run.GetFullIonPressure(time = self.time)[3]
        elif self.field == "piyz" : rawdata = self.run.GetFullIonPressure(time = self.time)[4]
        elif self.field == "pizz" : rawdata = self.run.GetFullIonPressure(time = self.time)[5]
        elif self.field == "ncez" : rawdata = self.run.GetNCEZ(time = self.time)[0]
        elif self.field == "jcbx" : rawdata = self.run.GetJcrossB(time = self.time)[0]
        elif self.field == "jcby" : rawdata = self.run.GetJcrossB(time = self.time)[1]
        elif self.field == "jcbz" : rawdata = self.run.GetJcrossB(time = self.time)[2]
        elif self.field == "jxby" : rawdata = self.run.GetJcrossB(time = self.time)[3]
        elif self.field == "jxbz" : rawdata = self.run.GetJcrossB(time = self.time)[4]
        elif self.field == "jybx" : rawdata = self.run.GetJcrossB(time = self.time)[5]
        elif self.field == "jybz" : rawdata = self.run.GetJcrossB(time = self.time)[6]
        elif self.field == "jzbx" : rawdata = self.run.GetJcrossB(time = self.time)[7]
        elif self.field == "jzby" : rawdata = self.run.GetJcrossB(time = self.time)[8]
        else : print "invalid field..."

        # .. set xdata, ydata & zdata
        if   self.plane == "xy" :
          self.xdata = np.linspace(0+self.shifts[0], self.run.l[0]+self.shifts[0], self.run.n[0]+1)
          self.ydata = np.linspace(0+self.shifts[1], self.run.l[1]+self.shifts[1], self.run.n[1]+1)
          self.zdata = rawdata[:, :, self.cut]
        elif self.plane == "xz" :
          self.zdata = rawdata[:, self.cut, :]
          self.xdata = np.linspace(0+self.shifts[0], self.run.l[0]+self.shifts[0], self.run.n[0]+1)
          self.ydata = np.linspace(0+self.shifts[1], self.run.l[2]+self.shifts[1], self.run.n[2]+1)
        elif self.plane == "yz" :
          self.zdata = rawdata[self.cut, :, :]
          self.xdata = np.linspace(0+self.shifts[0], self.run.l[1]+self.shifts[0], self.run.n[1]+1)
          self.ydata = np.linspace(0+self.shifts[1], self.run.l[2]+self.shifts[1], self.run.n[2]+1)
        else : print "still wrong direction"

        # .. set extent if not imposed
        if self.extent == None : self.extent = [[np.min(self.xdata), np.max(self.xdata)], [np.min(self.ydata), np.max(self.ydata)]]
        else : self.extent = [[self.extent[0][0]+self.shifts[0], self.extent[0][1]+self.shifts[0]], [self.extent[1][0]+self.shifts[1], self.extent[1][1]+self.shifts[1]]]

        # .. set bounds
        if self.bounds is None : self.bounds = [np.min(self.zdata), np.max(self.zdata)]


    def Shape1D(self):

        """
        reshape "field" in "direction" sliced at "cut" on "domain" @ "time"
        """

        # .. fill "rawdata" @ the given "time"
        if   self.field == "db"   : rawdata = self.run.GetDivB(time = self.time)[0]
        elif self.field == "bt"   : rawdata = self.run.GetModulusB(time = self.time)[0]
        elif self.field == "pm"   : rawdata = self.run.GetMagneticPressure(time = self.time)[0]
        elif self.field == "n0"   : rawdata = self.run.GetElectronDensity(time = self.time)[0]
        elif self.field == "n1"   : rawdata = self.run.GetProtonDensity(time = self.time)[0]
        elif self.field == "n2"   : rawdata = self.run.GetAlphaDensity(time = self.time)[0]
        elif self.field == "v0"   : rawdata = self.run.GetModulusElectronVelocity(time = self.time)[0]
        elif self.field == "v1"   : rawdata = self.run.GetModulusProtonVelocity(time = self.time)[0]
        elif self.field == "v2"   : rawdata = self.run.GetModulusAlphaVelocity(time = self.time)[0]
        elif self.field == "vi"   : rawdata = self.run.GetModulusIonVelocity(time = self.time)[0]
        elif self.field == "p0"   : rawdata = self.run.GetElectronKineticPressure(time = self.time)[0]
        elif self.field == "p1"   : rawdata = self.run.GetProtonKineticPressure(time = self.time)[0]
        elif self.field == "p2"   : rawdata = self.run.GetAlphaKineticPressure(time = self.time)[0]
        elif self.field == "pi"   : rawdata = self.run.GetIonKineticPressure(time = self.time)[0]
        elif self.field == "pk"   : rawdata = self.run.GetTotalKineticPressure(time = self.time)[0]
        elif self.field == "t0"   : rawdata = self.run.GetElectronTemperature(time = self.time)[0]
        elif self.field == "t1"   : rawdata = self.run.GetProtonTemperature(time = self.time)[0]
        elif self.field == "t2"   : rawdata = self.run.GetAlphaTemperature(time = self.time)[0]
        elif self.field == "ti"   : rawdata = self.run.GetIonTemperature(time = self.time)[0]
        elif self.field == "r1"   : rawdata = self.run.GetProtonRamPressure(time = self.time)[0]
        elif self.field == "r2"   : rawdata = self.run.GetAlphaRamPressure(time = self.time)[0]
        elif self.field == "ri"   : rawdata = self.run.GetIonRamPressure(time = self.time)[0]
        elif self.field == "jde"  : rawdata = self.run.GetJdotE(time = self.time)[0]
        elif self.field == "psi"  : rawdata = self.run.GetPsiAngle(time = self.time)[0]
        elif self.field == "j0"   : rawdata = self.run.GetModulusElectronCurrent(time = self.time)[0]
        elif self.field == "j1"   : rawdata = self.run.GetModulusprotonCurrent(time = self.time)[0]
        elif self.field == "j2"   : rawdata = self.run.GetModulusAlphaCurrent(time = self.time)[0]
        elif self.field == "jt"   : rawdata = self.run.GetModulusTotalCurrent(time = self.time)[0]
        elif self.field == "bx"   : rawdata = self.run.GetB(time = self.time)[0]
        elif self.field == "by"   : rawdata = self.run.GetB(time = self.time)[1]
        elif self.field == "bz"   : rawdata = self.run.GetB(time = self.time)[2]
        elif self.field == "bxy"  : rawdata = self.run.GetB(time = self.time)[3]
        elif self.field == "ax"   : rawdata = self.run.GetA(time = self.time)[0]
        elif self.field == "ay"   : rawdata = self.run.GetA(time = self.time)[1]
        elif self.field == "az"   : rawdata = self.run.GetA(time = self.time)[2]
        elif self.field == "ex"   : rawdata = self.run.GetE(time = self.time)[0]
        elif self.field == "ey"   : rawdata = self.run.GetE(time = self.time)[1]
        elif self.field == "ez"   : rawdata = self.run.GetE(time = self.time)[2]
        elif self.field == "exy"  : rawdata = self.run.GetE(time = self.time)[3]
        elif self.field == "jx"   : rawdata = self.run.GetJ(time = self.time)[0]
        elif self.field == "jy"   : rawdata = self.run.GetJ(time = self.time)[1]
        elif self.field == "jz"   : rawdata = self.run.GetJ(time = self.time)[2]
        elif self.field == "jxy"  : rawdata = self.run.GetJ(time = self.time)[3]
        elif self.field == "vix"  : rawdata = self.run.GetIonVelocity(time = self.time)[0]
        elif self.field == "viy"  : rawdata = self.run.GetIonVelocity(time = self.time)[1]
        elif self.field == "viz"  : rawdata = self.run.GetIonVelocity(time = self.time)[2]
        elif self.field == "v0x"  : rawdata = self.run.GetElectronVelocity(time = self.time)[0]
        elif self.field == "v0y"  : rawdata = self.run.GetElectronVelocity(time = self.time)[1]
        elif self.field == "v0z"  : rawdata = self.run.GetElectronVelocity(time = self.time)[2]
        elif self.field == "v0xy" : rawdata = self.run.GetElectronVelocity(time = self.time)[3]
        elif self.field == "v1x"  : rawdata = self.run.GetProtonVelocity(time = self.time)[0]
        elif self.field == "v1y"  : rawdata = self.run.GetProtonVelocity(time = self.time)[1]
        elif self.field == "v1z"  : rawdata = self.run.GetProtonVelocity(time = self.time)[2]
        elif self.field == "v2x"  : rawdata = self.run.GetAlphaVelocity(time = self.time)[0]
        elif self.field == "v2y"  : rawdata = self.run.GetAlphaVelocity(time = self.time)[1]
        elif self.field == "v2z"  : rawdata = self.run.GetAlphaVelocity(time = self.time)[2]
        elif self.field == "j0x"  : rawdata = self.run.GetElectronCurrent(time = self.time)[0]
        elif self.field == "j0y"  : rawdata = self.run.GetElectronCurrent(time = self.time)[1]
        elif self.field == "j0z"  : rawdata = self.run.GetElectronCurrent(time = self.time)[2]
        elif self.field == "j1x"  : rawdata = self.run.GetProtonCurrent(time = self.time)[0]
        elif self.field == "j1y"  : rawdata = self.run.GetProtonCurrent(time = self.time)[1]
        elif self.field == "j1z"  : rawdata = self.run.GetProtonCurrent(time = self.time)[2]
        elif self.field == "j2x"  : rawdata = self.run.GetAlphaCurrent(time = self.time)[0]
        elif self.field == "j2y"  : rawdata = self.run.GetAlphaCurrent(time = self.time)[1]
        elif self.field == "j2z"  : rawdata = self.run.GetAlphaCurrent(time = self.time)[2]
        elif self.field == "idx"  : rawdata = self.run.GetVcrossB(time = self.time)[0]
        elif self.field == "idy"  : rawdata = self.run.GetVcrossB(time = self.time)[1]
        elif self.field == "idz"  : rawdata = self.run.GetVcrossB(time = self.time)[2]
        elif self.field == "idxy" : rawdata = self.run.GetVcrossB(time = self.time)[3]
        elif self.field == "hax"  : rawdata = self.run.GetHall(time = self.time)[0]
        elif self.field == "hay"  : rawdata = self.run.GetHall(time = self.time)[1]
        elif self.field == "haz"  : rawdata = self.run.GetHall(time = self.time)[2]
        elif self.field == "haxy" : rawdata = self.run.GetHall(time = self.time)[3]
        elif self.field == "gp0x" : rawdata = self.run.GetDivElectronPressure(time = self.time)[0]
        elif self.field == "gp0y" : rawdata = self.run.GetDivElectronPressure(time = self.time)[1]
        elif self.field == "gp0z" : rawdata = self.run.GetDivElectronPressure(time = self.time)[2]
        elif self.field == "gp0xy": rawdata = self.run.GetDivElectronPressure(time = self.time)[3]
        elif self.field == "gp1x" : rawdata = self.run.GetDivProtonPressure(time = self.time)[0]
        elif self.field == "gp1y" : rawdata = self.run.GetDivProtonPressure(time = self.time)[1]
        elif self.field == "gp1z" : rawdata = self.run.GetDivProtonPressure(time = self.time)[2]
        elif self.field == "gp2x" : rawdata = self.run.GetDivAlphaPressure(time = self.time)[0]
        elif self.field == "gp2y" : rawdata = self.run.GetDivAlphaPressure(time = self.time)[1]
        elif self.field == "gp2z" : rawdata = self.run.GetDivAlphaPressure(time = self.time)[2]
        elif self.field == "gpix" : rawdata = self.run.GetDivIonPressure(time = self.time)[0]
        elif self.field == "gpiy" : rawdata = self.run.GetDivIonPressure(time = self.time)[1]
        elif self.field == "gpiz" : rawdata = self.run.GetDivIonPressure(time = self.time)[2]
        elif self.field == "rex"  : rawdata = self.run.GetResistivity(time = self.time)[0]
        elif self.field == "rey"  : rawdata = self.run.GetResistivity(time = self.time)[1]
        elif self.field == "rez"  : rawdata = self.run.GetResistivity(time = self.time)[2]
        elif self.field == "hvx"  : rawdata = self.run.GetHyperViscosity(time = self.time)[0]
        elif self.field == "hvy"  : rawdata = self.run.GetHyperViscosity(time = self.time)[1]
        elif self.field == "hvz"  : rawdata = self.run.GetHyperViscosity(time = self.time)[2]
        elif self.field == "p0xx" : rawdata = self.run.GetFullElectronPressure(time = self.time)[0]
        elif self.field == "p0xy" : rawdata = self.run.GetFullElectronPressure(time = self.time)[1]
        elif self.field == "p0xz" : rawdata = self.run.GetFullElectronPressure(time = self.time)[2]
        elif self.field == "p0yy" : rawdata = self.run.GetFullElectronPressure(time = self.time)[3]
        elif self.field == "p0yz" : rawdata = self.run.GetFullElectronPressure(time = self.time)[4]
        elif self.field == "p0zz" : rawdata = self.run.GetFullElectronPressure(time = self.time)[5]
        elif self.field == "p1xx" : rawdata = self.run.GetFullprotonPressure(time = self.time)[0]
        elif self.field == "p1xy" : rawdata = self.run.GetFullprotonPressure(time = self.time)[1]
        elif self.field == "p1xz" : rawdata = self.run.GetFullprotonPressure(time = self.time)[2]
        elif self.field == "p1yy" : rawdata = self.run.GetFullprotonPressure(time = self.time)[3]
        elif self.field == "p1yz" : rawdata = self.run.GetFullprotonPressure(time = self.time)[4]
        elif self.field == "p1zz" : rawdata = self.run.GetFullprotonPressure(time = self.time)[5]
        elif self.field == "p2xx" : rawdata = self.run.GetFullIonPressure(time = self.time)[0]
        elif self.field == "p2xy" : rawdata = self.run.GetFullIonPressure(time = self.time)[1]
        elif self.field == "p2xz" : rawdata = self.run.GetFullIonPressure(time = self.time)[2]
        elif self.field == "p2yy" : rawdata = self.run.GetFullIonPressure(time = self.time)[3]
        elif self.field == "p2yz" : rawdata = self.run.GetFullIonPressure(time = self.time)[4]
        elif self.field == "p2zz" : rawdata = self.run.GetFullIonPressure(time = self.time)[5]
        elif self.field == "pixx" : rawdata = self.run.GetFullIonPressure(time = self.time)[0]
        elif self.field == "pixy" : rawdata = self.run.GetFullIonPressure(time = self.time)[1]
        elif self.field == "pixz" : rawdata = self.run.GetFullIonPressure(time = self.time)[2]
        elif self.field == "piyy" : rawdata = self.run.GetFullIonPressure(time = self.time)[3]
        elif self.field == "piyz" : rawdata = self.run.GetFullIonPressure(time = self.time)[4]
        elif self.field == "pizz" : rawdata = self.run.GetFullIonPressure(time = self.time)[5]
        elif self.field == "ncez" : rawdata = self.run.GetNCEZ(time = self.time)[0]
        elif self.field == "jcbx" : rawdata = self.run.GetJcrossB(time = self.time)[0]
        elif self.field == "jcby" : rawdata = self.run.GetJcrossB(time = self.time)[1]
        elif self.field == "jcbz" : rawdata = self.run.GetJcrossB(time = self.time)[2]
        elif self.field == "jxby" : rawdata = self.run.GetJcrossB(time = self.time)[3]
        elif self.field == "jxbz" : rawdata = self.run.GetJcrossB(time = self.time)[4]
        elif self.field == "jybx" : rawdata = self.run.GetJcrossB(time = self.time)[5]
        elif self.field == "jybz" : rawdata = self.run.GetJcrossB(time = self.time)[6]
        elif self.field == "jzbx" : rawdata = self.run.GetJcrossB(time = self.time)[7]
        elif self.field == "jzby" : rawdata = self.run.GetJcrossB(time = self.time)[8]
        else : print "invalid field..."

        # .. set xdata & ydata
        if self.direction == "x" :
          self.xdata = np.linspace(0+self.shifts[0], self.run.l[0]+self.shifts[0], self.run.n[0]+1)
          self.ydata = rawdata[:, self.idx[0], self.idx[1]]
        elif self.direction == "y" :
          self.xdata = np.linspace(0+self.shifts[0], self.run.l[1]+self.shifts[0], self.run.n[1]+1)
          self.ydata = rawdata[self.idx[0], :, self.idx[1]]
        elif self.direction == "z" :
          self.xdata = np.linspace(0+self.shifts[0], self.run.l[2]+self.shifts[0], self.run.n[2]+1)
          self.ydata = rawdata[self.idx[0], self.idx[1], :]
        else : print "still wrong direction"

        # .. set extent if not imposed
        if self.extent == None : self.extent = [[np.min(self.xdata), np.max(self.xdata)], [np.min(self.ydata), np.max(self.ydata)]]
        else : self.extent = [[self.extent[0][0]+self.shifts[0], self.extent[0][1]+self.shifts[0]], [self.extent[1][0]+self.shifts[1], self.extent[1][1]+self.shifts[1]]]


    def MakeLabels(self):

        """
        set the labels for "field" @ "time" : "path"
        """

        # .. set the x & y labels for 2d
        if self.plane is not None:

          if self.plane == "xy" :
              self.xlabel = r'$x \, \Omega_p/c$'
              self.ylabel = r'$y \, \Omega_p/c$'

          elif self.plane == "yx" :
              self.xlabel = r'$y \, \Omega_p/c$'
              self.ylabel = r'$x \, \Omega_p/c$'

          elif self.plane == "xz" :
              self.xlabel = r'$x \, \Omega_p/c$'
              self.ylabel = r'$z \, \Omega_p/c$'

          elif self.plane == "zx" :
              self.xlabel = r'$z \, \Omega_p/c$'
              self.ylabel = r'$x \, \Omega_p/c$'

          elif self.plane == "yz" :
              self.xlabel = r'$y \, \Omega_p/c$'
              self.ylabel = r'$z \, \Omega_p/c$'

          elif self.plane == "zy" :
              self.xlabel = r'$z \, \Omega_p/c$'
              self.ylabel = r'$y \, \Omega_p/c$'
          else : print "wrong plane..."

        # .. set the x & y labels for 1d
        if self.direction is not None:

          if   self.direction == "x" : self.xlabel = r'$x \, \Omega_p/c$'
          elif self.direction == "y" : self.xlabel = r'$y \, \Omega_p/c$'
          elif self.direction == "z" : self.xlabel = r'$z \, \Omega_p/c$'
          else : print "wrong direction..."

          if self.ylabel is None : self.ylabel = ''

