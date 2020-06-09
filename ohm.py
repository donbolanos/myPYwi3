
import collections
from scipy.ndimage import gaussian_filter1d as gf
import numpy as np
import scipy.ndimage as ndimage
#import pressure
import sys
#hekle_path = '/home/sbolanos/Documents/code/'
#sys.path.append(hekle_path+"/mypywi")
#sys.path.append(hekle_path+"/mypywi/fields")
import fields.fields as fields





#----------------------------------------------------------
#----------------------------------------------------------
def hyperresistivity(run, time, smooth='no', sigma=1):
    """returns the hyper-resistive term
    @return: numpy array with hyperresistive term
    Exemple  :
    Creation : 2012-08-15 14:52:15.062321
    """

    # checks whether 'time' is one scalar or a time interval
    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    hprsty = np.zeros(shape=run.vfield_shape, dtype=run.dtype)
    nu     = -run.GetHyperResistivity()

    for t in time:
        J = run.GetJ(t)
        LapJ = fields.vectlaplacian(J, run.dl)
        hprsty += nu * LapJ

    hprsty /= ntot

    if smooth.lower() == 'yes':
        for c in range(3):
            hprsty[...,c] = gf(hprsty[...,c],sigma=sigma,order=0)

    return hprsty
#==========================================================







#----------------------------------------------------------
#----------------------------------------------------------
def electron_inertia(run, time, smooth='no', sigma=1):
    """returns the electron inertia term
    @return: numpy array with electron inertia term
    Exemple  :
    Creation : 2012-08-15 14:52:15.062321
    """
    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    einertia = np.zeros(shape=run.vfield_shape, dtype=run.dtype)

    for t in time:
        Ve  = run.GetVe(time)
        vgv = fields.VgradV(Ve, run.dl)
        q   = run.GetCharge('electrons')
        m_e = run.GetMass('electrons')

        einertia += m_e/q * vgv


    einertia /= ntot

    if smooth.lower() == 'yes':
        for c in range(3):
            einertia[...,c] = gf(inertia[...,c], sigma=sigma, order=0)

    return einertia
#==========================================================










#----------------------------------------------------------
#----------------------------------------------------------
def resistivity(run, time, smooth='no',sigma=1):
    """@todo: returns the resistive term
    @return: 2D numpy array with etaJ
    Exemple  :
    Creation : 2012-08-15 14:47:43.852335
    """


    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    etaJ = np.zeros(shape=run.vfield_shape, dtype=run.dtype)

    for t in time:
        etaJ   += run.GetJ(t)

    etaJ /= ntot

    if smooth.lower() == 'yes':
        for c in range(3):
            etaJ[...,c] = gf(etaJ[...,c], sigma=sigma, order=0)

    etaJ *= run.GetResistivity()

    return etaJ
#==========================================================







#----------------------------------------------------------
#----------------------------------------------------------
def ion_ideal(run, time, smooth='no', sigma = 1):
    """calculates the ion ideal electric field
    @return: numpy array 2D containing -VixB
    Exemple  :
    Creation : 2012-08-15 09:48:31.278583
    """


    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    mVixB= np.zeros(shape=run.vfield_shape, dtype=run.dtype)

    for t in time:
        Vi = run.GetVi(t)
        B  = run.GetB(t)
        mVixB += -np.cross(Vi,B)

    mVixB /= ntot

    if smooth.lower() == 'yes':
        for c in range(3):
            mVixB[...,c] = gf(mVixB[...,c], sigma=sigma, order=0)

    return -np.cross(Vi,B)

#==========================================================







#----------------------------------------------------------
#----------------------------------------------------------
def hall(run, time, smooth='no', sigma=1,silent='no'):
    """calculates the hall term 

    @return: numpy 2D array containing the Hall term

    Exemple  : 

    Creation : 2012-08-15 10:04:33.349386

    """

    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    ht  = np.zeros(shape=run.vfield_shape, dtype=run.dtype)


    for t in time:
        J = run.GetJ(t) # silent=silent implemented initially, it seems useless
        B = run.GetB(t) # silent=silent implemented initially, it seems useless
        e = run.GetCharge('proton')
        n = run.GetN(t, 'electron')

        for i in range(3):
            ht[..., i] += np.cross(J,B)[..., i]/(n*e)

    ht /= ntot

    if smooth.lower() == 'yes':
        for c in range(3):
            ht[...,c] = gf(ht[...,c], sigma=sigma, order=0)

    return ht
#==========================================================









#----------------------------------------------------------
#----------------------------------------------------------
def electron_ideal(run, time, smooth='no', sigma=1):
    """calculates the electron ideal electric field
    @return: numpy 2D array containing the Hall term
    Exemple  :
    Creation : 2012-08-15 10:04:33.349386
    """

    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    ei   = np.zeros(shape=run.vfield_shape, dtype=run.dtype)

    for t in time:
        Ve = run.GetVe(t)
        B  = run.GetB(t)
        ei = -np.cross(Ve,B)

    ei /= ntot

    if smooth.lower() == 'yes':
        for c in range(3):
            ei[...,c] = gf(ei[...,c], sigma=sigma, order=0)

    return ei
#==========================================================








#----------------------------------------------------------
#----------------------------------------------------------
def electron_pressure(run, time, smooth='no', sigma=1):
    """calculates the electron pressure term
    this method is private
    @param silent 'no' or 'yes'
    @return: 2D numpy array containing the electron pressure term
    Exemple  :
    Creation : 2012-08-15 10:21:33.373750
    """

    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    mdPeone = np.zeros(shape=run.vfield_shape, dtype=run.dtype)


    for t in time:
        Pe    = run.GetPe(t)
        divPe = fields.divT(Pe, run.dl)
        n     = run.GetNe(t)

        for c in np.arange(3):
	    mdPeone[...,c] += divPe[...,c]/n

    mdPeone /= ntot

    mdPeone /= run.GetCharge('electrons')

    if smooth.lower() == 'yes':
        for c in range(3):
            mdPeone[...,c] = gf(mdPeone[...,c], sigma=sigma, order=0)


    return mdPeone
#==========================================================









#----------------------------------------------------------
# #----------------------------------------------------------
# def electron_pressure_gyro(run, time, smooth='no', sigma=1, silent='no'):
#     """calculates the contribution of the gyrotropic electron pressure
#     to the electron pressure term in Ohm's law


#     @param silent 'no' or 'yes'

#     @return: 2D numpy array containing the electron pressure term

#     Exemple  : 

#     Creation : 2012-08-15 10:21:33.373750

#     """

#     if isinstance(time, collections.Iterable) == False:
#         time = [time]
#         ntot = 1
#     else:
#         ntot = time.size

#     shape = (3,run.ncells[0]+1, run.ncells[1]+1)

#     mdPeone = np.zeros(shape,'float32',order='F')

#     for t in time:
#         Pe    = run.GetPe(t, silent=silent)
#         B     = run.GetB(t, silent=silent)
#         P_ng,P_pp,P_para,P_perp = pressure.nongyro_decomposition(Pe,B)

#         divPeng = run.divT(P_pp)
#         n     = run.GetN(t, species='electrons', silent=silent)

#         mdPeone += divPeng/n

#     mdPeone /= ntot

#     mdPeone /= run.GetCharge('electrons')

#     if smooth.lower() == 'yes':
#         for c in range(3):
#          mdPeone[c,:,:] = ndimage.gaussian_filter(mdPeone[c,:,:],\
#                             sigma=sigma, order=0)


#     return mdPeone
# #==========================================================






# #----------------------------------------------------------
# #----------------------------------------------------------
# def electron_pressure_nongyro(run, time, smooth='no', sigma=1):
#     """calculates the contribution of the non-gyrotropic electron pressure
#     to the electron pressure term in Ohm's law


#     @param silent 'no' or 'yes'

#     @return: 2D numpy array containing the electron pressure term

#     Exemple  :

#     Creation : 2012-08-15 10:21:33.373750

#     """

#     if isinstance(time, collections.Iterable) == False:
#         time = [time]
#         ntot = 1
#     else:
#         ntot = time.size

#     mdPeone = np.zeros(shape=run.vfield_shape,dtype=run.dtype)

#     for t in time:
#         Pe    = run.GetPe(t)
#         B     = run.GetB(t)
#         P_ng,P_pp,P_para,P_perp = pressure.nongyro_decomposition(Pe,B)


#         divPeng = run.divT(P_ng)
#         n     = run.GetN(t, species='electrons', silent=silent)

#         mdPeone += divPeng/n

#     mdPeone /= ntot

#     mdPeone /= run.GetCharge('electrons')

#     if smooth.lower() == 'yes':
#         for c in range(3):
#          mdPeone[c,:,:] = ndimage.gaussian_filter(mdPeone[c,:,:],\
#                             sigma=sigma, order=0)


#     return mdPeone
# #==========================================================


















