


import numpy as np
from scipy.ndimage import gaussian_filter1d as gf



#----------------------------------------------------------
#----------------------------------------------------------
def deriv(scalarfield, dl, order=1, axis=1):
    """
    @todo: Brief Docstring for deriv
    Longer description here

    @param scalarfield @todo
    @param order @todo
    @param *dl @todo

    @return: @todo

    Exemple  :

    Creation : 2015-03-02

    """
    if len(scalarfield.shape) == 1:
        ret = deriv1D(scalarfield, order, dl)

    elif len(scalarfield.shape) == 2:
        if axis > 1:
            raise ValueError("wrong axis (%d) for array len(shape) (%d)" \
                    % (axis,len(scalarfield.shape)))
        ret = deriv2D(scalarfield, order, axis, dl)

    elif len(scalarfield.shape) == 3:
        if axis > 2:
            raise ValueError("wrong axis (%d) for array len(shape) (%d)" \
                    % (axis,len(scalarfield.shape)))
        ret = deriv2D(scalarfield, order, axis, dl)

    else:
        raise ValueError("bad shape of the scalar field")

    return ret
#==========================================================





#----------------------------------------------------------
#----------------------------------------------------------
def deriv1D(scalarfield, order, dl):
    """
    @todo: Brief Docstring for deriv1D
    Longer description here

    @param scalarfield @todo
    @param order @todo
    @param dl @todo

    @return: @todo

    Exemple  :

    Creation : 2015-03-02

    """
    return gf(scalarfield, order=order,sigma=1, axis=0)/dl**order
#==========================================================








#----------------------------------------------------------
#----------------------------------------------------------
def deriv2D(scalarfield, order, axis, dl):
    """
    @todo: Brief Docstring for deriv1d
    Longer description here

    @param scalarfield @todo
    @param order @todo
    @param dl @todo

    @return: @todo

    Exemple  :

    Creation : 2015-03-02

    """
    return gf(scalarfield, order=order,sigma=1, axis=axis)/dl**order
#==========================================================







#----------------------------------------------------------
#----------------------------------------------------------
def deriv3D(scalarfield, order, axis, dl):
    """
    @todo: Brief Docstring for deriv1d
    Longer description here

    @param scalarfield @todo
    @param order @todo
    @param dl @todo

    @return: @todo

    Exemple  :

    Creation : 2015-03-02

    """
    return gf(scalarfield, order=order,sigma=1, axis=axis)/dl**order
#==========================================================







#----------------------------------------------------------
#----------------------------------------------------------
def div(vectorfield, dl):
    """
    @todo: Brief Docstring for div
    Longer description here

    @param vectorfield @todo

    @return: @todo

    Exemple  :

    Creation : 2015-03-02

    """
    ndims = len(vectorfield.shape)-1

    if ndims == 1:
        return deriv(vectorfield[...,0], dl[0], 1, 0)

    elif ndims == 2:
        ret = deriv(vectorfield[...,0], dl[0], 1, 0) \
            + deriv(vectorfield[...,1], dl[1], 1, 1)

    elif ndims == 3:
        ret = deriv(vectorfield[...,0], dl[0], 1, 0) \
            + deriv(vectorfield[...,1], dl[1], 1, 1) \
            + deriv(vectorfield[...,2], dl[2], 1, 2)


    return ret
#==========================================================






#----------------------------------------------------------
#----------------------------------------------------------
def curl(vectorfield, dl):
    """
    @todo: Brief Docstring for curl
    Longer description here

    @param vectorfield @todo
    @param *dl @todo

    @return: @todo

    Exemple  :

    Creation : 2015-03-02

    """

    ndims = len(vectorfield.shape)-1
    shape = vectorfield.shape
    dtype = vectorfield.dtype

    res = np.zeros(shape=shape, dtype=dtype)

    if ndims == 1:
        res[...,0] = 0
        res[...,1] = -deriv(vectorfield[...,2], dl, 1, 0)
        res[...,2] =  deriv(vectorfield[...,1], dl, 1, 0)


    elif ndims == 2:
        res[...,0] =  deriv(vectorfield[...,2], dl[1], 1, 1)

        res[...,1] = -deriv(vectorfield[...,2], dl[0], 1, 0)

        res[...,2] =  deriv(vectorfield[...,1], dl[0], 1, 0) \
                     -deriv(vectorfield[...,0], dl[1], 1, 1)


    elif ndims == 3:
        res[...,0] =  deriv(vectorfield[...,2], dl[1], 1, 1) \
                     -deriv(vectorfield[...,1], dl[2], 1, 2)

        res[...,1] = -deriv(vectorfield[...,2], dl[0], 1, 0) \
                     +deriv(vectorfield[...,0], dl[2], 1, 2)

        res[...,2] =  deriv(vectorfield[...,1], dl[0], 1, 0) \
                     -deriv(vectorfield[...,0], dl[1], 1, 1)
    return res
#==========================================================








#----------------------------------------------------------
#----------------------------------------------------------
def scallaplacian(scalarfield, dl):
    """
    @todo: Brief Docstring for scallaplacian
    Longer description here

    @param scalarfield @todo
    @param *dl) @todo

    @return: @todo

    Exemple  :

    Creation : 2015-03-02

    """
    ndims = len(scalarfield.shape)

    if ndims == 1:
        ret = deriv(scalarfield, dl, 2, 0)

    elif ndims == 2:
        ret = deriv(scalarfield, dl[0], 2, 0) \
            + deriv(scalarfield, dl[1], 2, 1) \

    elif ndims == 3:
        ret = deriv(scalarfield, dl[0], 2, 0) \
            + deriv(scalarfield, dl[1], 2, 1) \
            + deriv(scalarfield, dl[2], 2, 2)

    else:
        raise ValueError("dimension (%d) is to big" % ndims)

    return ret
#==========================================================







#----------------------------------------------------------
#----------------------------------------------------------
def vectlaplacian(vectorfield, dl):
    """
    @todo: Brief Docstring for vectlaplacian
    Longer description here

    @param vectorfield @todo
    @param *dl @todo

    @return: @todo

    Exemple  :

    Creation : 2015-03-02

    """
    ndims = len(vectorfield.shape) - 1
    shape = vectorfield.shape
    dtype = vectorfield.dtype

    ret   = np.zeros(shape=shape, dtype=dtype)

    if ndims == 1:
        ret[...,0] = scallaplacian(vectorfield[...,0], dl)
        ret[...,1] = scallaplacian(vectorfield[...,1], dl)
        ret[...,2] = scallaplacian(vectorfield[...,2], dl)


    elif ndims == 2:
        ret[...,0] = scallaplacian(vectorfield[...,0], dl)
        ret[...,1] = scallaplacian(vectorfield[...,1], dl)
        ret[...,2] = scallaplacian(vectorfield[...,2], dl)


    elif ndims == 3:
        ret[...,0] = scallaplacian(vectorfield[...,0], dl)
        ret[...,1] = scallaplacian(vectorfield[...,1], dl)
        ret[...,2] = scallaplacian(vectorfield[...,2], dl)

    return ret

#==========================================================





#----------------------------------------------------------
#----------------------------------------------------------
def dot(V1,V2):
    """@todo: Brief Docstring for dot
    Longer description here

    @param V1 @todo
    @param V2 @todo

    @return: @todo

    Exemple  :

    Creation: 2015-03-03

    """
    return V1[...,0]*V2[...,0] + V1[...,1]*V2[...,1] + V1[...,2]*V2[...,2]
#==========================================================







#----------------------------------------------------------
#----------------------------------------------------------
def VgradV(V, dl):
    """
    @todo: Brief Docstring for VgradV
    Longer description here

    @param V @todo
    @param dl @todo

    @return: @todo

    Exemple  :

    Creation : 2015-03-03

    """
    ndims = len(V.shape)-1
    shape = V.shape
    dtype = V.dtype

    vgv  = np.zeros(shape=shape, dtype=dtype)

    if ndims == 1:
        for c in np.arange(3):
            vgv[...,c] = V[...,0]*deriv(V[...,c], dl[0], order=1, axis=0)

    elif ndims == 2:
        for c in np.arange(3):
            vgv[...,c] = V[...,0]*deriv(V[...,c], dl[0], order=1, axis=0) \
                       + V[...,1]*deriv(V[...,c], dl[1], order=1, axis=1)

    elif ndims == 3:
        for c in np.arange(3):
            vgv[...,c] = V[...,0]*deriv(V[...,c], dl[0], order=1, axis=0) \
                       + V[...,1]*deriv(V[...,c], dl[1], order=1, axis=1) \
                       + V[...,2]*deriv(V[...,c], dl[2], order=1, axis=2)


#==========================================================






#----------------------------------------------------------
#----------------------------------------------------------
def divT(T,dl):
    """
    @todo: Brief Docstring for divT
    Longer description here

    @param T @todo
    @param dl @todo

    @return: @todo

    Exemple  :

    Creation : 2015-03-03

    """

    ndims = len(T.shape)-2
    divt  = np.zeros(shape=T.shape[:-1], dtype=T.dtype)

    if ndims == 1:
        divt[...,0] = deriv(T[...,0,0], dl[0], order=1, axis=0)
        divt[...,1] = deriv(T[...,0,1], dl[0], order=1, axis=0)
        divt[...,2] = deriv(T[...,0,2], dl[0], order=1, axis=0)

    elif ndims == 2:

        divt[...,0] = deriv(T[...,0,0], dl[0], order=1, axis=0) \
                    + deriv(T[...,0,1], dl[1], order=1, axis=1)

        divt[...,1] = deriv(T[...,0,1], dl[0], order=1, axis=0) \
                    + deriv(T[...,1,1], dl[1], order=1, axis=1)

        divt[...,2] = deriv(T[...,0,2], dl[0], order=1, axis=0) \
                    + deriv(T[...,1,2], dl[1], order=1, axis=1)

    elif ndims == 3:

        divt[...,0] = deriv(T[...,0,0], dl[0], order=1, axis=0) \
                    + deriv(T[...,0,1], dl[1], order=1, axis=1) \
                    + deriv(T[...,0,2], dl[2], order=1, axis=2)

        divt[...,1] = deriv(T[...,0,1], dl[0], order=1, axis=0) \
                    + deriv(T[...,1,1], dl[1], order=1, axis=1) \
                    + deriv(T[...,1,2], dl[2], order=1, axis=2)

        divt[...,2] = deriv(T[...,0,2], dl[0], order=1, axis=0) \
                    + deriv(T[...,1,2], dl[1], order=1, axis=1) \
                    + deriv(T[...,2,2], dl[2], order=1, axis=2)

    return divt

#==========================================================







#----------------------------------------------------------
#----------------------------------------------------------
def trace(T):
    """
    calculates the Trace of a tensor field

    @param T tensor
    @return: trace of the tensor field

    Exemple  : tr = fields.trace(Pe)

    Creation : 2015-03-05
    """
    return T[...,0] + T[...,1] + T[...,2]
#==========================================================




#----------------------------------------------------------
#----------------------------------------------------------
def norm(V):
    """
    calculates the norm of a vector field

    @param V the vector field

    @return: sqrt(Vx**2 + Vy**2 + Vz**2)

    Exemple  :

    Creation : 2015-03-05

    """
    return np.sqrt(V[...,0]**2 + V[...,1]**2 + V[...,2]**2)
#==========================================================







