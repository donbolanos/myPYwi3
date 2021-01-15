import numpy as np

# this module read the potental Az store in folder Az and does stuff with it
# REMARQUE SUR LES VARIABLES : 
# --> position positionRef, un couple de coordonnees [x,y]
#     si rien, on prend le milieu
# --> path : chemin jusqu'au dossier ou l'on contient tous les fichiers 'Az_t_???.txt'

def readAz(path,time,Nx,Ny,position = None, positionRef = None):
    fname = path + 'Az_t_' + str(time) +'.txt'
    data = np.loadtxt(fname , skiprows = 1)

    x  = data[:,0]
    y  = data[:,1]
    Az = data[:,2]

    if len(position) != 2 :
        position = np.zeros(2)
        position[0] = (Nx-1)/2
        position[1] = (Ny-1)/2

    if len(positionRef) != 2 :
        positionRef = np.zeros(2)
        positionRef[0] = 1
        positionRef[1] = (Ny-1)/2

    far = Az[positionRef[0]*Ny + positionRef[1]]
    az  = Az[position[0]*Ny + position[1]]

    return far - az

def getAzTemp(path,timegroup,Nx,Ny,position = None, positionRef = None):

    if len(position) != 2 :
        position = np.zeros(2)
        position[0] = (Nx-1)/2
        position[1] = (Ny-1)/2

    if len(positionRef) != 2 :
        positionRef = np.zeros(2)
        positionRef[0] = 1
        positionRef[1] = (Ny-1)/2

    az = []

    for time in timegroup : 
        tmp = readAz(path,time,Nx,Ny,position,positionRef)
        az.append(tmp)

    az = np.array(az)
    
    return az

def getEzTemp(path,timegroup,Nx,Ny,position = None, positionRef = None, returnPot = False):

    if len(position) != 2 :
        position = np.zeros(2)
        position[0] = (Nx-1)/2
        position[1] = (Ny-1)/2

    if len(positionRef) != 2 :
        positionRef = np.zeros(2)
        positionRef[0] = 1
        positionRef[1] = (Ny-1)/2

    az = getAzTemp(path,timegroup,Nx,Ny,position,positionRef)
    #SORTING 
    t,a = zip(*sorted(zip(timegroup,az)))
    dt = np.array(t[1:]) - np.array(t[0:-1])    

    ez = (np.array(a[1:]) - np.array(a[:-1]))/dt
    if returnPot :
        return t[:-1]+dt/2,ez,a
    elif returnPot == False :
        return t[:-1]+dt/2,ez
    else :
        print('What are you doing, motherfucker ?')