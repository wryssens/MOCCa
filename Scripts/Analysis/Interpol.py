#===============================================================================
# Python library that uses SplineInter by B. Avez to generate a 2D energy 
# surface from MOCCa output files. 
#===============================================================================
import numpy as np
import matplotlib.pyplot as plt 
import os


def Inter(datafile, points, COLUMN=2, SYM=1, execloc='$HOME/Documents/Codes/inter/SplineInter/exe/'):
       
       if(COLUMN < 2) : 
        print 'Cannot plot a surface of the X/Y coordinates.'
        return
       elif(COLUMN != 2) :
        # Make a new file that reorders columns
        ref = np.loadtxt(datafile)
        new = ref
        new[:,2]      = ref[:,COLUMN]
        new[:,COLUMN] = ref[:,2]
        np.savetxt('tmp',new)
        toplot = 'tmp'
       else:
        toplot =datafile
       
       os.system('cp %s/SplineInter.exe inter.exe'%execloc)
       num_lines = sum(1 for line in open(datafile))
       
       options=""
       if(SYM == 1) :
        options = "+++"
       os.system('./inter.exe 2 %d %d "tps" %s  < %s > inter.out'%(num_lines, points, options, toplot))
       
       data      = np.loadtxt(datafile)
       X = np.loadtxt('data.x', skiprows=1)
       Y = np.loadtxt('data.y', skiprows=1)
       Z = np.loadtxt('data.z', skiprows=1)
       
       #os.system('rm data.x')
       #os.system('rm data.y')
       #os.system('rm data.z')
       os.system('rm inter.exe')
       return (data, X,Y,Z)
       
def SurfPlot( data, X,Y,Z, AXIS=None, LEVELS=[], PLOTDATA=-1, SYMX=1, SYMY=1, LABEL='', CMAP=plt.cm.jet_r, NORM=1):

    #Default to current axis
    if AXIS is None:
        AXIS = plt.gca()

    minZ     = np.min(Z)
    coordmin = (np.argmin(Z)/len(Z[0]),np.argmin(Z)%len(Z[0]))
    xmin = X[coordmin]
    ymin = Y[coordmin]

    if(SYMX == 1 and xmin < 0):
        xmin = -xmin
    if(SYMY == 1 and ymin < 0):
        ymin = -ymin

    if(NORM == 1):
        off = np.min(Z)
    else:
        off = 0

    #Plot on the axis
    if(len(LEVELS)>0):
        contour = AXIS.contourf(X,Y,Z - off, levels=LEVELS, cmap=CMAP)
    else:
        contour = AXIS.contourf(X,Y,Z - off, cmap=CMAP)

    AXIS.plot(xmin,ymin,'kd')
    AXIS.set_xlim(min(data[:,0]),max(data[:,0]))
    AXIS.set_ylim(min(data[:,1]),max(data[:,1]))
    cbar    = plt.colorbar(contour, label=LABEL)

    if(PLOTDATA> 0):
        plt.plot(data[:,0],data[:,1], 'kx')

    return(xmin,ymin,np.min(Z))
