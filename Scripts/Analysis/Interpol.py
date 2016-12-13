#===============================================================================
# Python library that uses SplineInter by B. Avez to generate a 2D energy 
# surface from MOCCa output files. 
#===============================================================================
import numpy as np
import matplotlib.pyplot as plt 
import os


def Inter(datafile, points, SYM=1, execloc='$HOME/Documents/Codes/inter/SplineInter/exe/'):
       
       os.system('cp %s/SplineInter.exe inter.exe'%execloc)
       num_lines = sum(1 for line in open(datafile))

       options=""
       if(SYM == 1) :
        options = "+++"
       os.system('./inter.exe 2 %d %d "tps" %s  < %s > inter.out'%(num_lines, points, options, datafile))
       
       data      = np.loadtxt(datafile)
       X = np.loadtxt('data.x', skiprows=1)
       Y = np.loadtxt('data.y', skiprows=1)
       Z = np.loadtxt('data.z', skiprows=1)
       
       os.system('rm data.x')
       os.system('rm data.y')
       os.system('rm data.z')
       os.system('rm inter.exe')
       return (data, X,Y,Z)
       
def SurfPlot( data, X,Y,Z, AXIS=None, LEVELS=[], PLOTDATA=-1, SYMX=1, SYMY=1):
    
    #Default to current axis
    if AXIS is None:
        AXIS = plt.gca()
    
    minZ     = np.min(Z)
    coordmin = c=(np.argmin(Z)/len(Z[0]),np.argmin(Z)%len(Z[0]))
    xmin = X[coordmin]
    ymin = Y[coordmin]
    
    if(SYMX == 1 and xmin < 0):
        xmin = -xmin
    if(SYMY == 1 and ymin < 0):
        ymin = -ymin
    print xmin, ymin, coordmin
    #Plot on the axis
    if(len(LEVELS)>0):
        contour = AXIS.contourf(X,Y,Z - np.min(Z), levels=LEVELS)
    else:
        contour = AXIS.contourf(X,Y,Z - np.min(Z))
        
    AXIS.plot(xmin,ymin,'kd')    
    AXIS.set_xlim(min(data[:,0]),max(data[:,0]))
    AXIS.set_ylim(min(data[:,1]),max(data[:,1]))
    cbar    = plt.colorbar(contour)
    
    if(PLOTDATA> 0):
        plt.plot(data[:,0],data[:,1], 'kx')
    

