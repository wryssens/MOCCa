#===============================================================================
# Python library that uses SplineInter by B. Avez to generate a 2D energy 
# surface from MOCCa output files. 
#===============================================================================
import numpy as np
import matplotlib.pyplot as plt 
import os
import glob

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
       
       os.system('rm data.x')
       os.system('rm data.y')
       os.system('rm data.z')
       os.system('rm inter.exe')
       return (data, X,Y,Z)
       
def SurfPlot( data, X,Y,Z, AXIS=None, LEVELS=[], PLOTDATA=-1, SYMX=1, SYMY=1, LABEL='', CMAP=plt.cm.jet_r, NORM=1, PLOTMIN=1):

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

    print 'surface minimum at ', xmin, ymin
    if(PLOTMIN == 1):
        AXIS.plot(xmin,ymin,'ks')

    if(PLOTDATA> 0):
        AXIS.plot(data[:,0],data[:,1], 'kx')

    return(xmin,ymin,np.min(Z), contour)

#===============================================================================
def BetaGammaFULL(prefix, A, title='', AX=plt.gca(), PLOTDATA=0, LEVELS=None,
                  CMAP=plt.cm.jet_r,
                  execloc  ='$HOME/Documents/Codes/inter/SplineInter/exe/',
                  scriptloc='/home/wryssens/Documents/Codes/MOCCa/Scripts/Plotting/BGscripts'):
    #---------------------------------------------------------------------------
    # Use matplotlib plot something in the full beta-gamma plane
    #---------------------------------------------------------------------------
    # Get the data in terms of Q20, Q22
    data  = np.loadtxt(prefix + '.t.qlm.tab', skiprows=1)
    Edata = np.loadtxt(prefix + '.e.tab', skiprows=1)
    
    N     = len(data[:,1])    
    with open('input', 'w') as f:
     for i in range(N):
      Q20 = data[i,1]
      Q22 = data[i,3]
      
      f.write('%10.6e %10.6e %10.6e \n'%(Q20,Q22,Edata[i,1]))
      
      ##########################################################################
      # Reflection here
      Q     = np.sqrt(Q20**2 + 2*Q22**2)
      gamma = np.arctan2(np.sqrt(2)*Q22,Q20)
           
      # 0 -> -pi/3
      if(abs(gamma * 180/np.pi) > 1  and Q > 0.1):
        f.write('%10.6e %10.6e %10.6e \n'%(Q20,-Q22,Edata[i,1]))
        
      refQ20 = Q * np.cos(2*np.pi/3 - gamma)
      refQ22 = 1.0/np.sqrt(2) * Q * np.sin(2*np.pi/3 - gamma)
      if(Q > 0.1):
        f.write('%10.6e %10.6e %10.6e \n'%(refQ20,refQ22,Edata[i,1]))
        
      refQ20 = Q * np.cos(2*np.pi/3 + gamma)
      refQ22 = 1.0/np.sqrt(2) * Q * np.sin(2*np.pi/3 + gamma)
      if(abs(gamma - np.pi/3) > 0.1 and abs(gamma) > 0.1 and Q > 0.1):
        f.write('%10.6e %10.6e %10.6e \n'%(refQ20,refQ22,Edata[i,1]))
        
        
      refQ20 = Q * np.cos(- 2*np.pi/3 + gamma)
      refQ22 = 1.0/np.sqrt(2) * Q * np.sin(- 2*np.pi/3 + gamma)
      if ( Q > 0.1):
        f.write('%10.6e %10.6e %10.6e \n'%(refQ20,refQ22,Edata[i,1]))
        
      refQ20 = Q * np.cos(- 2*np.pi/3 - gamma)
      refQ22 = 1.0/np.sqrt(2) * Q * np.sin(- 2*np.pi/3 - gamma)
      if( abs(gamma) > 0.1 and Q > 0.1):
        f.write('%10.6e %10.6e %10.6e \n'%(refQ20,refQ22,Edata[i,1]))
    
    ############################################################################  
    # Interpolate the data in iq1 and iq2
    ntps=0
    with open('input', 'r') as f:
        for i in f:
            ntps = ntps+1
    os.system('cp %s/SplineInter.exe .'%execloc)
    os.system('./SplineInter.exe 2 %d %d "tps"   < %s > inter.out'%(ntps,201,'input'))

    X = np.loadtxt('data.x', skiprows=1)
    Y = np.loadtxt('data.y', skiprows=1)
    Z = np.loadtxt('data.z', skiprows=1)
    #---------------------------------------------------------------------------
    # Retransform the data to Q0 and gamma
    R0    = 1.2 * A**(1.0/3.0)
    Q     = 4.0*np.pi/(3.0*R0**2*A) * np.sqrt((X**2+2*Y**2))
    gamma = np.arctan2(np.sqrt(2)*Y, X) 

    points = np.loadtxt('input')

    Qdata     = np.sqrt((points[:,0]**2 + 2*points[:,1]**2))
    Qdata     = Qdata * 4.0 * np.pi/(3.0 * R0**2 * A ) 
    gammadata = np.arctan2(np.sqrt(2)*points[:,1],points[:,0]) 

    #---------------------------------------------------------------------------
    # Plot the contour
    if(LEVELS != None):    
        contour = AX.contour(gamma, Q, Z - np.amin(Z), levels=LEVELS,cmap=CMAP)
    else:   
        contour = AX.contour(gamma, Q, Z - np.amin(Z), cmap=CMAP)
    if(PLOTDATA == 1):
        AX.plot( gammadata, Qdata,'ro')

    os.system('rm  input *.z data.x data.y inter.out SplineInter.exe')
    
    return(contour)

#===============================================================================
#def BetaGammaAlt(datafile, A, title='', AX=plt.gca(), PLOTDATA=0, LEVELS=None,
#                                     CMAP=plt.cm.jet_r,
#                                     execloc  ='$HOME/Documents/Codes/inter/SplineInter/exe/',
#                                     scriptloc='/home/theorie/ryssens/Documents/Codes/MOCCa/Scripts/Plotting/BGscripts'):

#    ############################################################################
#    #First, extend the data to the full beta-gamma plane
#    os.system('cp %s/surf.awk .'%scriptloc)
#    try:
#        essai = np.loadtxt(datafile)
#    except IOError:
#        print 'No datafile %s'%datafile
#        return
#    
#    os.system('awk -f surf.awk "mass=%g" < %s '%(A,datafile))
#    
#    #Combine all of the 6 files
#    os.system('cat ESPDATA1.bg ESPDATA2.bg ESPDATA3.bg ESPDATA4.bg ESPDATA5.bg ESPDATA6.bg >  input')
#    
#    # Interpolate
#    ntps=0
#    with open('input', 'r') as f:
#        for i in f:
#            ntps = ntps+1
#    
#    check = np.loadtxt('input')
#    for i in range(len(check[:,0])):
#        for j in range(i+1,len(check[:,0])):
#            if(abs(check[i,0] - check[j,0]) < 0.01 and abs(check[i,1]-check[j,1]) < 0.01):
#                print 'Duplicate', i,j
#    os.system('cp %s/SplineInter.exe .'%execloc)
#    os.system('./SplineInter.exe 2 %d %d "tps" "++" < %s > inter.out'%(ntps, 501, 'input'))
#    
#    # Transfer to (Q20,Q22)
#    data = np.loadtxt('inter.out', skiprows=25)
#    
#    Q20  = 0.5 *np.sqrt(5.0/(16.0 * np.pi)) * ( 2 * data[:,0] + data[:,1] )
#    Q22  = 1.5 * data[:,1] * np.sqrt(5.0/(4.0 * np.pi)) /(2 * np.sqrt(6.0))
#    
#    #Transfor to (Q, gamma)
#    R0    = 1.2 * A**(1.0/3.0)
#    Q     = np.sqrt((Q20**2+2*Q22**2)) #4.0*np.pi/(3.0*R0**2*A) *
#    gamma = np.arctan2(np.sqrt(2)*Q22, Q20) 
#    
#    #---------------------------------------------------------------------------
#    # Plot the contour
#    if(LEVELS != None):    
#        contour = AX.contour(gamma, Q, Z - np.amin(Z), levels=LEVELS,cmap=CMAP)
#    else:   
#        contour = AX.contour(gamma, Q, Z - np.amin(Z), cmap=CMAP)
#    if(PLOTDATA == 1):
#        AX.plot( gammadata, Qdata,'ro')

#    os.system('rm  *.z data.x data.y inter.out SplineInter.exe')
#    
#    return(contour)
#===============================================================================
def BetaGamma(datafile, A, title='', figname  = 'BG.out.eps', MAXB=0.7, EMAX = 10, 
                                     Q1MIN    = 10**8 , Q2MIN = 10**8, PLOTDATA=1, 
                                     LEGEND   = 1,
                                     execloc  ='$HOME/Documents/Codes/inter/SplineInter/exe/',
                                     scriptloc='/home/theorie/ryssens/Documents/Codes/MOCCa/Scripts/Plotting/BGscripts'):

    #First, extend the data to the full beta-gamma plane
    os.system('cp %s/surf.awk .'%scriptloc)
    try:
        essai = np.loadtxt(datafile)
    except IOError:
        print 'No datafile %s'%datafile
        return
    
    os.system('awk -f surf.awk "mass=%g" < %s '%(A,datafile))
    
    #Combine all of the 6 files
    os.system('cat ESPDATA1.bg ESPDATA2.bg ESPDATA3.bg ESPDATA4.bg ESPDATA5.bg ESPDATA6.bg >  input')
    
    # Interpolate
    ntps=0
    with open('input', 'r') as f:
        for i in f:
            ntps = ntps+1
    
    check = np.loadtxt('input')
    for i in range(len(check[:,0])):
        for j in range(i+1,len(check[:,0])):
            if(abs(check[i,0] - check[j,0]) < 0.01 and abs(check[i,1]-check[j,1]) < 0.01):
                print 'Duplicate', i,j
    os.system('cp %s/SplineInter.exe .'%execloc)
    os.system('./SplineInter.exe 2 %d %d "tps" "++" < %s > inter.out'%(ntps, 501, 'input'))
    
    try:
        data = np.loadtxt('data.z', skiprows=1)
    
        minx = 1000000000
        for i in range(len(data[:,0])):
            if data[i,0] < Q1MIN :
                if data[i,1] < Q2MIN :
                    if data[i,2] < minx:
                        minx = data[i,2]
        
        print minx, np.where(data[:,2] == minx), data[np.where(data[:,2] == minx),0]
        
        # rename generic filenames produced by the interpolation code
        os.system('cp data.z ESPDATA.input.z  ')    # data for solid  isolines at integer energies in MeV
        os.system('cp data.z ESPDATA.input.2.z')    # data for dotted isolines at  half-integer energies

        if(PLOTDATA == 1):
            PD = ''
        else:
            PD = '!'

        if(title !='') :
            NC = ''
        else :
            NC = '!'
            
        if(LEGEND == 0):
            LEG = '!'
        else:
            LEG = ''

        with open("%s/betagamma.template"%(scriptloc), 'r') as template:
            with open('BG.gle', 'w') as file: 
                for line in template:
                    line = line.replace('$NAMESTRING', title)
                    line = line.replace('$XMIN', str(minx))
                    line = line.replace('$MAXB', str(MAXB))
                    line = line.replace('$EMAX', str(EMAX))
                    line = line.replace('$NC', str(NC))
                    line = line.replace('$PD', str(PD))
                    line = line.replace('$LEGEND', str(LEG))
                    file.write(line)
                    
        os.system('cp -r %s .'%scriptloc)
        os.system('gle BG.gle')
        os.system('mv BG.eps %s'%figname)
        os.system('rm  input  *.dat *.gle *.z data.x data.y ')
        os.system('rm -r BGscripts/')
    except IOError:
        print 'Interpolation failed.'
    
    os.system('rm  *.bg  ESP.xflr6 SplineInter.exe' )
    os.system('rm surf.awk')

    
