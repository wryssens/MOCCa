#-------------------------------------------------------------------------------
# Python module for plotting data generated by MOCCa and the analysis scripts
#    MOCCa.spectra.awk
#    AxialSort.awk
#
# Matplotlib is our plotting engine of choice
#
# Files from a set of runs should be organized as 
#     PREFIX        e.tab
#     PREFIX       ef.tab
#     PREFIX         .tab
#
# etc...
#-------------------------------------------------------------------------------

import numpy as np
from scipy.interpolate import interp1d
from tokenizer import *
import matplotlib.pyplot as plt
import glob
import sys
import math

def Determinedata(PREFIX, XARG, YARG,PC,SC, SORT):
    
    altx = 0
    alty = 0
    absy = 0
    fac  = 1
    divy = 0
    if(XARG=='Q20') :
        xlabel =r'$\langle \hat{Q}_{20} \rangle$ (fm$^2$)'
        xfname =PREFIX + '.t.qlm.tab'
        xcolumn=1
        if(PC != 1) :
            xcolumn = 3
    elif(XARG=='B20') :
        xlabel =r'$\beta_{20}$ '
        xfname =PREFIX + '.t.qlm.tab'
        xcolumn=2
        if(PC != 1) :
            xcolumn = 4
        if(SC != 1) :
            xcolumn = 4
        if(SC !=1 and PC != 1):
            xcolumn = 6
    elif(XARG=='B30') :
        if(PC == 1) :
            print "Can't plot Q30 when parity is conserved."
            sys.exit()
        xlabel =r'$\beta_{30}$ '
        xfname =PREFIX + '.t.qlm.tab'
        xcolumn=8
        if (SC != 1 ):
            xcolumn=12
    elif(XARG=='B32') :
        if(PC == 1) :
            print "Can't plot Q30 when parity is conserved."
            sys.exit()
        xlabel =r'$\beta_{32}$ '
        xfname =PREFIX + '.t.qlm.tab'
        xcolumn=10
        if (SC != 1 ):
            xcolumn=16
    elif(XARG=='OmZ') :
        xlabel =r'$\hbar \omega_{z}$ (MeV) '
        xfname =PREFIX + '.e.tab'
        xcolumn=14
    elif(XARG=='OmX') :
        xlabel =r'$\omega_{x}$ (MeV $\hbar^{-1}$) '
        xfname =PREFIX + '.e.tab'
        xcolumn=10
    elif(XARG=='JX') :
        xlabel =r'$\langle J_{x}\rangle$ ($\hbar$) '
        xfname =PREFIX + '.e.tab'
        xcolumn=11
    elif(XARG=='JZ') :
        xlabel =r'$\langle \hat{J}_{z} \rangle$ ($\hbar$)'
        xfname =PREFIX + '.e.tab'
        xcolumn=15
    elif(XARG=='THETAZ'):
        xlabel =r'$\theta_z$'
        xfname =PREFIX + '.e.tab'
        xcolumn= 14
        altx   = 10
    elif(XARG=='NX'):   
        xlabel =r'$nx$'
        xfname =PREFIX + '.calc.tab'
        xcolumn=3
    elif(XARG=='DX'):   
        xlabel =r'$dx$ (fm)'
        xfname =PREFIX + '.calc.tab'
        xcolumn=6
    elif(XARG=='ALIY'):
        xlabel =r'$\alpha_x (^{\circ})$'
        xfname =PREFIX + '.block.tab'
        xcolumn=1
    else :
        print 'XARG not recognized'
        return

    if(YARG=='E') :
        ylabel =r'E (MeV)'
        yfname =PREFIX + '.e.tab'
        ycolumn=1
        derivY = 0
    elif(YARG=='Routh') :
        ylabel =r'Routhian (MeV)'
        yfname =PREFIX + '.e.tab'
        ycolumn=4
        derivY = 0
    elif(YARG=='E-Routh') :
        ylabel =r'Routh - E (MeV)'
        yfname =PREFIX + '.e.tab'
        ycolumn=1
        derivY = 0
        alty   = 4
    elif(YARG=='MAGIC') :
        ylabel =r'$|E - E_{spwf}|$ (MeV)'
        yfname =PREFIX + '.e.tab'
        ycolumn=1
        alty   =3
        derivY =0
        absy   =1
    elif(YARG=='EFD') :
        ylabel =r'E (MeV)'
        yfname =PREFIX + '.e.tab'
        ycolumn=2
        derivY = 0
    elif(YARG=='B40') :
        ylabel =r'$\beta_{40}$ '
        yfname =PREFIX + '.t.qlm.tab'
        if(PC == 1) :
            ycolumn=6
        else:
            ycolumn=12
        derivY = 0
    elif(YARG=='B42') :
        ylabel =r'$\beta_{42}$ '
        yfname =PREFIX + '.t.qlm.tab'
        if(PC == 1) :
            ycolumn=8
        else:
            ycolumn=14
        derivY = 0
    elif(YARG=='JZ') :
        ylabel =r'$\langle \hat{J}_{z} \rangle$ ($\hbar$)'
        yfname =PREFIX + '.e.tab'
        ycolumn=15
        derivY = 0
    elif(YARG=='JX') :
        ylabel =r'$\langle \hat{J}_{x} \rangle$ ($\hbar$)'
        yfname =PREFIX + '.e.tab'
        ycolumn= 11
        derivY = 0
    elif(YARG=='OmX') :
        ylabel =r'$\omega_{x}$ (MeV $\hbar^{-1}$) '
        yfname =PREFIX + '.e.tab'
        ycolumn=10
        derivY = 0
    elif(YARG=='I2Z') :
        # Dynamical moment of inertia. 
        ylabel =r'$\mathcal{I}^{(2)}$ ($\hbar^2$ MeV$^{-1}$)'
        yfname =PREFIX + '.e.tab'
        ycolumn=15
        derivY = 1
    elif(YARG=='I2X') :
        # Dynamical moment of inertia.
        ylabel =r'$\mathcal{I}^{(2)}$ ($\hbar^2$ MeV$^{-1}$)'
        yfname =PREFIX + '.e.tab'
        ycolumn= 11
        derivY = 1
    elif(YARG=='I1X') :
        # Dynamical moment of inertia.
        ylabel =r'$\mathcal{I}^{(1)}$ ($\hbar^2$ MeV$^{-1}$)'
        yfname =PREFIX + '.e.tab'
        ycolumn= 11
        divy = 1
        derivY = 0

    elif(YARG=='B20') :
        ylabel =r'$\beta_{20}$ '
        yfname =PREFIX + '.t.qlm.tab'
        ycolumn=2
        if(PC != 1) :
            ycolumn = 4
        if(SC != 1) :
            ycolumn = 4
        if(SC !=1 and PC != 1):
            ycolumn = 6
        derivY=0
    elif(YARG=='B2') :
        ylabel =r'$\beta_{2}$ '
        yfname =PREFIX + '.t.ql.tab'
        ycolumn=2
        if( PC!= 1):
            ycolumn = 4
        derivY=0
        fac = np.sqrt(5/( 16 * np.pi))

    elif(YARG=='B22') :
        yfname =PREFIX + '.t.qlm.tab'
        ycolumn=4

        if(PC != 1) :
            ycolumn = 6
        if(SC != 1) :
            ycolumn = 6
        if( SC !=1 and PC != 1):
            ycolumn = 10
        derivY = 0
        ylabel=r'$\beta_{22}$'
    elif(YARG=='B30') :
        ylabel =r'$\beta_{30}$ '
        yfname =PREFIX + '.t.qlm.tab'
        ycolumn=8
        if (SC != 1 ):
            ycolumn=12
        derivY=0
        
        if(PC == 1) :
            ycolumn=-1
    elif(YARG=='ABSB30'):
        ylabel =r'$\beta_{30}$ '
        yfname =PREFIX + '.t.qlm.tab'
        ycolumn=8
        if (SC != 1 ):
            ycolumn=12
        derivY=0

        if(PC == 1) :
            ycolumn=-1 
        absy=1 
    elif(YARG=='B32') :
        if(PC == 1) :
            print "Can't plot Q30 when parity is conserved."
            sys.exit()
        ylabel =r'$\beta_{32}$ '
        yfname =PREFIX + '.t.qlm.tab'
        ycolumn=10
        if (SC != 1 ):
            ycolumn=16
        derivY=0
    elif(YARG=='RRMS') :
        ylabel =r'$\langle r^2 \rangle$ (fm$^2$)'
        yfname =PREFIX + '.e.tab'
        ycolumn=8
        derivY=0
    elif(YARG=='DNP') :
        ylabel =r'$\langle \Delta N \rangle$ '
        yfname =PREFIX + '.ef.tab'
        ycolumn=8
        derivY=0
    elif(YARG=='DNN') :
        ylabel =r'$\langle \Delta N \rangle$'
        yfname =PREFIX + '.ef.tab'
        ycolumn=7
        derivY=0
    elif(YARG=='ESO') :
        ylabel =r'$\langle E_{\rm spin orbit} \rangle$'
        yfname =PREFIX + '.edecomp.tab'
        ycolumn=12
        derivY=0   
    elif(YARG=='ETENSOR') :
        ylabel =r'$\langle E_{\rm tensor} \rangle$'
        yfname =PREFIX + '.edecomp.tab'
        ycolumn=15
        derivY=0    
    
    else :
        print 'YARG not recognized'
        return


    if(SORT=='JX') :
        sortfname =PREFIX + '.e.tab'
        sortcolumn=11
    elif(SORT == ''):
        sortfname  = xfname
        sortcolumn = xcolumn        

    return(xlabel, ylabel, xfname,yfname,sortfname,xcolumn,ycolumn,sortcolumn,derivY,altx,alty, divy, absy, fac)

def MOCCaPlot(XARG, YARG, PREFIX,  PC=1,  SC=1, PC2=1, SC2=1, 
              AXIS     =None,  INTERPOLATE=-1, INTERMIN=None,    
              INTERMAX =None,  LABEL=None, NORMY=1, SPHERNORM=0, 
              LINESTYLE='-' ,  MARKER='', COLOR='', OFFSET=None, 
              XMIN     =None,  XMAX=None, MINRANGE=None, MAXRANGE=None,
              LINEWIDTH=1.0, INTERKIND='cubic', SORT='', INVERTINTERPOL=0,     
              YMAX=None):
    #===========================================================================
    # Function that plots two values obtained in a set of data files, labelled
    # by PREFIX, onto the axes passed into the routine.
    # Currently recognized options for XARG and YARG
    # 'E'   the energy from functional
    # 'Q20' <Q20_t> in fm^2
    #===========================================================================

    #Default to current axis
    if AXIS is None:
        AXIS = plt.gca()

    if(isinstance(PREFIX, list)):
        for i in range(len(PREFIX)):
            (xlabel, ylabel, xfname,yfname,sortname,xcolumn, ycolumn, sortcolumn,derivY,altx,alty, divy, absy, fac) = Determinedata(PREFIX[i], XARG, YARG,PC[i],SC[i], SORT)
    
            dataX=np.loadtxt(xfname,skiprows=1)
            dataY=np.loadtxt(yfname,skiprows=1)
            dataS= np.loadtxt(sortname, skiprows=1)
            
            
            tempx    = dataX[:,xcolumn]
            tempsort = dataS[:,sortcolumn]
            if( ycolumn != -1): 
                tempy=dataY[:,ycolumn]
            else:
                tempy=np.zeros_like(tempx)
            
            if(XMIN != None):
                indices = []
                for j in range(len(tempx)):
                    if(tempx[j] < XMIN[i]):
                        indices.append(j)
                tempx = np.delete(tempx,indices)
                tempy = np.delete(tempy,indices)
                tempsort = np.delete(tempsort,indices)     
            if(XMAX != None):
                indices = []
                for j in range(len(tempx)):
                    if(tempx[j] > XMAX[i]):
                        indices.append(j)
                tempx = np.delete(tempx,indices)
                tempy = np.delete(tempy,indices)
                tempsort = np.delete(tempsort,indices)     
            if(YMAX != None):
                indices = []
                for j in range(len(tempy)):
                    if(tempy[j] > YMAX[i]):
                        indices.append(j)
                tempx = np.delete(tempx,indices)
                tempy = np.delete(tempy,indices)
                tempsort = np.delete(tempsort,indices)            
            
            if(i==0) :
                xdata=tempx
                ydata=tempy
                sortdata = tempsort
            else:
                xdata    = np.append(xdata, tempx)
                ydata    = np.append(ydata, tempy)
                sortdata = np.append(sortdata, tempsort)
            
            ydata= fac * ydata
            if(absy == 1):
                ydata = abs(ydata)
            
    else:
        (xlabel, ylabel, xfname,yfname,sortname,xcolumn, ycolumn, sortcolumn, derivY,altx,alty, divy, absy,fac) = Determinedata(PREFIX, XARG, YARG,PC,SC, SORT)
    
        dataX=np.loadtxt(xfname,skiprows=1)
        dataY=np.loadtxt(yfname,skiprows=1)
        sortdata = np.loadtxt(sortname, skiprows=1)
        xdata=dataX[:,xcolumn]
        sortdata = sortdata[:,sortcolumn]
        if( ycolumn != -1): 
                tempy=dataY[:,ycolumn]
        else:
                tempy=np.zeros_like(tempx)

        ydata = tempy
        if(XMIN != None):
            indices = []
            for j in range(len(xdata)):
                if(xdata[j] < XMIN):
                    indices.append(j)
            xdata = np.delete(xdata,indices)
            ydata = np.delete(ydata,indices)
        if(XMAX != None):
            indices = []
            for j in range(len(xdata)):
                if(xdata[j] > XMAX):
                    indices.append(j)
            xdata = np.delete(xdata,indices)
            ydata = np.delete(ydata,indices)
    
        if(YMAX != None):
            indices = []
            for j in range(len(tempy)):
                if(tempy[j] > YMAX):
                    indices.append(j)
            tempx = np.delete(tempx,indices)
            tempy = np.delete(tempy,indices)
    
        if(altx != 0):
            altxdata = dataX[:, altx]
            xdata = np.arctan2(  altxdata, xdata)*180/np.pi
        if(alty != 0):
            altydata = dataY[:, alty]
            ydata    =  altydata - ydata
        
        if(absy == 1):
            ydata = abs(ydata)
        
        ydata = fac * ydata
        
        if(divy !=0):
            ydata = ydata/xdata
        
    # Sort along sort-data for surety
    sortdata, xdata, ydata = zip(*sorted(zip(sortdata,xdata, ydata)))
    
    xdata = list(xdata)
    ydata = list(ydata)
    
    if(derivY == 1) :
        #=====================================================
        #Flag that tells us to take the finite difference of Y
        # The dynamical moment of inertia is a good example
        deriv = np.zeros((len(ydata)-1))
        for i in range(0,len(ydata)-1):
            deriv[i] = ( ydata[i+1] - ydata[i])/(xdata[i+1] - xdata[i])
        ydata = abs(deriv)
        xdata = xdata[0:-1]
        
    if(INTERPOLATE > 0):
        if(INVERTINTERPOL == 0):
            f     = interp1d(sortdata, ydata,kind=INTERKIND)
            g     = interp1d(sortdata, xdata,kind=INTERKIND)
            if(INTERMIN == None and INTERMAX == None):  
                interx= np.arange(min(sortdata), max(sortdata), INTERPOLATE)
            else:
                interx= np.arange(INTERMIN, INTERMAX, INTERPOLATE)

            ydata = f(interx)
            xdata = g(interx)
        else:
            f     = interp1d(ydata, sortdata, kind=INTERKIND)
            g     = interp1d(xdata, sortdata, kind=INTERKIND)
            if(INTERMIN == None and INTERMAX == None):  
                interx= np.arange(min(ydata), max(ydata), INTERPOLATE)
            else:
                interx= np.arange(INTERMIN, INTERMAX, INTERPOLATE)

            xdata = f(interx)
            ydata = g(interx)

    
    if(MINRANGE == None):
        ymin = min(ydata)
    else:
        ymin = + 100000000
        for i in range(len(ydata)):
            if( MINRANGE[0] < xdata[i] < MINRANGE[1]):
                if( ydata[i] < ymin):
                    ymin = ydata[i]
    
    if(NORMY == 1) :
        ydata = ydata - ymin
    
    if(SPHERNORM==1) :
        spher = 1000000.0
        loc   = 10
        for i in range(len(xdata)):
            if(abs(xdata[i]) < spher ):
                spher = abs(xdata[i])
                loc   = i
        spher = ydata[loc]
        for i in range(len(ydata)):
            ydata[i] = ydata[i] - spher

    if(OFFSET != None):
        ydata = ydata -OFFSET
            
 

    if(COLOR != ""):
    	AXIS.plot(xdata,ydata, label=LABEL, linestyle=LINESTYLE, marker=MARKER, color=COLOR, linewidth=LINEWIDTH)
    else:
	    AXIS.plot(xdata,ydata, label=LABEL, linestyle=LINESTYLE, marker=MARKER, linewidth=LINEWIDTH)
    AXIS.set_xlabel(xlabel)
    AXIS.set_ylabel(ylabel)
    
    if(MINRANGE == None):
       minsearch = [min(xdata), max(xdata)] 
    else:
       minsearch = MINRANGE
       
    if(MAXRANGE == None):
       maxsearch = [min(xdata), max(xdata)] 
    else:
       maxsearch = MAXRANGE

    if(INTERPOLATE > 0.01):               
        maskmin = (xdata <= minsearch[1]) &  (xdata >= minsearch[0])
        maskmax = (xdata <= maxsearch[1]) &  (xdata >= maxsearch[0])

        emax = np.max(ydata[maskmax])
        emin = np.min(ydata[maskmin])

        minx = xdata[maskmin][np.argmin(ydata[maskmin])]
        maxx = xdata[maskmax][np.argmax(ydata[maskmax])]
    else:
        xdata = np.array(xdata)
        ydata = np.array(ydata) 

        maskmax = (xdata <= maxsearch[1]) &  (xdata >= maxsearch[0])
        maskmin = (xdata <= minsearch[1]) &  (xdata >= minsearch[0])
        
        argmin = np.argmin(ydata[maskmin])
        minx   = xdata[maskmin][argmin]
        emin = min(ydata[maskmin])
        argmax = np.argmax(ydata[maskmax])
        maxx = xdata[maskmax][argmax]
        emax = max(ydata[maskmax])
    

    return (minx, ymin, maxx, emax)

#################################################################################
def mini(PREFIX, XARG, YARG, PC=1, SC=1, XRANGE=[]):
    
    if(XARG=='Q20') :
        xfname =PREFIX + '.t.qlm.tab'
        xcolumn=1
        if(PC != 1) :
            xcolumn = 3
    elif(XARG=='B20') :
        xfname =PREFIX + '.t.qlm.tab'
        xcolumn=2

        if(PC != 1) :
            xcolumn = 4
        if(SC != 1) :
            xcolumn = 4
        if( SC !=1 and PC != 1):
            xcolumn = 6
    elif(XARG=='B30') :
        if(PC == 1) :
            print "Can't plot Q30 when parity is conserved."
            sys.exit()
        xlabel =r'$\beta_{30}$ '
        xfname =PREFIX + '.t.qlm.tab'
        xcolumn=8
        if (SC != 1 ):
            xcolumn=12
            
    if(YARG=='E') :
        ylabel =r'$E$ (MeV)'
        yfname =PREFIX + '.e.tab'
        ycolumn=1
        derivY=0
    
    elif(YARG=='RRMS') :
        ylabel =r'$\langle r^2 \rangle$ (fm$^2$)'
        yfname =PREFIX + '.e.tab'
        ycolumn=8
        derivY=0
        
    elif(YARG=='OmX') :
        ylabel =r'$\omega_{x}$ (MeV $\hbar^{-1}$) '
        yfname =PREFIX + '.e.tab'
        ycolumn=10
        derivY = 0
    elif(YARG=='JX') :
        ylabel =r'$\langle \hat{J}_{x} \rangle$ ($\hbar$)'
        yfname =PREFIX + '.e.tab'
        ycolumn= 11
        derivY = 0
    
    efname = PREFIX + '.e.tab'
    
    dataE=np.loadtxt(efname,skiprows=1)
    dataX=np.loadtxt(xfname,skiprows=1)
    dataY=np.loadtxt(yfname,skiprows=1)

    xdata=dataX[:,xcolumn]
    ydata=dataY[:,ycolumn]
    edata=dataE[:,1]
    
    if(len(XRANGE) != 0):
        interx= np.arange(XRANGE[0], XRANGE[1], (max(xdata) - min(xdata))/1000)

    else:
        interx= np.arange(min(xdata), max(xdata)- (max(xdata) - min(xdata))/100, (max(xdata) - min(xdata))/1000)
    f     = interp1d(xdata, edata,kind='cubic')
    g     = interp1d(xdata, ydata)
    
    E     = f(interx)
    emin  = np.argmin(E) 
    minx  = interx[emin]
    miny  = g(minx)
    
    return (minx, miny)
################################################################################
def Nilsson(PREFIX, BASIS, ISO, PAR, SIG, KMAX=0, AXIS=None, INTERPOLATE=-1, 
            PLOTDATA=-1, MARKER='', FERMIWINDOW=-1, XARG='B20', PC=1, SC=1, HF=0
            , LINESTYLES=['-.','-', '--'], LINEWIDTH=1, DASHES=None, SIMPLEX=0, 
            SIMPLEXSORT=0, MARKERSIZE=1):

    #Default to current axis
    if AXIS is None:
        AXIS = plt.gca()

    
    (xlabel, ylabel, xfname,yfname,sortfname,xcolumn,ycolumn,sortcolumn,derivY,altx,alty, divy, absy, fac)= Determinedata(PREFIX, XARG,'E' ,PC,SC, SORT='')
    
    dataX=np.loadtxt(xfname,skiprows=1)
    xdata = dataX[:,xcolumn]
    
    fermifname=PREFIX + '.ef.tab'
    if(ISO == 'neutron'):
        fermicolumn=1
    else :
        fermicolumn=2
    fermidata = np.loadtxt(fermifname, skiprows=1)
    fermi     = fermidata[:,fermicolumn]

    xdata, fermi = zip(*sorted(zip(xdata, fermi)))
    AXIS.plot(xdata, fermi, 'k'+LINESTYLES[0], label='$\epsilon_{\\rm fermi}$')
    xdata = dataX[:,xcolumn]

    minx = 10000000.0
    minj = 0
    for j in range(len(xdata)):
        if(abs(xdata[j]) < abs(minx)):
            minj = j
            minx = xdata[j]
    fermispher = fermi[minj]
    
    yplus = max(fermi)
    ymin  = min(fermi)
    
    shells = {}
    
    if(FERMIWINDOW > 0):
        AXIS.set_ylim((ymin-FERMIWINDOW, yplus + FERMIWINDOW))

    colors   = ['b', 'r', 'c', 'g', 'k', 'm', 'y', 'salmon', 'y' ]  
    
    for P in PAR:
        shells[P] = {}
        for S in SIG:
            if(SIMPLEXSORT == 0):
                fnametemp=PREFIX + '.' + BASIS + '.' + ISO + '.' + 'par=' + P + '.' + 'sig=' + S 
                
                if (P == '-1'):
                    linestyle = LINESTYLES[2]
                else:
                    linestyle = LINESTYLES[1]
                
            else:
                fnametemp=PREFIX + '.' + BASIS + '.' + ISO + '.' + 'sx=' + S 
                
                if (S == '-1'):
                    linestyle = LINESTYLES[2]
                else:
                    linestyle = LINESTYLES[1]
            
           
            if(KMAX !=0):
                #===============================================================
                # Axial case
                for K in range(1,KMAX+1,2):
                    
                    shells[P][K]= []
                    fname = fnametemp + '.k=%d.tab'%K
                    c = colors[(K-1)/2]
                    try:
                        tokens = tokenizer(fname)
                        tokens.next()
                        spwfs = [np.loadtxt(A) for A in tokens]
                    except IOError:
                        #This happens if KMAX is too big
                        break
                    for i in range(len(spwfs)):
                        spwf = spwfs[i]

                        try:
                            ydata = spwf[:,7]
                            indexes = [int(x)-1 for x in spwf[:,1]]
                        except IndexError:
                            # Just ignore spwfs that are a single point
                            continue
                        
                        xdata   = dataX[indexes,xcolumn]
                        if(PLOTDATA == 1):
                            AXIS.plot(xdata, ydata, c +'x')
                        if(INTERPOLATE > 0):
                            try:
                                f     = interp1d(xdata, ydata, kind='cubic')
                            except:
                                continue
                            interx= np.arange(min(xdata), max(xdata) - (max(xdata) - min(xdata))/100, INTERPOLATE) 
                            ydata = f(interx)
                            xdata = interx
                            
                     
                        if((spwf[0,3] == -1)):
                            m  = (5,1)
                            me = 20
                            ms = 10
                        else:
                            m  = MARKER
                            me = 1
                            ms = MARKERSIZE
                        lw = LINEWIDTH
                        
                        xdata, ydata = zip(*sorted(zip(xdata, ydata)))
                        if(i == 0 and P == PAR[0] ):
                            if(DASHES != None):
                                AXIS.plot(xdata, ydata, color=c, linestyle=linestyle, label=r'$J_z = \frac{%d}{2}$'%K, marker=m,linewidth=lw, dashes=DASHES, markevery=me, ms=ms)
                            else :
                                AXIS.plot(xdata, ydata, color=c, linestyle=linestyle, label=r'$J_z = \frac{%d}{2}$'%K, marker=m,linewidth=lw, markevery=me, ms=ms)
                        
                        else:
                            if(DASHES != None):
                                AXIS.plot(xdata, ydata, color=c, linestyle=linestyle, marker=m, linewidth=lw, dashes=DASHES, markevery=me, ms=ms)
                            else:
                                AXIS.plot(xdata, ydata, color=c, linestyle=linestyle, marker=m, linewidth=lw, markevery=me,ms=ms)
                                
                        # Find the point of the spwf closest to x = 0 and add it
                        # to shells
                        minx = 10000000.0
                        minj = 0
                        for j in range(len(xdata)):
                            if(abs(xdata[j]) < abs(minx)):
                                minj = j
                                minx = xdata[j]
                        if(abs(minx) < 0.05 ):
                            shells[P][K].append(ydata[minj])
                            
            else:
               #################################################################
               # NONAXIAL CASE
                fname = fnametemp + '.tab'
                c='k'
                try:
                    tokens = tokenizer(fname)
                    tokens.next()
                    try:
                        spwfs = [np.loadtxt(A) for A in tokens]
                    except ValueError:
                        print tokens.next()
                        exit()
                except IOError:
                    break
                for i in range(len(spwfs)):
                    spwf = spwfs[i]
                    try:
                        ydata = spwf[:,7]
                        
                        if(SIMPLEX == 1):
                            if (spwf[3,12] < 0):
                                linestyle = LINESTYLES[2]
                            else:
                                linestyle = LINESTYLES[1]
                        indexes = [int(x)-1 for x in spwf[:,1]]
                        if(HF == 1):
                            occupation = spwf[4,6]
                            if(occupation == 1):
                                c = 'r'
                            else:
                                c = 'k'
                    except IndexError:
                        # Just ignore spwfs that are a single point
                        continue
                        
                    if((spwf[0,3] == -1)):
                        m  = (5,1)
                        me = 72
                        ms = 10
                    else:
                        m  = MARKER
                        me = 1
                        ms = MARKERSIZE
                        lw = LINEWIDTH
                            
                        
                    xdata   = dataX[indexes,xcolumn]
                    if(PLOTDATA == 1):
                        AXIS.plot(xdata, ydata, c +'x')
                    if(INTERPOLATE > 0):
                        try:
                            f     = interp1d(xdata, ydata, kind='cubic')
                        except:
                            continue
                        interx= np.arange(min(xdata), max(xdata) - (max(xdata) - min(xdata))/100, INTERPOLATE) 
                        ydata = f(interx)
                        xdata = interx
                    if(i == 0 and P == PAR[0] ):
                        AXIS.plot(xdata, ydata, c+linestyle, marker=m, ms=ms, markevery=me)
                    else:
                        AXIS.plot(xdata, ydata, c+linestyle, marker=m, ms=ms, markevery=me)
               
    AXIS.set_xlabel(xlabel)
    AXIS.set_ylabel(r'E (MeV)')
    
    # Return the shell structure at the spherical point (or at least as close)
    # as possible.
    return (fermispher,shells)

################################################################################
def Qps(PREFIX, PAR, ISO, AXIS=None, INTERPOLATE=-1, PLOTDATA=-1, MARKER=''):
    
    #Default to current axis
    if AXIS is None:
        AXIS = plt.gca()

    xfname=PREFIX + '.e.tab'
    dataX=np.loadtxt(xfname,skiprows=1)
    xdata = dataX[:,11]
    
    
    for P in PAR:
        if (P == '-1'):
            linestyle = '--'
        else:
            linestyle = '-'

        fname =PREFIX + '.%s.qp.'%ISO + 'P=' + P + '.tab'
        efname=PREFIX + '.ef.tab'
       
        l2data = np.loadtxt(efname, skiprows=1)
        l2     = l2data[:,3]

#        try:
        tokens = tokenizer(fname)
        tokens.next()
        spwfs = [np.loadtxt(A) for A in tokens]
#        except IndexError:
#            continue

        for i in range(len(spwfs)):
            spwf = spwfs[i]
            try:
                ydata = spwf[:,3]
                indexes = [int(x)-1 for x in spwf[:,1]]
            except IndexError:
                # Just ignore spwfs that are a single point
                continue
            xdata   = dataX[indexes,11]
            
            xdata, ydata = zip(*sorted(zip(xdata,ydata)))
            xdata = np.asarray(xdata)
            ydata = np.asarray(ydata)

            if(PLOTDATA == 1):
                AXIS.plot(xdata, ydata, c +'x')
            if(INTERPOLATE > 0):
                try:
                    f     = interp1d(xdata, ydata, kind='cubic')
                except:
                    continue
                interx= np.arange(min(xdata), max(xdata) - (max(xdata) - min(xdata))/100, INTERPOLATE) 
                ydata = f(interx)
                xdata = interx
            if(i == 0 and P == PAR[0] ):
                AXIS.plot(xdata,  ydata, linestyle=linestyle, marker=MARKER)
                AXIS.plot(xdata, -ydata, linestyle=linestyle, marker=MARKER)
            else:
                AXIS.plot(xdata,  ydata, linestyle=linestyle, marker=MARKER)
                AXIS.plot(xdata, -ydata, linestyle=linestyle, marker=MARKER)

    AXIS.set_ylabel(r'$E_{qp}$ (MeV)')

