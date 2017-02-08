import numpy as np
from tokenizer import tokenizer

def Sort(fname, KMAX=15):
    #===========================================================================
    # STEP 1
    #   Read a table of spwfs generated by MOCCa.spectra.awk
    #   into an array of all of the spwfs
    tokens = tokenizer(fname)
    A = tokens.next() # Skip the first line
    
    with open(fname, 'r') as f:
        header = f.next()

    Ecolumn     = 7
    Jcolumn     = 10
    humancolumn = 2
    prefix  = fname.replace('.tab', '') 
    
    # This is a list of length NWT of a set of wavefunctions at every 
    # deformation
    spwfs  =[np.loadtxt(A) for A in tokens]
    
    #Define some dimensions
    nwt          = len(spwfs)
    spwf         = spwfs[0]
    deformations = len(spwf[:,0])
    info         = len(spwf[0,:])
    
    #Make new array out of all the data
    # Dimension 1: deformation
    # Dimension 2: values of the SPWF
    # Dimension 3: number of SPWF in the deformation block
    spwfs        = np.dstack(spwfs)
    
    #===========================================================================
    # STEP 2
    #   Tentatively sort all the points in the graph by their <Jz> value. 
    #   This uses Michaels definition of 'confidence intervals' for their
    #   <Jz> value. 
    K = np.zeros((deformations,nwt))
    for j in range(deformations):
        for i in range(nwt):
            kkk = spwfs[j,Jcolumn,i]
            hk  = spwfs[j,humancolumn,i]
            
            if(hk == 0) : 
                if (  -0.3 < kkk and kkk <  1.2 ):
                    twok =  1.0
                if (   1.2 < kkk and kkk <  3.2 ):
                    twok =  5.0
                if (   3.4 < kkk and kkk <  5.2 ):
                    twok =  9.0
                if (   5.4 < kkk and kkk <  7.2 ):
                    twok = 13.0
                if (   7.4 < kkk and kkk <  9.2 ):
                    twok = 17.0
                if (   9.4 < kkk and kkk < 11.2 ):
                    twok = 21.0
                if (  -2.5 < kkk and kkk < -0.4 ):
                    twok =  3.0
                if (  -4.2 < kkk and kkk < -2.6 ):
                    twok =  7.0
                if (  -6.2 < kkk and kkk < -4.4 ):
                    twok = 11.0
                if (  -8.2 < kkk and kkk < -6.7 ):
                    twok = 15.0
                if ( -10.2 < kkk and kkk < -8.7 ):
                    twok = 19.0
            
                #Assign to the mesh
                K[j,i] = twok
            else:
                K[j,i] = hk
    #===========================================================================
    # STEP 3
    #    Construct chains throughout the data
    #
#    Taken = []
#    Chains={}
#    for k in range(1,KMAX,2):
#        Chains[k] = []    
#        for i in range(nwt):
#            Chains[k].append([])
#            for j in range(deformations):
#                for ii in range(nwt): 
#                    if(k == K[j,ii] and (ii,j) not in Taken):
#                        Chains[k][i].append(ii)
#                        Taken.append((ii,j))
#                        break
    #===========================================================================
    # STEP 4
    #   Detect jumps and remove
#    for k in range(1,KMAX,2):
#        for i in range(len(Chains[k])):
#            chain = Chains[k][i]
#            for j in range(1,len(chain)-1):
#                ind  = chain[j]
#                indp = chain[j+1]
#                indm = chain[j-1]
#                Ehere = spwfs[j,Ecolumn,ind]
#                Eplus = spwfs[j+1,Ecolumn,indp]
#                Emin  = spwfs[j-1,Ecolumn,indm]
#                if( abs(Ehere - 0.5*(Eplus + Emin)) > 5 ):
#                    chain[j] = -1
#    
  
    #===========================================================================
    # STEP 5
    #    Output into gle readable tables
    for k in range(1,KMAX,2):
        with open('%s.k=%d.tab'%(prefix,k), 'w') as f:
            # Write the header
            f.write(header)
            for i in range(len(Chains[k])):
                chain = Chains[k][i]
                for j in range(len(chain)):
                   ind = chain[j]
                   if (ind != -1) :
                      for l in range(info):
                        f.write(' %8.3f '%spwfs[j,l,ind])
                   f.write('\n')
                f.write('*\n')