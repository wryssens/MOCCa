#------------------------------------------------------------------------------------
# Copy the human assignments from an already processed file of spwfs to another one.
#
#------------------------------------------------------------------------------------

import numpy as np
from tokenizer import *
import os

def copyassignments( fexample, ftarget ):

    column=2

    ex =  open(fexample, 'r') 
    ta =  open(ftarget , 'r')
    out=  open('temp', 'w')
    
    
    for tarline in ta:
        
        exline = ex.next()
        
        tasp = tarline.split()
        exsp =  exline.split()
        
        if(len(exsp) > 1) :
            tasp[column] = exsp[column]
        
        for i in range(len(tasp)):
            out.write(tasp[i] + ' ') 
        out.write('\n')
        
    ex.close()
    out.close()
    ta.close()
    os.system('mv temp ' + ftarget )
