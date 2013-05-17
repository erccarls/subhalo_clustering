import numpy as np
from scipy import weave

x = np.array(((1,2),(3,4)))



#support_code = "#include <stdio.h>"
#code = """
#    
#    std::cout << x(0,1) << x(1,1)  << std::endl;
#    
#"""
#weave.inline(code,['x'],support_code = support_code, compiler='gcc', type_converters= weave.converters.blitz) 

#print weave.inline(code,['radius','BGTemplate','size','x','y','start'],support_code = support_code, compiler='gcc')


a = np.random.poisson([5,10,50])
print a
print np.append(np.array([0,]),np.cumsum(a))