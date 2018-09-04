#! C:\Python27\python
#
# function g=MSCD_wavelet_kernel_function(a,b,t2,x)
#
# Creates band pass wavelet filter kernel
#
# Inputs : 
# - a, which is the parameter Alpha, coding the behavior of the wavelet kernel at low eigenvalue
# - b, which is the parameter Beta, coding the behavior of the wavelet kernel at high eigenvalue
# - t2, the parameter deciding where does the wavelet kernel start decaying as a power law
# - x is the real value at which one wants to evaluate g
#
# Output :
# - g, which is in fact the value of the function at x: g(x)

import numpy as np

def MSCD_wavelet_kernel_function(a,b,t2,x):
    t1=1
    g=np.zeros(x.shape(0))
    a1=(a*t2 - b*t1)/(t1*t2*(t1 - t2)^2)
    a2=-(2*a*t2^2 - 2*b*t1^2 + a*t1*t2 - b*t1*t2)/(t1*t2*(t1 - t2)^2)
    a3=(- b*t1^3 - 2*b*t1^2*t2 + 2*a*t1*t2^2 + a*t2^3)/(t1*t2*(t1 - t2)^2)
    a4=(b*t1^2 - a*t2^2 - 2*t1*t2 + t1^2 + t2^2)/(t1^2 - 2*t1*t2 + t2^2)
    
    g1=np.nonzero(x>=0 & x<t1)[0]
    g2=np.nonzero(x>=t1 & x<t2)[0]
    g3=np.nonzero(x>=t2)[0]
    g(g1)=x(g1)#.^a*t1^(-a)






