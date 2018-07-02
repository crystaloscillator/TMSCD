#! C:\Python27\python

# sgwtChebyOp : Chebyshev polynomial of Laplacian applied to vector
#
# def sgwtChebyOp(f,L,coef,arng):
#
# Compute polynomials of laplacian (in Chebyshev
# basis) applied to input.
#
# Coefficients for multiple polynomials may be passed as a array of array. 
#
# f- input vector
# L - graph laplacian 
# coef - Chebyshev coefficients. If c is a plain array, then they are
#     coefficients for a single polynomial. If c is a array of array, 
#     then it contains coefficients for multiple polynomials, such 
#     that c[j][1+k] is k'th Chebyshev coefficient the j'th polynomial.
# arng - interval of approximation
#
# Outputs:
# rr - result.  r is carray of vectors size of f. 
#    


import numpy as np



def sgwtChebyOp(f,L,coef,arng):
  # to be decided later
  # if ~iscell(c)
  #   r=sgwt_cheby_op(f,L,{c},arange);
  #  r=r{1};
  #  return;
  #end
  Nscales=coef.shape[0]
  M=np.zeros(Nscales)
  for j in range(0,Nscales):
      M[j]=np.size(coef[j])
  maxM=int(np.max(M))
  a1=(arng[1]-arng[0])/2
  a2=(arng[1]+arng[0])/2
  twfOld=f
  twfCur=(L.dot(f)-a2*f)/a1
  rr=np.zeros((Nscales,Nscales))
  for j in range(0,Nscales):
    rr[j]=np.array(.5*coef[j][0]*twfOld + coef[j][1]*twfCur)
  for k in range(1,maxM-1):
    twfNew = (2./a1)*(L.dot(twfCur)-a2*twfCur)-twfOld
    for j in range(0,Nscales):
      if 1+k<=M[j]:
          rr[j]=rr[j]+coef[j,k+1]*twfNew
    twfOld=twfCur
    twfCur=twfNew;
  return rr
          
print('Hello')
L=np.array([[1.,2.,3.], [2.,5.,6.],[3.,6.,7.]])
f=np.array([1.,2.,3.])
c=np.array([[1.,2.,3.],[4.,2.,5],[3.,1.,4.]]) 
arng=np.array([-1.,1.])
r=sgwtChebyOp(f,L,c,arng)




