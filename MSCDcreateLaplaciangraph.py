#! C:\Python27\python

import numpy as np
import scipy.sparse as ss

#checking whether A is triungular
def isTriu(A):
   B=ss.triu(A)==A  
   return B.all()

#checking whether A is  symmetric
def isSym(A):
   D=A==A.transpose()
   return D.all()

def MSCDCreateLaplaciangraph(A):
   if isTriu(A):  
        A=A+A.transpose()
   elif ~isSym(A):  
          print('MSCD_create_laplaciangraph: A should be tri sup or sym')
   Q=np.diag(np.sum(A,0))-A 
   return Q

def MSCDCreateNormLaplaciangraph(A):
   if isTriu(A):  
        A=A+A.transpose()
   elif ~isSym(A):  
          print('MSCD_create_normlaplaciangraph: A should be tri sup or sym')
   DEGDem=ss.csr_matrix(np.diag(1./np.sqrt(sum(A,0))))
   Q=np.eye(A.shape[0])-DEGDem*A*DEGDem
   # make sure the Laplacian is exactly symmetrical:
   Qbis=ss.triu(Q,1)
   Qbis=Qbis+Qbis.transpose()
   Ln=Qbis+np.diag(np.diag(Q))
   Ln[np.isnan(Ln)] = 0 
   return Ln

def MSCDCreateRwLaplaciangraph(A):
   if isTriu(A):  
        A=A+A.transpose()
   elif ~isSym(A):  
          print('MSCD_create_rwlaplaciangraph: A should be tri sup or sym')
   DI = 1./sum(A,0)
   Lrw=np.eye(A.shape[0])-ss.csr_matrix(np.diag(DI))*A;
   Lrw[np.isnan(Lrw)] = 0 
   return Lr

if __name__ == '__main__':    
    A=np.array([[1.,2.,3.], [2.,5.,6.],[3.,6.,7.]])
    L=MSCDCreateLaplaciangraph(A)
   
