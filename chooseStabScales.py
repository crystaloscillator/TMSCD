#! C:\Python27\python

import numpy as np
import scipy.sparse as ss
import partAgreeCoef_ARonly as pa

def chooseStabScales(stab,N,COM_RV):

    Nscales=stab.shape[0]
    # find local minima:
    aux=np.diff(1-stab)
    L=-1
    stabScales=np.zeros(Nscales)
    for i in range(1,Nscales-1):
      if (np.sign(aux[i])>=0) and (np.sign(aux[i-1])==-1):
          L=L+1
          stabScales[L]=i
    if 1-stab[0]<=1-stab[1]:
      L=L+1
      stabScales[L]=1
      
    # do not consider partitions with more than N/2 communities:
    notAnalyzable=np.nonzero(np.max(COM_RV,0)>N/2)
    #notAnalyzable=notAnalyzable[1]
    stabScales=np.setdiff1d(stabScales,notAnalyzable)
    if stabScales.shape[0]==0: #isempty(stabScales)
      if 1-stab[-1]<=1-stab[-2]:
        stabScales=np.zeros(1)
        stabScales[0]=Nscales
        L=L+1
    if stabScales.shape[0]==0: #isempty(stabScales)
      print('No local minima of 1-stab has less than N/2 communities')
      return stabScales
    # look at similarity between detected poartitions:
    M=stabScales.shape[0]
    simi=np.zeros((M,M))
    for i in range(0,M):
      for j in range(i+1,M):
        simi[i,j]=pa.partAgreeCoef_AR(COM_RV[:,stabScales[i]],COM_RV[:,stabScales[j]]);
    I,J=np.nonzero(simi==1)
    K=np.unique([I,J])
    simi=np.array([simi,np.zeros((0,M))])
    simi=simi+simi.transpose()
    # if two (or more) partitions are exactly the same, 
    # keep the one with the largest stability:
    lookedAt=[]
    stabScalesRemove=[]
    KL=K.shape[0]
    for l in range(0,KL):
      k=K[l]
      if ~isMember(k,lookedAt):
        i=np.nonzero(simi[:,k]==1)
        ind=np.min(1-stab[stabScales[k,i]]) #[~,ind]  i'
        if ind==1:
            stabScalesRemove=[stabScalesRemove,stabScales[i]]
        else:
            stabScalesRemove=[stabScalesRemove,stabScales[[k,np.setdiff1d(i,i[ind-1])]]];
        lookedAt=np.unique([lookedAt,i,k])

    stabScales=np.setdiff1d(stabScales,stabScalesRemove)
    return stabScales

#return >=0 if k is memeber of sVec
#else -1
def isMember(k,sVec):
    l=sVec.shape(0)
    for i in range(0,l):
        if k==sVec[i]:
            return i
    return -1

if __name__ == '__main__':
    stab= np.array([1.,2.,3.])
    N=3 
    comRV=np.array([[1.,2.,3.], [4.,5.,6.],[7.,8.,9.]])
    resS=chooseStabScales(stab,N,comRV)
    print(resS)



