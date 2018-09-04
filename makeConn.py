#! C:\Python27\python

import numpy as np

def makeConn(A,CC):
    n=CC.shape[0]
    max=int(np.max(CC))
    Anew=np.zeros((max,max))
    for i in range(0,n):
        for j in range(i+1,n):
          aa=CC[i,:]
          aaa=aa[aa>0]
          if aaa.shape[0]==1:
              a=aaa
          else:
              a=rndSmpl(aaa) #np.random.randint(1,5) #random_sample(3) #rand(aaa) #,1)
          bb=CC[j,:]
          bbb=bb[bb>0]
          if bbb.shape[0]==1:
              b=bbb
          else:
              b=rndSmpl(bbb) #,1)
          Anew[a,b]=1
          Anew[b,a]=1
    sing=[]
    for i in range(0,n):
        a=CC[i,:]
        b=a[a>0]
        if b.shape[0]==1:
            sing=np.array(sing, b)
    return Anew,sing       

def rndSmpl(vec):
    n=vec.shape[0]
    ind=np.random.randint(0,n-1)
    return vec[ind]

if __name__ == '__main__':
    flagC=0
   # C=np.array([[1.,2.,3.], [4.,5.,6.],[7.,8.,9.]])
    C=np.array([[1,2,3], [4,5,6],[7,8,9]])
    A=np.array([[1.,6.,3.], [4.,5.,7.],[6.,8.,2.]])
    CC,s=makeConn(A,C)
    print(CC)
