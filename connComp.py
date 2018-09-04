#! C:\Python27\python

import numpy as np
import scipy.sparse as ss
import sys

def isEmpty(A):
    return A.shape[0]==0

def connComp(C,flagCCsingle):

# syntax: CC = conn_comp (C)
# returns the matrix CC. Each line of CC is a connected component with the
# list of the connected nodes of that component. Returns CC=[] if there is
# no connection in the matrix (no CC of length bigger than 1)
# matrix C is square

    ACT=np.array(np.arange(1,C.shape[0]))
    #if len(sys.argv)==1:
    #  flagCCsingle=0  # if ==1 returns CC of length 1 as well
    I,J=np.nonzero(C!=0)
    connected=np.union1d(I,J)   #list of connected vertices
    if isEmpty(connected): #if there are no connections, return empty matrix
       CC=np.array(np.arange(1,C.shape[0]))
       return CC
    NCC=0; #initialise number of CC

    # As soon as I put a connected vertice in a connected component I delete it
    # from the vector "connected". This is why there is this loop "while"

    while not isEmpty(connected):
      NCC=NCC+1
      i=connected[0]
      CCaux=np.array([i])  #initialise this new Connected Compound
      k=0
      neighbour=np.union1d(J[np.nonzero(I==i)],I[np.nonzero(J==i)]) # compute the list of neighbours
      CCaux=np.concatenate((CCaux, np.setdiff1d(neighbour,CCaux)),0) # I add the NEW vertices that I just found to the CC
      kL=CCaux.shape[0]
      while k < kL-1:
        j=CCaux[k+1]  #for each of those added neighbours
        neighbour=np.union1d(J[np.nonzero(I==j)],I[np.nonzero(J==j)]) #I compute there neighbours
        CCaux=np.concatenate((CCaux, np.setdiff1d(neighbour,CCaux)),0) #and add the new ones
        k=k+1
      connected=np.setdiff1d(connected,CCaux); # I don't have to look the vertices in connected AND CCaux
      #ccL=CCaux.shape[0]
      #CC=np.zeros((NCC,ccL))
      #CC[NCC,1:ccL]=CCaux #I save the connected compound I just computedconnected=union(I,J);
      CC=CCaux
    connected=np.union1d(I,J)
    if flagCCsingle==1:
      notConnected=np.setdiff1d(ACT,connected)
      #ncL=notConnected.shape[0]
      #CC[NCC+1:NCC+ncL,1]=notConnected
      CC=notConnected
    return CC

if __name__ == '__main__':
    flagC=0
    C=np.array([[1.,2.,3.], [4.,5.,6.],[7.,8.,9.]])
    CC=connComp(C,flagC)
    print(CC)

     
