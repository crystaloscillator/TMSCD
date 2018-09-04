#def DMSCD_partition_stab2(G,sA,jstab,lapFlag,filterFlag,Nscales,Nsignals,njobs):
#Nsccales - number of scales
#G-cell where each entry is network at one time point
#sA - supra Adjacency matrix
#Nsignals - number of signals used to approximate wavelets
#jstab-number of times process should be repeated to obtain instability - 6 is recommended for very large networks
#lapFlag - the type of laplacian network to be used
#filterFlag - filter kernel function to be used 'wavelet' or 'wavelet2'
#njobs - number of cores you want to run the pairwise distance computation on

import numpy as np
import scipy.sparse as ss
import json
from ftplib import FTP
import os


def DMSCDpartitionStab(sNum,G,sA,jstab,lapFlag,filterFlag,Nscales,Nsignals,njobs):
    T=np.size(G)
    n=G[0].shape[0]
    N=sA.shape[0]
    if T==1:
        if lapFlag=='lap':
            QA=MSCDCreateLaplaciangraph(G)
        elif lapFlag=='normlap':
            QA=MSCDCreateNormLaplaciangraph(G)
        elif lapFlag=='rwlap':
            QA=MSCDCreateRwLaplaciangraph(G)
        else :
            print('MSCD_partition_RV : lap_flag should be equal to ''lap'', ''normlap'', or ''rwlap''')
    # estimate largest eigenvalue and the three smallest : 
    #opts.isreal=1
    #opts.issym=1
    #opts.tol=1e-3
        V,D = np.linalg.eig(QA,1)  #[V,D] = eigs(QA,1,'SA'); [V,D] = eigs(QA,1,'SA'); v1=[]
        v1[0,:]=V
    else:
        for t in range(0,T):
            if lapFlag=='lap':
                QA=MSCDCreateLaplaciangraph(G[t])
            elif lapFlag=='normlap':
                QA=MSCDCreateNormLaplaciangraph(G[t])
            elif lapFlag=='rwlap':
                QA=MSCDCreateRwLaplaciangraph(G[t])
            else :
                print('MSCD_partition_RV : lap_flag should be equal to ''lap'', ''normlap'', or ''rwlap''')
            V,D = np.linalg.eig(QA,1)  #[V,D] = eigs(QA,1,'SA'); [V,D] = eigs(QA,1,'SA'); v1=[]
            v1[t,:]=V
    #Setting up parameters and obtaining wavelets
    if sNum == 1 or sNum == 2 :
        tg,cg,g0,arange,L=DMSCDpartitionRVnoStab(2,Nscales,sA,T,v1,Nsignals,lap_flag,filter_flag)
    else:
       tg,cg,g0,arange,L=DMSCDpartitionRVnoStab(3,Nscales,sA,T,v1,Nsignals,lap_flag,filter_flag) 
    for j in range(0,jstab):
       print('J iteration:'+j)
       if sNum == 1 or sNum == 3:
           DMSCDgetWavelets(N,Nsignals,L,cg,arange,Nscales);
           #Now we will use python to produce fast agglomerative 
           #clustering due to the dymensionality of the data
           p=' '
           p1=str(T)
           p2=str(n)
           p3=str(Nscales)
           p4='correlation'
           p5=str(njobs)
           p01='python'
           p02=str(FTP.pwd()) +'Clustering_Benchmarks.py'
           commandStr=p01+p+p02+p+p1+p+p2+p+p3+p+p4+p+p5
           status=execfile(commandStr)
           for t in range(0,Nscales):
               f=open('WsigTree.mat'+t,'r')
               tree=json.load(f)
               f.close()
               tree=tree.tree
               tree=int64(tree)
               tree[:,0]=tree[:,0]+1
               tree[:,1]=tree[:,1]+1
               COMcomplete[j,:,t]=MSCDcutDendrogram2(tree)
               os.system('rm Wsig*.mat')

       elif sNum == 2 :
           sig=np.random(N,Nsignals)
           Wsig=sgwtChebyOp(sig,L,cg,arange); # wavelet coefficients
           for t in range(0,Nscales):
               print('scale '+t)
               v=Wsig[t]
               f=open('Wsig.mat','r+')
               json.dump(v,f)
               f.close()
           #Now we will use python to produce fast agglomerative 
           #clustering due to the dymensionality of the data
               p=' '
               p1=str(T)
               p2=str(n)
               p3=str(Nscales)
               p4='correlation'
               p5=str(njobs)
               p01='python'
               p02=str(FTP.pwd()) +'Clustering_Benchmarks2.py'
               commandStr=p01+p+p02+p+p1+p+p2+p+p3+p+p4+p+p5
               status=execfile(commandStr)
               f=open('WsigTree.mat','r')
               tree=json.load(f)
               f.close()
               tree=tree.tree
               tree=int64(tree)
               tree[:,0]=tree[:,0]+1
               tree[:,1]=tree[:,1]+1
               COMcomplete[j,:,t]=MSCDcutDendrogram2(tree)
               os.system('rm Wsig.mat')
               os.system('rm WsigTree.mat')
       else :
           print('sNum must be 1,2 or 3')
    if jstab>1:
        jj=0
        for j in range(0,jstab): 
            for k in range(j+1,jstab):
                jj=jj+1
                for t in range(0,Nscales):
                    stab[jj,t]=PartAgreeCoef_AR(COMcomplete[j,:,t],COMcomplete[k,:,t])
    if jstab>2:
        stab=np.mean(stab)
    return tg,COMcomplete,stab

        

