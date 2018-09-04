#! C:\Python27\python

#
# def DMSCDpartitionRVnoStab(sNum,Nscales,A,Nsignals,lap_flag,filter_flag)
#
# Computes the multi scale community structure of a graph 
# but not the stability of each scale. 
#
# Inputs : 
# - number of script - must be 2 or 3
# - Nscales, the number of scales to scan between the
# (automatically computed) scale boundaries.
# - A the adjacency matrix (possibly weighted) of the graph. 
# A should be either symmetrical or triangular sup.
# - the number of random vectors Nsignals to use to compute
# the partitions at each scale.
# - lap_flag, the flag that enables to choose between 
# the classical combinatorial Laplacian L=D-A [if lap_flag='lap'], 
# the normalized Laplacian D^(-0.5)LD^(-0.5) [if normlap_flag='normlap'], 
# or the "random walk" Laplacian D^(-1)L [if normlap_flag='rwlap']; 
# - filter_flag is optional. default is filter_flag='wavelet'. But one can
# use filter_flag='scaling_function'if needed.
#
# Outputs :
# - tg is the vector of scales.
# - COM_RV is a NxNscales matrix where the i-th column
# corresponds to the community structure found at scale tg(i).


import numpy as np
import scipy.sparse as ss
import MSCDcreateLaplaciangraph as lp

def DMSCDpartitionRVnoStab(sNum,Nscales,A,T,sVec,Nsignals,lapFlag,filterFlag):
    N=A.shape[0]/T
    # create chosen Laplacian : 
    if lapFlag=='lap':
      L=ss.csr_matrix(lp.MSCDCreateLaplaciangraph(A))
    elif lapFlag=='normlap':
      L=ss.csr_matrix(lp.MSCDCreateNormLaplaciangraph(A))
    elif lapFlag=='rwlap':
      Lnorm=ss.csr_matrix(lp.MSCDCreateNormLaplaciangraph(A))
      L=ss.csr_matrix(lp.MSCDCreateRwLaplaciangraph(A))
    else:
      print('MSCD_partition_RV : lap_flag should be equal to ''lap'', ''normlap'', or ''rwlap''')
    # estimate largest eigenvalue and the three smallest : 
    #opts.isreal=1
    #opts.issym=1
    #opts.tol=1e-3
    if lapFlag=='rwlap':
        #lmax=np.linalg.eigh(Lnorm) #,1,'LA') #,opts) # Lnorm and Lrw have same spectrum. But eigs is 
                                                      # more stable with symmetric matrices                                                     
        D,V=np.linalg.eig(Lnorm.toarray()) #V,D=np.linalg.eigs(Lnorm,(T+2),'SA') #,opts)
    else:
        #lmax=np.linalg.eig(L) #,1,'LA') #,opts) 
        D,V=np.linalg.eig(L.toarray()) #,(T+2),'SA') #,opts)
    lmax=D.shape[0]
    if sNum == 2:
        #D=np.diag(D)
        #Select lambda where small changes are detected
        V1=np.zeros((T,T*N))
        Norm=np.zeros((T))
        beta=np.zeros((T,T))
        for j in range(0,T):
            #V1[j,(1+(j-1)*N):(j*N)]=sVec[j,:]
            V1[j,:]=sVec[j,:]
        for j in range(0,T):
            #[beta(:,j),sigma(j)]=mvregress(V1',V(:,j));
            beta[:,j]=np.linalg.lstsq(V1.transpose(),V[:,j])[0]
            Norm[j]=np.linalg.norm(V1.transpose()*beta[:,j]-V[:,j])
        ind=np.nonzero(Norm>0.9)[0]
        ind=ind[1]
        lambda2=D[ind]
        lambda3=D[ind+1]
    elif sNum == 3:
        lambda2=D[2]
        lambda3=D[3]
    else:
        print('sNum must be 2 or 3')
    bet=1/np.log10(lambda3/lambda2) #this is needed only for 'wavelet' Trembplay2014 filter function
    arange=[0, lmax]
    #Now we can create appropriate scale selection
    # scale interval and filter parameters :
    if filterFlag=='wavelet':
        smin=1/lambda2
        t2=1/lambda2 # this is if one wants to scan very small scales
        smax=t2/lambda2
        tg=np.logspace(np.log10(smin),np.log10(smax),Nscales)
    elif filterFlag=='wavelet2':
        smin=1/lambda2
        t2=1/lambda2
        smax=t2/lambda2
        tg=np.logspace(np.log10(smin),np.log10(smax),Nscales);
    m=100 # Chebychev approximation coefficient
   #cg,g0=MSCDgetInfoG(bet,t2,tg,m,arange,lambda2,lambda3,lmax,smin,smax,filter_flag);
    cg=0
    g0=0
    return tg,cg,g0,arange,L


    
if __name__ == '__main__':           
    print('Hello')
    A=np.array([[1.,2.,3.], [2.,5.,6.],[3.,6.,7.]])
    sVec=np.array([[1.,2.,3.], [1.,4.,6.],[3.,6.,7.]])
    tg,cg,g0,arange,L=DMSCDpartitionRVnoStab(2,3,A,3,sVec,3,'normlap','wavelet')
    pass