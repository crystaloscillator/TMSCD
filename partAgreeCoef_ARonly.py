#! C:\Python27\python

# the following function is used to calculate the 
# adjusted Rand partition comparison coefficients
import numpy as np
import scipy.sparse as ss

def partAgreeCoef_AR(c1,c2):
    c1=c1-np.min(c1)+1
    c2=c2-np.min(c2)+1
    n=np.size(c1); 
    ng1=int(np.max(c1))
    ng2=int(np.max(c2))
    one=np.ones(n)
    confmat=ss.csr_matrix((one, (c1, c2)), shape=(ng1+1, ng2+1)).toarray()
    coltot=np.sum(confmat,0)
    rowtot=sum(confmat.transpose(),0).transpose()
    nis=sum(rowtot**2)  #sum of squares of sums of rows
    njs=sum(coltot**2)  #sum of squares of sums of columns
    t1=nPerm(n,2)		#total number of pairs of entities
    t2=np.sum(np.sum(confmat**2)) #sum over rows & columnns of nij^2
    t3=.5*(nis+njs)
    #Expected index (for adjustment)
    nc=(n*(n**2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1))
    A=t1+t2-t3	#no. agreements
    if t1==nc:
        res=0			#avoid division by zero; if k=1, define Rand = 0
    else:
        res=(A-nc)/(t1-nc)		#adjusted Rand - Hubert & Arabie 1985
   
    return res

#This is the number of combinations of n items taken k at a time
def nPerm(n,k):
    np=1
    for i in range(1,n+1):
        np=np*i
    nk=1
    for i in range(1,n-k+1):
        nk=nk*i
    kp=1
    for i in range(1,k+1):
        kp=kp*i
    return np/(nk*kp)

#npp=nPerm(3,2)
if __name__ == '__main__':
    c1=np.array([2.,2.,2.])
    c2=np.array([2.,2.,2.])
    r=partAgreeCoef_AR(c1,c2)
pass
