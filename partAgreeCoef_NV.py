#! C:\Python27\python

# the following function is used to calculate the 
#  partition comparison coefficients

import numpy as np
import scipy.sparse as ss

def partAgreeCoef_NV(c1,c2):
    c1=c1-np.min(c1)+1
    c2=c2-np.min(c2)+1
    n=np.size(c1); 
    ng1=int(np.max(c1))
    ng2=int(np.max(c2))
    one=np.ones(n)
    confmat=ss.csr_matrix((one, (c1, c2)), shape=(ng1+1, ng2+1)).toarray()
    coltot=np.sum(confmat,0)
    rowtot=sum(confmat.transpose(),0).transpose()
    H1=-np.sum((rowtot/n)*np.log2((rowtot/n)))
    H2=-np.sum((coltot/n).transpose()*np.log2((coltot/n).transpose()))
    indmat=(rowtot/n)*(coltot/n)
    nozeromat=(confmat/n)+(confma==0)
    H12=-np.sum(np.sum((confmat/n)*np.log2(nozeromat)))
    MI=H1+H2-H12
    VI=H1+H2-2*MI
    res=VI/log2(n)

