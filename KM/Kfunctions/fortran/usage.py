import fortran_k as FKM
#import _kernels as CKM
import timeit
import time
import sys
import math as mth
import numpy as np


def py_Bernoulli(x,degre):
    """
    Trying to solve BN1 problem
    """
    y=None
    if degre == 0: y = 1.
    if degre == 1: y = x - 0.5
    if degre == 2: y = x * x - x + 1./6.
    if degre == 3: y = x * (x * (x - 1.5) + 0.5)
    if degre == 4: y = x * (x * (x * (x - 2.0) + 1.0)) - 1./30.
    if degre == 5: y = x * (x * (x * (x * (x - 2.5) + 5./3.)) - 1./6.)
    if degre == 6: y = x * (x * (x * (x * (x * (x - 3.0) + 2.5)) - 0.5)) + 1./42.
    return y
def py_KernelBernoulli(x, y, degre):
    #double v, s, tgamma(double), pow(double, double), floor(double);
    s=0
    for d in range (degre+1):
        s = s + py_Bernoulli(x, d) * py_Bernoulli(y, d) / mth.factorial(d)**2
        #print 'py',d,s,mth.factorial(d)**2
    #return s
    v=abs(x-y)# parte entera o abs?
    #v=x-y-mth.floor(x-y)
    #print 'coso py', mth.factorial(2. * degre )
    return s + ( ((-1)**(degre+1.)) /mth.factorial(2. * degre ))* py_Bernoulli(v, 2 * degre)

def K_prod(x,y,k,kernel,p=1):
    for xi,yi,ki in zip(x,y,k):
        #print zip(x,y,k)
        #sys.exit()
        p=p*kernel(xi,yi,ki)
    return p

def calculK(X,k,Kprod,kernel,regul):
    K=np.zeros([len(X),len(X)],order="C")
    for i,x in enumerate(X):
        for j,y in enumerate(X):
            K[i][j]=Kprod(x,y,k,kernel)
    K=K +regul*np.eye(len(K))
    #invK=np.linalg.inv(K)
    return K

def predict(k,x,X,kernel,alpha,norm_xsmean, xsmax,Kprod):
    """
    Making the predcition a static method allows to use it outside the class
    """
    ret_val=np.inner(alpha,[Kprod(x,xi,k,kernel) for xi in X])
    return (ret_val+norm_xsmean)*xsmax
x=0.2
y=0.3
degre=int(3)
#print 'here', FKM.bernoulli(x,degre)
#print py_KernelBernoulli(x, y, degre)
#print py_KernelBernoulli(x, y, degre)
#print FKM.kprod_bernoulli(x, y, degre)
# sys.exit()

x=np.array((0.5,0.3,0.25),order='F')
y=np.array((0.5,0.4,0.1),order='F')
k=np.array((2,2,2),order='F')

print 'starting'
N=1E2
X=[x for _ in range(int(N))]
X=np.array(X,order='F')
alpha=[1.0 for _ in range(int(N))]
alpha=np.array(alpha,order='F')
#print FKM.predict2(alpha,int(len(alpha)),x,k,int(len(k)))
print y ,int(len(y))
print 'ciao antonio'
print FKM.predict3(alpha,int(len(alpha)))
print FKM.predict2(y,k,int(len(y)))
print FKM.predict4(alpha,y,k,int(len(y)),int(len(alpha)))
print int(len(y)),int(len(alpha))
print FKM.predict5(alpha,y,k,3,100)
sys.exit()
#kernel=FKM.kernel_bernoulli
# print calculK(X,k,K_prod,kernel,0)
# print FKM.calcul_k('KBN',len('KBN'),X,k,0.0,int(len(k)),int(len(X)))

NN=int(1E2)
beg=time.time()
i=0
for _ in range(NN):
    [K_prod(x,y,k,kernel) for x in  X]
    #calculK(X,k,K_prod,kernel,0)
    #np.linalg.pinv(calculK(X,k,K_prod,kernel,0))
pyt=time.time()-beg

# beg=time.time()
# i=0
# while i<N:
#     CKM._KernelBernoulli(x, y, degre)
#     i=i+1
# ct=time.time()-beg

ft=0
beg=time.time()
for _ in range(NN):
    [FKM.kprod_spline(x,y,k,len(x)) for x in  X]
    #FKM.calcul_k('KBN',len('KBN'),X,k,0.0,int(len(k)),int(len(X)))
    #np.linalg.pinv(FKM.calcul_k('KBN',len('KBN'),X,k,0.0,int(len(k)),int(len(X))))
ft=time.time()-beg

# ftv=0
# beg=time.time()
# P=FKM.vect_kprod_spline(X, y, k, len(y) ,len(X))
# ftv=time.time()-beg

# x=np.array((0.1,0.3,0.25),order='F')
# y=np.array((0.5,0.4,0.1),order='F')
# k=np.array((2,2,2),order='F')

#X=[x]

# print K_prod(x,y,k,CKM._KernelSpline),K_prod(x,y,k,CKM._KernelBernoulli)
# print FKM.vect_kprod_bernoulli('KSP',X,y,k, int(len(y)) ,int(len(X)))[0]
# print FKM.vect_kprod_bernoulli('KBN',X,y,k, int(len(y)) ,int(len(X)))[0]




print pyt,ft
# print timeit.timeit(CKM._Bernoulli(x,degre))
print 'sucess'