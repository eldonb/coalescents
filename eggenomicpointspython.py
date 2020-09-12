import sys
import ctypes
import numpy as np
from scipy import stats
from multiprocessing import Pool
import statistics
import itertools

## compute statistic (X_I, X_J, S) 
##g++ -Wall -fPIC -shared -Ofast -o xibetagenomicpython.so xibetagenomicpython.cpp -lm -lgsl -lgslcblas

def runxi(ah):
    c_double_p = ctypes.POINTER(ctypes.c_double)
    dll_dll = ctypes.CDLL("./genomeg.so")
    xibetagenomic_func = dll_dll.xibetagenomic
    xibetagenomic_func.argtypes=[ctypes.c_double, ctypes.c_double,  c_double_p]
    ## dimension of statistic (Z_I, Z_J, S) = 3
    dimstat = 3
    dataout=(ctypes.c_double * dimstat)()
    trials=10000
    m = np.zeros((trials,dimstat))
    for j in range(trials):
        xibetagenomic_func(ah[0], ah[1],  dataout)
        for i in range(dimstat):
            m[j,i] = dataout[i]
    ## k=stats.gaussian_kde(m.T, bw_method=0.01)
    ## svar= k.pdf(gagnapunktur)
    ## aaa is I = 1, J = 15, ... n-1,  K = 1, ... n-1 and L=20
    ## thiegAB is I = 2...4, J = 15 ...n-2, K = Imin ... Jmax, and L=20
    ## thigenomegBB is excluding singleton and L = 40
    skra="thigenomegBB" + '_' + str(ah[0]) + "_" + str(ah[1]) + "_res.out" 
    f=open(skra, "a")
    for j in range(trials):
        f.write( str(m[j,0]) + ' ' + str(m[j,1]) + ' ' + str(m[j,2]) + "\n")
    f.close()
    del m

valpha=[]
for i in range(0,11,5):
    ## alpha in 1, 1.01, ...,  2.00
    for j in range(200, 510, 10):
        ## theta in  100,200, ..., 1000
        valpha.append([float(i)/10.0, float(j)])

for i in list( itertools.chain( range(2,12, 2), range(20,110, 10), range(200, 1050,50 ))):
    for j in range(200, 510, 10):
        valpha.append([ float(i),  float(j)])


def b(x):
    return x[0] + x[1]

if __name__ == '__main__':
##p=multiprocessing.Pool(processes=6)
    with Pool(processes=50) as p:
        p.map(runxi, valpha)
##    p.close()
del valpha
