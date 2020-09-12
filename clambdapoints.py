import sys
import ctypes
import numpy as np
from scipy import stats
from multiprocessing import Pool
import statistics
import itertools


##constants.hpp  hjhj      thigchgl2lg03grisjad_resout0  thigchgl2lg4.out       thilg4_resout0            xibetahomolmergerspython.hpp  xibetapython.so
##BJARKI@MFN0449:~/verk/sweepstakes> less thigchgl2lg03grisjad_resout0 | tee
##[1.1136842105263158, 0.06943849583796058, 333.7719298245614, 39.56377960230557, 0.28990510916300405]BJARKI@MFN0449:~/verk/sweepstakes> 
##BJARKI@MFN0449:~/verk/sweepstakes> 

## compute statistic (X_I, X_J, S) 
##g++ -Wall -fPIC -shared -Ofast -o xibetagenomicpython.so xibetagenomicpython.cpp -lm -lgsl -lgslcblas

def runxi(ah):
    c_double_p = ctypes.POINTER(ctypes.c_double)
    dll_dll = ctypes.CDLL("./clambda.so")
    xibetagenomic_func = dll_dll.lambdaijk
    xibetagenomic_func.argtypes=[ctypes.c_double, ctypes.c_double,  c_double_p]
    ## dimension of statistic (Z_I, Z_J, S) = 3
    dimstat = 4
    dataout=(ctypes.c_double * dimstat)()
    trials=3
    m = np.zeros((trials,dimstat))
    for j in range(trials):
        xibetagenomic_func(ah[0], ah[1],  dataout)
        for i in range(dimstat):
            m[j,i] = dataout[i]
    ## k=stats.gaussian_kde(m.T, bw_method=0.01)
    ## svar= k.pdf(gagnapunktur)
    ## bbb for betaxi, egg for exponential growth
    skra="omegalambdabeta" + '_' + str(ah[0]) + "_" + str(ah[1]) + "_res.out" 
    f=open(skra, "a")
    for j in range(trials):
        f.write( str(m[j,0]) + ' ' + str(m[j,1]) + ' ' + str(m[j,2]) + "\n")
    f.close()
    del m

valpha=[]
for i in range(11):
    ## alpha in 1, 1.01, ...,  2.00
    ## mle  theta = 44 for south-coast; mle theta = 31 for thistilfj
    for j in range(40, 51, 1):
        ## theta in  100,200, ..., 1000
        valpha.append([1.0 + float(i)/100.0, float(j)])

def b(x):
    return x[0] + x[1]

if __name__ == '__main__':
##p=multiprocessing.Pool(processes=6)
    with Pool(processes=10) as p:
        p.map(runxi, valpha)
##    p.close()
del valpha
