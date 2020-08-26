import linecache
import sys
import ctypes
import numpy as np
from scipy import stats
import multiprocessing
import statistics


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

## num_lines = sum(1 for line in open('myfile.txt'))

def rmle(fskra):
    ## fskra is list of file names with simulated points
    ## check if need this adjustment
    ## sfs now contains random samples of (Z_I, Z_J, S) with S > 5000
    modelfiles=[]
    with open(fskra, "r") as f:
        for l in f:
            modelfiles.append(l)
    f.close()
    ## run over modelfiles and compute density for each
    ##m = np.zeros((10000,3))
    ## there is only one point  for i in range(len(sfs)):
    ## gagnapunktur=sfs[i]
    gagnapunktur=[ 0.4313527, 0.3129360 ]
    mlea = 0.0
    mleh=0.0
    lmax=0.0
    l=0.0
    for z in modelfiles:
        ##r =  slod + z
        m=np.loadtxt(z.strip())
        k=stats.gaussian_kde(m.T, bw_method=0.01)
        del m
        l=float(k.pdf(gagnapunktur))
        del k
        r=z.strip()
        alpha = float( (r.partition("_")[2]).partition("_")[0])
        theta= float(  ((r.partition("_")[2]).partition("_")[2]).partition("_")[0] )
        mlea = alpha if l > lmax else mlea
        mleh = theta if l > lmax else mleh
        lmax = l if l > lmax else lmax
    f=open("mleresoutBBB", "a")
    f.write( str([mlea, mleh,   np.log( lmax if lmax > 0.0 else np.exp(-100.0))]) + "\n")
    f.close()

##rmle(1, 180, 5000, "thisthomolegA/eggfiles", "/home/bjarki/verk/mfn/truncation/wfiles/testcode/thisthomolegA/")
## check lagmark when excluding singletons and anti-singletons
## count lines in file with data points
##num_lines = sum(1 for line in open('data'))
rmle("./aaafiles")
