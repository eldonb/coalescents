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
## setwd('/home/bjarki/verk/mfn/truncation/wfiles/testcode/')
## thiaaagenomicfragments.txt

## mm <- alllgpunktar(180, 20,  5000)

def rmle(gskra, fskra):
    ## check if need this adjustment
    sfs=[]
    with open(gskra) as f:
        for l in f:
          sfs.append([float(x) for x in l.split()])
    f.close()
    ## sfs now contains random samples of (Z_I, Z_J, S) with S > 5000
    gl= len(sfs)
    sennileiki = []
    alpha=[]
    theta=[]
    modelfiles=[]
    with open(fskra, "r") as f:
        for l in f:
            modelfiles.append(l)
    f.close()
    mleskra=""
    maxsennileiki = 0.0
    d=0.0
    ## run over modelfiles and compute density for each
    ##m = np.zeros((10000,3))
    for i in range(len(sfs)):
        gagnapunktur=sfs[i]
        mleskra=""
        d=0.0
        maxsennileiki=0.0
        for z in modelfiles:
            ##r =  slod + z
            m=np.loadtxt(z.strip())
            k=stats.gaussian_kde(m.T, bw_method=0.01)
            del m
            d=float(k.pdf(gagnapunktur))
            mleskra=z.strip() if d > maxsennileiki else mleskra
            maxsennileiki= d if d > maxsennileiki else maxsennileiki
            ## collect the alpha and theta from
            if mleskra.count("") > 1:
                ## add statistics for genomic point i
                alpha.append( float( (mleskra.partition("_")[2]).partition("_")[0]) )
                theta.append( float(  ((mleskra.partition("_")[2]).partition("_")[2]).partition("_")[0] ) )
                sennileiki.append( np.log( maxsennileiki if maxsennileiki > 0.0 else np.exp(-100.0) ) )
    f=open("thiaaamlxibetaresout", "a")
    f.write("genomic points statistics :\n" )
    if len(alpha) > 1:
        f.write(str([statistics.mean(alpha), statistics.stdev(alpha), statistics.mean(theta), statistics.stdev(theta), statistics.mean(sennileiki)]) + "\n" )
    else:
        f.write( str([alpha, theta, sennileiki])  + "\n")
    f.write( "close output for lgroup " + str(j) + "\n" )
    f.close()

##rmle(1, 180, 5000, "thisthomolegA/eggfiles", "/home/bjarki/verk/mfn/truncation/wfiles/testcode/thisthomolegA/")
## check lagmark when excluding singletons and anti-singletons
## count lines in file with data points
##num_lines = sum(1 for line in open('data'))
##for a in range(1,24):
##    rmle(a, 180,  "aaafiles", "./")
rmle("./thiaaagenomicfragments.txt", "./aaafiles")
