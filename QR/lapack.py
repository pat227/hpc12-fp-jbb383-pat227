#!/usr/bin/env python
from numpy import *
from scipy import *
from time import *
#I confirmed that lapack libary is being used via 
# import numpy.distutils.system_info as sysinfo
# sysinfo.get_info('lapack')
# {'libraries': ['lapack'], 'library_dirs': ['/usr/lib'], 'language': 'f77'}

def main(args):
    if len(args) != 3:
        print "Usage: python lapack.py n rank"
        print "   n -> # of time to compute QR decomposition"
        print "   rank -> the size of the square matrices to decompose"
        print "   verbosity -> 0,1"
        return
    n = int(args[0])
    r = int(args[1])
    verbose = int(args[2])
#    example()
    print("Timing QR decompositions...")
    starttime = time()
    for x in range(n):
        a = random.randn(r, r)
        qr(a, verbose)
    endtime = time()
    elapsed = (endtime-starttime)
    print "Total Time elapsed:", elapsed, " s"
    print "QR decomps / s: ", (n/elapsed)
    print "Gflops /s: ", (n*r*r/elapsed)

def example():
    print("Example QR decomposition for a 3x3 matrix...")
    a = array([[12.0,-51.0,4.0],[6.0,167.0,-68.0],[-4.0,24.0,-41.0]])
    b = identity(3, float)
    c = identity(3, float)
    print "A is:"
    print(a)
    print "B is:"
    print(b)
    print "C is:"
    print(c)
    #the key line:
    q, r = linalg.qr(a, mode='full')
    print "QR decomp of A:"
    print "Q:"
    print(q)
    print "R:"
    print(r)


#define function we can time
def qr(a, v):
    if(v == 1):
        print "A is:"
        print(a)
        #the key line:
    q, r = linalg.qr(a, mode='full')  
    if(v == 1):    
        print "QR decomp of A:"
        print "Q:"
        print(q)
        print "R:"
        print(r)


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])    
