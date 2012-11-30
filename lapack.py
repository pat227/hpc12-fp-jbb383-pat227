#!/usr/bin/env python
from numpy import *
from scipy import *

#I confirmed that lapack libary is being used via 
# import numpy.distutils.system_info as sysinfo
# sysinfo.get_info('lapack')
# {'libraries': ['lapack'], 'library_dirs': ['/usr/lib'], 'language': 'f77'}

def main(n):
    example()
    print("Timing QR decompositions for 3x3s")
    starttime = time()
    
    endtime = time()


def example():
    print("Example QR decomposition for a 3x3 matrix...")
    a = array([[12.0,-51.0,4.0],[6.0,167.0,-68.0],[-4.0,24.0,-41.0]])
    b = identity(3, float)
    c = identity(3, float)
    print("A is:")
    print(a)
    print("B is:")
    print(b)
    print("C is:")
    print(c)
    #the key line:
    q, r = linalg.qr(a, mode='full')
    print("QR decomp of A:")
    print("Q:")
    print(q)
    print("R:")
    print(r)






#define function we can time
def qr(a):
    verbose = false
    if(verbose):
        print("A is:")
        print(a)
        #the key line:
    q, r = linalg.qr(a, mode='full')  
    if(verbose):    
        print("QR decomp of A:")
        print("Q:")
        print(q)
        print("R:")
        print(r)



if __name__ == "__main__":
    import sys
    sys.argv[1]    
