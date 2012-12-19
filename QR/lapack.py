#!/usr/bin/env python
import numpy
#from scipy import *
import time
#I confirmed that lapack libary is being used via 
# import numpy.distutils.system_info as sysinfo
# sysinfo.get_info('lapack')
# {'libraries': ['lapack'], 'library_dirs': ['/usr/lib'], 'language': 'f77'}

def main(args):
    if len(args) != 3:
        print "Usage: python lapack.py n rank"
        print "   n -> # of times to compute QR decomposition"
        print "   rank -> the size of the square matrices to decompose"
        print "            as 2^(r-1)"
        print "   verbosity -> 0,1"
        return
    n = int(args[0])
    r = int(args[1])
    verbose = int(args[2])
    if(n < 1):
        print "n must be > 0; exercise caution"
    if(r < 1):
        print "rank must be > 0; exercise caution"
#    example()
    print("Timing QR decompositions...")
    sizes = [2,4,8,16,32,64,128,256,512,1024,2048]
    file = open("lapack_time.txt", 'a')
    file2 = open("lapack_gflps.txt", 'a')
    gflps = 0.0
    for x in range(n):
        file.write('\n')
        file2.write('\n')
        for s in range(r):
            file.write('\n')
            file2.write('\n')
            for s2 in range(r):
                size = sizes[s]
                size2 = sizes[s2]
                starttime = time.time()
                a = numpy.random.randn(size, size2)
                qr(a, verbose)
                endtime = time.time()
                elapsed = (endtime-starttime)
                print "m:", size, "n:", size2, "time:", elapsed,
                file.write('\n')
                file.write(str(size))
                file.write(' ')
                file.write(str(size2))
                file.write(' ')
                file.write(str(elapsed))
#                print "QR decomps / s: ", (r/elapsed)
                gflps = size * size2 * size2 / elapsed / 1000000000
                print "Gflops /s: ", gflps
                file2.write('\n')
                file2.write(str(size))
                file2.write(' ')
                file2.write(str(size2))
                file2.write(' ')
                file2.write(str(gflps))
    file.close()
    file2.close()

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
    q, r = numpy.linalg.qr(a, mode='full')
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
    q, r = numpy.linalg.qr(a, mode='full')  
    if(v == 1):    
        print "QR decomp of A:"
        print "Q:"
        print(q)
        print "R:"
        print(r)


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])    
