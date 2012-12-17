#!/usr/bin/env python
import numpy
import scipy
import time
import os
import subprocess
#output is suitable for input to gnu plot: (x,y,z) triples that may need sorting
#output is done by the called c - routine itself...this is just like a shell script

def main(args):
    if len(args) != 3:
        print "Usage: python nonlapacktimings.py n rank verbosity"
        print "   n         -> the maximum number of iterations to be used, specified as the exponent of 10^(n-1)"
        print "   rank      -> the upper bound of the square matrix size to decompose as exponent of 2^(r-1)"
        print "   algorithm -> one or more of (a,b,c,d)"
        print " Computes the QR decompositions of random square matrices for a range of iterations and sizes"
        print " using the following optional methods: "
        print "    a)Householder reflectors       b)WY"
        print "    c)Blocked-QR utilizing WY"
#        print " NOTE: on a quad-core Q6600 intel cpu, with rank=50, 100,000 iterations take about 74 seconds."
        return
    n = int(args[0])
    r = int(args[1])
    methods = str(args[2])
    print("Timing QR decompositions...")
    if( n > 7 or n < 1):
        print "N must be a value of 1 through 7 (corresponding to 10^(N-1)); exercise caution in selection."
        return
    if( r > 9 or r < 1):
        print "R must be a value of 1 through 9 (corresponding to 2^(R-1)); exercise caution in selection."
        return
    sizes = [2,4,8,16,32,64,128,256,512,1024]
    iterations = [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000]
    if(methods.find("a") > -1):
        print "Starting iterations up to:", iterations[n], " with sizes up to:", sizes[r]
        for xr in range(r):
            rlocal = sizes[xr]
            subprocess.call(["./householder", str(rlocal), str(rlocal), str(iterations[n-1]), str(0)] )

    if(methods.find("b") > -1):
        print "Starting iterations up to:", iterations[n-1], " with sizes up to:", sizes[r-1]
        for xr in range(r):
            rlocal = sizes[xr]
            for xr2 in range(8):
                rlocal2 = sizes[xr2]
                if(rlocal < 256):
                    subprocess.call(["./wy", str(rlocal), str(rlocal2), str(iterations[n-1]), str(0)] )
                else:
                    pass

    if(methods.find("c") > -1):                
        for xr in range(r):
            rlocal = sizes[xr]
            subprocess.call(["./BlockedQR", str(rlocal), str(rlocal), str(iterations[n-1]), str(0)] )
  
#    if(methods.find("d") > -1):
#        print "Starting iterations up to:", iterations[n], " with sizes up to:", sizes[r]
#        for xr in range(r):
#            rlocal = sizes[xr]
#            subprocess.call(["./BlockedQR2", str(rlocal), str(rlocal), str(iterations[n-1]), str(0)] )

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])    
