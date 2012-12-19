#!/usr/bin/env python
#import numpy
#import scipy
#import time
import os
import subprocess
#output is suitable for input to gnu plot: (x,y,z) triples that may need sorting
#output is done by the called c - routine itself...this is just like a shell script

def main(args):
    if len(args) != 3:
        print "Usage: python nonlapacktimings.py n rank verbosity"
        print "   n         -> the maximum number of iterations to be used, specified as the "
        print "                exponent of 10^(n-1)"
        print "   rank      -> the upper bound of the square matrix size to decompose as "
        print "                (almost) an exponent of 2^(r-1); we with a few cases of "
        print "                non-powers of 2 thrown in at 100, 500, and 1000."
        print "   algorithm -> one or more of (a,b,c,d)"
        print " Computes the QR decompositions of random square matrices for a range of "
        print " iterations and sizes using the following optional methods: "
        print "    a)Householder reflectors       b)WY"
        print "    c)Blocked-QR utilizing WY      d)BlockedQR2 "
        print "    e)Blocked-QR using optimized blocksize=8"
        print "    f)Blocked-QR2 using optimized blocksize=8"
        print "    g)Scaling measurements (different #s of OMP threads) for e (n and r fixed)"
        print "    h)Scaling measurements (different #s of OMP threads) for f (n and r fixed)"
        return
    n = int(args[0])
    r = int(args[1])
    methods = str(args[2])
    print("Timing QR decompositions...")
    if( n > 8 or n < 1):
        print "N must be a value of 1 through 8 (corresponding to 10^(N-1)); exercise caution in selection."
        return
    if( r > 14 or r < 1):
        print "R must be a value of 1 through 13 (almost corresponding to 2^(R-1)); exercise caution in selection."
        return
    sizes = [2,4,8,16,32,64,100,128,256,500,512,1000,1024,2048]
    iterations = [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000]

    if(methods.find("a") > -1):
        print "Starting iterations up to:", iterations[n-1], " with sizes up to:", sizes[r-1]
        for xr in range(r):
            rlocal = sizes[xr]
            subprocess.call(["/home/pat227/hpc-fall12/hpc12-proj-pat227-jbb383/QR/householder", str(rlocal), str(rlocal), str(iterations[n-1]), str(0)] )

    if(methods.find("b") > -1):
        print "Starting iterations up to:", iterations[n-1], " with sizes up to:", sizes[r-1]
        for xr in range(r):
            rlocal = sizes[xr]
            subprocess.call(["/home/pat227/hpc-fall12/hpc12-proj-pat227-jbb383/QR/wy", str(rlocal), str(rlocal), str(iterations[n-1]), str(0)] )

    if(methods.find("c") > -1):
        print "Starting iterations up to:", iterations[n-1], " with sizes up to:", sizes[r-1]
        for xr in range(r):
            rlocal = sizes[xr]
            for xr2 in range(13):
                rlocal2 = sizes[xr2]
                #this is an openmp sub call; be sure that environment variable OMP_NUM_THREAD=8 in the pbs script
#                environ = os.environ()
#                environ["OMP_NUM_THREADS"]="8" 
                subprocess.call(["/home/pat227/hpc-fall12/hpc12-proj-pat227-jbb383/QR/BlockedQR", str(rlocal), str(rlocal2), str(iterations[n-1]), str(0)])
  
    if(methods.find("d") > -1):
        print "Starting iterations up to:", iterations[n-1], " with sizes up to:", sizes[r-1]
        for xr in range(r):
            rlocal = sizes[xr]
            for xr2 in range(13):
                rlocal2 = sizes[xr2]
                #this is an openmp sub call; be sure that environment variable OMP_NUM_THREAD=8 in the pbs script
#                environ = os.environ()
#                environ["OMP_NUM_THREADS"]="8" 
                subprocess.call(["/home/pat227/hpc-fall12/hpc12-proj-pat227-jbb383/QR/BlockedQR2", str(rlocal), str(rlocal2), str(iterations[n-1]), str(0)])

    if(methods.find("e") > -1):
        print "Starting iterations up to:", iterations[n-1], " with sizes up to:", sizes[r-1]
        for xr in range(r):
            rlocal = sizes[xr]
            for xr2 in range(13):
                rlocal2 = sizes[xr2]
                #this is an openmp sub call; be sure that environment variable OMP_NUM_THREAD=8 in the pbs script
#                environ = os.environ()
#                environ["OMP_NUM_THREADS"]="8" 
                subprocess.call(["/home/pat227/hpc-fall12/hpc12-proj-pat227-jbb383/QR/BlockedQR_8", str(rlocal), str(rlocal2), str(iterations[n-1]), str(0)])
  
    if(methods.find("f") > -1):
        print "Starting iterations up to:", iterations[n-1], " with sizes up to:", sizes[r-1]
        for xr in range(r):
            rlocal = sizes[xr]
            for xr2 in range(13):
                rlocal2 = sizes[xr2]
                #this is an openmp sub call; be sure that environment variable OMP_NUM_THREAD=8 in the pbs script
#                environ = os.environ()
#                environ["OMP_NUM_THREADS"]="8" 
                subprocess.call(["/home/pat227/hpc-fall12/hpc12-proj-pat227-jbb383/QR/BlockedQR2_8", str(rlocal), str(rlocal2), str(iterations[n-1]), str(0)])

    if(methods.find("g") > -1):
        print "Starting scaling measurements for BlockedQR..."
        for x in range(8):
            rlocal = 2048
            rlocal2 = 2048
            #this is an openmp sub call; be sure that environment variable OMP_NUM_THREAD=x is set here 
            environ = os.environ
            environ["OMP_NUM_THREADS"]=str(8-x) 
            subprocess.call(["/home/pat227/hpc-fall12/hpc12-proj-pat227-jbb383/QR/BlockedQR_8_scaled", str(rlocal), str(rlocal2), str(iterations[n-1]), str(0)])
  
    if(methods.find("h") > -1):
        print "Starting scaling measurements for BlockedQR2..."
        for x in range(8):
            rlocal = 2048
            rlocal2 = 2048
            #this is an openmp sub call; be sure that environment variable OMP_NUM_THREAD=x is set here
            environ = os.environ
            environ["OMP_NUM_THREADS"]=str(8-x) 
            subprocess.call(["/home/pat227/hpc-fall12/hpc12-proj-pat227-jbb383/QR/BlockedQR2_8_scaled", str(rlocal), str(rlocal2), str(iterations[n-1]), str(0)])

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])    
