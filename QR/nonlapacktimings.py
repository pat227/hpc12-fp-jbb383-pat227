#!/usr/bin/env python
import numpy
import scipy
import time
import os
import subprocess

def main(args):
    if len(args) != 3:
        print "Usage: python nonlapacktimings.py n rank verbosity"
        print "   n         -> the maximum number of iterations to be used, specified as the exponent of 10^n"
        print "   rank      -> the upper bound of the square matrix size to decompose as exponent of 2^(r+1)"
        print "   verbosity -> 0,1"
        print " Computes the QR decompositions of random square matrices for a range of iterations and sizes"
        print " using a variety of methods: "
        print "    a) Householder reflectors      b)WY algo w/o blocking"
        print "    c)Blocked-QR upon method 'b'   d)Blocked-QR upon WY w/blocking"
#        print " NOTE: on a quad-core Q6600 intel cpu, with rank=50, 100,000 iterations take about 74 seconds."
        return
    n = int(args[0])
    r = int(args[1])
    verbose = int(args[2])
    #example()
    print("Timing QR decompositions...")
    if( n > 7 or n < 0):
        print "N must be a value of 1 through 7 (corresponding to 10^n); exercise caution in selection."
        return
    if( r > 9 or r < 0):
        print "R must be a value of 0 through 9 (corresponding to 2^1 - 2^10)"
        return
    sizes = [2,4,8,16,32,64,128,256,512,1024]
    iterations = [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000]
    elapsed = 0.0;
    print "Starting iterations up to:", iterations[n], " with sizes up to:", sizes[r]
    for xr in range(r):
        rlocal = sizes[xr]
        sizeargs = str(rlocal) + " " + str(rlocal)
        print "Using args:" + sizeargs
        for xn in range(n):
            starttime = time.time()
            for x in range(iterations[xn]):
                subprocess.call(["wy", "3 3"])
            endtime = time.time()
            elapsed = (endtime-starttime)
            print "Matrix r:", rlocal, "Iters:", iterations[xn], "Time:", elapsed, "s  WY decomps/s:", (iterations[xn] / elapsed), "GB/s:", (rlocal * rlocal * 8 * iterations[xn] / elapsed / 1000000000)
    #output suitable for input to gnu plot (x,y,z) triples that will need sorting
    #need multiple files for different series (iterations,time), (size,GB/s)
            filename1 = "WY_time.txt"
            outfile1 = open(filename1, 'a+')
    #outfile.write("#Matrix size:" + str(r) + "^2\n")
    #outfile.write("#Iterations, seconds \n")
            outfile1.write(str(rlocal) + " " + str(xn) + " " + str(elapsed) + "\n")
            filename2 = "WY_gb.txt"
            outfile2 = open(filename2, 'a+')
    #outfile.write("#Matrix size:" + str(r) + "^2\n")
    #outfile.write("#GBs, seconds \n")
            outfile2.write(str(rlocal) + " " + str(xn) + " " + str(rlocal * rlocal * 8 / elapsed / 1000000000) + " \n")
    outfile1.close()
    outfile2.close()

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])    
