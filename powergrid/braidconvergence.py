#!/usr/bin/env python

import os

# build the code
os.system("make")

# set executable 
exe = "./exciter-model-optim"

# set arguments
cfactor      = 2
maxoptimiter = 0

# set range for final time T and number of time steps
dt = 0.05
sstop = [  8,  16,  32,   64,  128,  256,   512]
ntime = [1600, 3200, 6400, 12800, 25600, 51200, 102400]


out = []
for i in range(len(sstop)):
   
    out.append("#sstop" + str(sstop[i]) + "\n")
    # set upt command string
    cmd = exe  + " -cf "  + str(cfactor) \
               + " -moi " + str(maxoptimiter) \
               + " -ntime " + str(ntime[i]) \
               + " -sstop " + str(sstop[i])

    # run the code
    print("\n Running ", cmd, ":\n")
    os.system(cmd)

    # Get convergence from 'braid.out.cycle'
    braidout = open("braid.out.cycle","r")
    lines = braidout.readlines()
    desiredlines = lines[1::9]
    for line in desiredlines:
        out.append(line.split(' ')[4])
        out.append("\n")
    braidout.close()

    out.append("\n\n\n")

# write convergence data to file
filename = "braidconvergence.dat"
print("Writing to ", filename)
output = open(filename, "w")
for item in out:
    output.write(item)
output.close()

