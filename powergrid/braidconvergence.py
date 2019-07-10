#!/usr/bin/env python

import os
import sys
import argparse

# build the code
os.system("make")

# set executable 
exe = "./exciter-model-optim"

# get cfactor and maxlevels from command line arguments 
parser = argparse.ArgumentParser()
parser.add_argument("--cfactor", "-cf", help="set coarsening factor, default 2")
parser.add_argument("--maxlevels", "-ml", help="set maxlevels, default 5")
parser.add_argument("--sstopmin", "-smin", help="set smallest sstop, default 8")
parser.add_argument("--sstopmax", "-smax", help="set biggest sstop, default 512")
parser.add_argument("--sstopfactor", "-sf", help="set factor for increasing sstop, 2")
args = parser.parse_args()

if args.maxlevels:
    maxlevels = int(args.maxlevels)
    print("setting maxlevels to ", maxlevels)
else:
    maxlevels = 5
    print("using default cfactor = ", cfactor)
if args.cfactor:
    cfactor = int(args.cfactor)
    print("setting cfactor to ", cfactor)
else:
    cfactor = 2
    print("using default cfactor = ", cfactor)
if args.sstopmin:
    sstopmin = int(args.sstopmin)
    print("setting sstopmin to ", sstopmin)
else:
    sstopmin = 8
    print("using default sstopmin = ", sstopmin)

if args.sstopmax:
    sstopmax = int(args.sstopmax)
    print("setting sstopmax to ", sstopmax)
else:
    sstopmax = 512
    print("using default sstopmax = ", sstopmax)

if args.sstopfactor:
    sstopfactor = int(args.sstopfactor)
    print("setting sstopfactor to ", sstopfactor)
else:
    sstopfactor= 2
    print("using default sstopfactor = ", sstopfactor)


# set dt 
dt = 0.05

out = []
#for i in range(len(sstop)):
sstop = sstopmin
while sstop <= int(sstopmax):
    # get ntime
    ntime = sstop / dt
   
    out.append("#sstop" + str(sstop) + "\n")
    # set upt command string
    cmd = exe  + " -cf "  + str(cfactor) \
               + " -ml " + str(maxlevels) \
               + " -moi 0"               \
               + " -ntime " + str(ntime) \
               + " -sstop " + str(sstop)

    # run the code
    print("\n Running ", cmd, ":\n")
    os.system(cmd)

    # Get convergence from 'braid.out.cycle'
    braidout = open("braid.out.cycle","r")
    lines = braidout.readlines()
    every=2*maxlevels-1
    desiredlines = lines[1::every] # start with 1st, every xth line, x=2*ml -1
    for line in desiredlines:
        out.append(line.split(' ')[4])
        out.append("\n")
    braidout.close()

    out.append("\n\n\n")

    # increase sstop
    sstop = sstop * sstopfactor

# write convergence data to file
filename = "braidconvergence.dat"
print("Writing to ", filename)
output = open(filename, "w")
for item in out:
    output.write(item)
output.close()

