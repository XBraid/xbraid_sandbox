from scipy import *
from matplotlib import pyplot as mpl
from os import sys
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from os import system

##
# Run like 
#  $ python viz.py
#
#  Vizualize the files from ex-burgers.c.  The output is assumed to have 
#  the following format.  The filenames are 
#       file_stem + '.' step_number + '.' + rank
#  So, the filename
#       ex-burgers.out.0000350.00000
#  is the output from step 350 from processor rank 0.  Right now,
#  there is no parallelism, so it's always rank 0.
#
#  Each output file has the format
#     
#     ntime_steps
#     tstart
#     tstop
#     nspace_points
#     xstart
#     xstop
#     x[0]
#     x[1]
#       .
#       .
#       .
#     x[k]
##

def update_rcparams(fig_width_pt=700.0, fontsize=32, fontfamily='sans-serif'):
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = 0.7#(sqrt(5)-1.0)/2.25         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]
    
    default_s = 'cm'
    default_ss = 'arial'
    if fontfamily == 'sans-serif':
        to_use = 'arial'
    else:
        to_use = 'cm'

    params = {'backend': 'ps',
              'font.family': fontfamily,    # can specify many more font and math font
              'font.serif':  'cm',          #properties, see ye ole info superhighway  
              'font.sans-serif' : 'arial',
              'axes.labelsize': fontsize,
              'text.fontsize': fontsize,
              'axes.titlesize' : fontsize,
              'legend.fontsize': 24,
              'xtick.labelsize': fontsize,
              'ytick.labelsize': fontsize,
              'text.usetex': True,
              'figure.figsize': fig_size,
              'lines.linewidth' : 3,
              'lines.markersize'  : 10,
              'xtick.major.size':  10, 
              'xtick.minor.size' : 0,
              'ytick.major.size':  10, 
              'ytick.minor.size' : 0,
              # Makes latex all sans-serif, you have two options to do this
              # 'text.latex.preamble' : '\usepackage{sfmath}'}
              # or
              'mathtext.default' : 'regular'} 


    return params
##
params = update_rcparams(fig_width_pt=850., fontfamily='sans-serif')
for key,entry in params.items():
    mpl.rcParams[key] = entry


# Set the braid iteration number and number of steps
rank = 0 
file_stem = 'ex-burgers.out.'

# Find out size of problem in space and time, and the 
# grid spacings
step = 0
fname = file_stem + "%07d"%step + '.' + "%05d"%rank
data = loadtxt(fname)
nsteps = int(data[0])
tstart = float(data[1])
tstop = float(data[2])
nspace = int(data[3])
xstart = float(data[4])
xstop = float(data[5])
mesh = linspace(xstart, xstop, nspace)
data = zeros((nsteps,data.shape[0]-6))

# Load space-time solution
for step in range(nsteps):
    fname = file_stem + "%07d"%step + '.' + "%05d"%rank
    data[step,:] = (loadtxt(fname))[6:]

# Plot IMSHOW 
mpl.figure(0)
mpl.imshow(data,origin='lower',extent=(xstart, xstop, tstart, tstop))
mpl.colorbar()
mpl.ylabel('time')
mpl.xlabel('space')
mpl.xlim([-2,3])
mpl.ylim([0,1.9])
mpl.xticks([-2, 0, 1.5, 3.0]) 
mpl.yticks([0.0, 1.9]) 
mpl.tight_layout(h_pad=0.05)
fname = 'burgers_imshow.pdf'
mpl.savefig(fname, dpi=500, format='pdf')
system('pdfcrop ' + fname + ' ' + fname) 


# Plot 3D Axis 
fig = mpl.figure(1)
ax = fig.gca(projection='3d')
X = linspace(xstart, xstop, nspace)
Y = linspace(tstart, tstop, nsteps)
X, Y = meshgrid(X, Y)
mpl.ylabel('time')
mpl.xlabel('space')
ax.set_xlim([-2,3])
ax.set_ylim([0,1.9])
ax.set_zlim([0.9,2.1])
ax.set_xticks([-2, 0, 2]) 
ax.set_yticks([0.0, 1.9]) 
ax.set_zticks([1, 1.5, 2]) 
ax.plot_trisurf(ravel(X), ravel(Y), ravel(data), linewidth=0.0, antialiased=True, cmap=cm.coolwarm)
ax.view_init(elev=11, azim=-97)
mpl.tight_layout(h_pad=0.05)
fname = 'burgers_3D.png'
mpl.savefig(fname, dpi=160, format='png')

mpl.figure(2)
mpl.plot(mesh, data[0,:], '-o')
mpl.ylabel('u')
mpl.xlabel('space')
mpl.xlim([-2,3])
mpl.ylim([0,2.1])
mpl.xticks([-2, 0, 1.5, 3.0]) 
mpl.yticks([0.0, 1.0, 2.0]) 
mpl.title('Step 0')
fname = 'burgers_step0.pdf'
mpl.savefig(fname, dpi=500, format='pdf')
system('pdfcrop ' + fname + ' ' + fname) 

mpl.figure(3)
mpl.plot(mesh, data[nsteps/2,:], '-o')
mpl.ylabel('u')
mpl.xlabel('space')
mpl.xlim([-2,3])
mpl.ylim([0,2.1])
mpl.xticks([-2, 0, 1.5, 3.0]) 
mpl.yticks([0.0, 1.0, 2.0]) 
mpl.title('Step %d'%(nsteps/2))
fname = 'burgers_step128.pdf'
mpl.savefig(fname, dpi=500, format='pdf')
system('pdfcrop ' + fname + ' ' + fname) 

mpl.figure(4)
mpl.plot(mesh, data[nsteps-1,:], '-o')
mpl.ylabel('u')
mpl.xlabel('space')
mpl.xlim([-2,3])
mpl.ylim([0,2.1])
mpl.xticks([-2, 0, 1.5, 3.0]) 
mpl.yticks([0.0, 1.0, 2.0]) 
mpl.title('Step %d'%(nsteps-1))
fname = 'burgers_step256.pdf'
mpl.savefig(fname, dpi=500, format='pdf')
system('pdfcrop ' + fname + ' ' + fname) 


