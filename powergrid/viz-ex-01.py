from scipy import loadtxt, ones_like
import matplotlib as mpl
import matplotlib.pyplot as plt
from os import system

def update_rcparams(fig_width_pt=700.0, fontsize=32, fontfamily='sans-serif'):
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = 0.7#(sqrt(5)-1.0)/2.25         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]

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
              'lines.linewidth' : 6,
              'lines.markersize'  : 13,
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


# Load data
savefiggers = True
system('cat ex-01-expanded.out.* > temp.out')
data = loadtxt('temp.out')

# Plot data
twos = 2.0*ones_like(data)
plt.plot(twos, '--b')
plt.plot(0.5*twos, '--b')
plt.plot(data, '-k')
plt.xlabel('Time Step')
plt.ylabel('Solution')
plt.title('Thermostat Problem as Nonlinear ODE')

if savefiggers:
    plt.tight_layout(h_pad=0.05)
    fname = 'thermostat_solution.pdf'
    # if using PNG, can set transparent=True in savefig
    plt.savefig(fname, dpi=1000, format='pdf')
    system('pdfcrop ' + fname + ' ' + fname) 

# Clean up
system('rm temp.out')

