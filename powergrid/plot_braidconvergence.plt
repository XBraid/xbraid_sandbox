set logscale y
set format y "%g"
set grid
set title 'excitermodel cf 2, ml 5, sinus, G0=2, tau=0.1'

set xlabel 'braid iteration'
set ylabel 'braid tnorm'


p \
    'braidconvergence.dat' index 0 w l t 'sstop 8' ,\
    'braidconvergence.dat' index 1 w l t 'sstop 16' ,\
    'braidconvergence.dat' index 2 w l t 'sstop 32' ,\
    'braidconvergence.dat' index 3 w l t 'sstop 64' ,\
    'braidconvergence.dat' index 4 w l t 'sstop 128',\
    'braidconvergence.dat' index 5 w l t 'sstop 256',\
    'braidconvergence.dat' index 6 w l t 'sstop 512',\

