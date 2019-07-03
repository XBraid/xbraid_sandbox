set logscale y
set format y "%g"
set grid
set title 'excitermodel with limits, ml5, cf 2'

set xlabel 'braid iteration'
set ylabel 'braid tnorm'

p \
    'excitermodel-limits-sstop10-ntime200.braid.out.ml5'    u 3:5 t 'sstop10' ,\
    'excitermodel-limits-sstop20-ntime400.braid.out.ml5'    u 3:5 t 'sstop20' ,\
    'excitermodel-limits-sstop40-ntime800.braid.out.ml5'    u 3:5 t 'sstop40' ,\
    'excitermodel-limits-sstop80-ntime1600.braid.out.ml5'   u 3:5 t 'sstop80' ,\
    'excitermodel-limits-sstop160-ntime3200.braid.out.ml5'  u 3:5 t 'sstop160',\
    'excitermodel-limits-sstop320-ntime6400.braid.out.ml5'  u 3:5 t 'sstop320',\

