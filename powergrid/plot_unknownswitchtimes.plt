set logscale y
set format y "%g"
set grid
set yrange[1e-8:1e+4]
set title 'optimiter 1, design=0.5: event every 0.5sec'

p \
    'unknownswitchtimes_ndisc3ntime200sstop4.out.cycle'        u 3:5 t 'ndisc3   sstop4   tstop2  ntime200',\
    'unknownswitchtimes_ndisc7ntime400sstop8.out.cycle'        u 3:5 t 'ndisc7   sstop8   tstop4  ntime400',\
    'unknownswitchtimes_ndisc15ntime800sstop16.out.cycle'      u 3:5 t 'ndisc15  sstop16  tstop8  ntime800',\
    'unknownswitchtimes_ndisc31ntime1600sstop32.out.cycle'     u 3:5 t 'ndisc31  sstop32  tstop16 ntime1600',\
    'unknownswitchtimes_ndisc63ntime3200sstop64.out.cycle'     u 3:5 t 'ndisc63  sstop64  tstop32 ntime3200',\
    'unknownswitchtimes_ndisc127ntime6400sstop128.out.cycle'   u 3:5 t 'ndisc127 sstop128 tstop64 ntime6400',\
