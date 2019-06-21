set logscale y
set format y "%g"
set grid
set yrange [1e-8:1e+3]

p \
'aAsSeparateVar2_tstop4ntime200.out.cycle'     u 3:5 t 'ndisc5   tstop4  ntime200',\
'aAsSeparateVar2_tstop8ntime400.out.cycle'     u 3:5 t 'ndisc10  tstop8  ntime400',\
'aAsSeparateVar2_tstop16ntime800.out.cycle'    u 3:5 t 'ndisc22  tstop16 ntime800',\
'aAsSeparateVar2_tstop32ntime1600.out.cycle'   u 3:5 t 'ndisc44  tstop32 ntime1600',\
'aAsSeparateVar2_tstop64ntime3200.out.cycle'   u 3:5 t 'ndisc88  tstop64 ntime3200',\
