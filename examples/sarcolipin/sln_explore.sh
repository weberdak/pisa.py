python ../../pisa.py \
       --peaks assignments.dat \
       --seq MGINTRELFLNFTIVLITVILMWLLVRSYQY \
       --flip 90.0 \
       --rho_start 6 \
       --fit_exclude 1 3 6 7 8 9 10 11 12 13 28 29 31 \
       --fit_tau 20.0 30.0 0.2 \
       --fit_rho0 35.0 55.0 1.0 \
       --fit_order 0.80 1.0 0.02 \
       --explore \
       --out_log sln_explore_log.dat \
       --out_wave sln_explore_wave.dat \
       --out_fit sln_explore_fit.dat