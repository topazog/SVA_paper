@echo off

del sva_tm.hds > nul
del sva_tm.bud > nul

ncpar2d_sva < ncpar2d_sva.in > nul

mf6 mfsim.nam > nul

mod2obs_dbl < concsim_.in > nul

mod2obs_dbl < p_concsim_.in > nul
