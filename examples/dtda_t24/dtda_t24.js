#!/bin/csh
#PBS -A AFMLW01543HNM
#PBS -M nguyenka
#PBS -l select=1:ncpus=128:mpiprocs=128
#PBS -q debug
#PBS -l place=scatter:excl
#PBS -l walltime=1:00:00
#PBS -N dtda_t24
#PBS -o dtda_t24.oe
#PBS -e dtda_t24.oe
#cp /p/home/nguyenka/azoulay/basis_sdd_dz /p/work1/nguyenka/dtda_t24/basis
#cp /p/home/nguyenka/azoulay/auxbasis_universal /p/work1/nguyenka/dtda_t24/auxbasis
cd /p/work1/nguyenka/dtda_t24
#calc_t19.py dtda_t24.xyz
setenv TURBODIR /p/home/nguyenka/TmoleX2024/TURBOMOLE
source $TURBODIR/Config_turbo_env.csh
set EXEC = /p/home/nguyenka/TmoleX2024/TURBOMOLE/smprun_scripts/ridft
#$EXEC -n 128 >& dtda_t24.out
$EXEC -n 128 >& dtda_t24_so.out

