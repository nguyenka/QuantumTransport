#!/bin/csh
#PBS -A AFMLW01543HNM
#PBS -M nguyenka
#PBS -l select=1:ncpus=128:mpiprocs=128
#PBS -q high
#PBS -l place=scatter:excl
#PBS -l walltime=2:00:00
#PBS -N cfe_ecp_t19
#PBS -o cfe_ecp_t19.oe
#PBS -e cfe_ecp_t19.oe
cp /p/home/nguyenka/QuantumTransport/examples/cfe_t19/t19_params.json /p/work1/nguyenka/cfe_ecp_t19 
cp /p/home/nguyenka/QuantumTransport/examples/cfe_t19/transport_params.json /p/work1/nguyenka/cfe_ecp_t19 
cd /p/work1/nguyenka/cfe_ecp_t19
t4qt cfe_ecp_t19.xyz -n 128 -e
#qt_dat cfe_ecp_t19.out 
qt cfe_ecp_t19.out -ciss
