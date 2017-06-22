#!/bin/sh
#PBS -N 8h4.5s4_NNN
#PBS -q long
#PBS -l nodes=1:ppn=1
#PBS -l walltime=168:00:00
#PBS -k eo

module load openmpi-1.6.4-gcc

cd WORKDIR

sh ../compile_xtc.sh fixed_pwII_NNN
./fixed_pwII_NNN.out /home/ssuruchi/Simulations/Superhydrophobic/Square_pillars/8H4.5S4W2/NNN/traj.xtc>NNN.dat
