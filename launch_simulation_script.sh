ne#!/bin/sh
#PBS -W group_list=cu_10184 -A cu_10184
#PBS -m n
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -l walltime=7:00:00:00

module load tools
module load gcc
module load intel/perflibs
module load R/4.1.0

Rscript /home/people/morytu/scripts/EVI2_EPICarray/simulation_script_general.R $1 $2 $3 $4 $5 $6
