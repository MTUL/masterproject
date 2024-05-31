#!/bin/sh
#PBS -W group_list=cu_10184 -A cu_10184
#PBS -N perm
#PBS -e perm.err
#PBS -o perm.log
#PBS -m n
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -l walltime=1:00:00

# 1: samplesize - Number of patients in simulated dataset (0 = use original dataset)
# 2: n_cpg - Number of CpG sites in simulated dataset (0 = use original dataset)
# 3: effect - Simulated treatment effect to be used (index)
# 4: nsim - Number of simulated datasets per simulation
# 5: nperm - Number of permutations per resampling procedure
# 6: res.dist - Residual distribution: 0: normal, 1: empirical residuals

# Statistical power for different g functions with either empirical or normal residuals in original dataset:
nsim="1000"
nperm="1000"
samplesize="0"
n_cpg="0"
for eff in {1..8}; do
      for resdist in {0..1}; do
            qsub -F "$samplesize $n_cpg $eff $nsim $nperm $resdist" /home/people/morytu/scripts/EVI2_EPICarray/launch_simulation_script.sh
            sleep 5
      done
done

# Investigating alternative method for CpG-specific treatment effects
nsim="1000"
nperm="1000"
samplesize="0"
n_cpg="720"
for eff in {1..8}; do
            for resdist in {0..1}; do
                  qsub -F "$samplesize $n_cpg $eff $nsim $nperm $resdist" /home/people/morytu/scripts/EVI2_EPICarray/launch_simulation_script.sh
            sleep 5
      done
done

# For resdist = empirical residuals, test how inclusion of additional CpG sites affects statistical power:
nsim="1000"
nperm="1000"
samplesize="0"
resdist="1"
for n_cpg in "10" "50" "100" "300" "1500" "3000" "5000" "10000"; do
      for eff in {1..8}; do
            qsub -F "$samplesize $n_cpg $eff $nsim $nperm $resdist" /home/people/morytu/scripts/EVI2_EPICarray/launch_simulation_script.sh
      sleep 5
      done
done

# Stratified simulations:
n_cpg="50"
nsim="1000"
nperm="1000"
for ptsubset in "mut" "WT"; do
      for eff in {1..8}; do
            qsub -F "$ptsubset $n_cpg $eff $nsim $nperm" /home/people/morytu/scripts/EVI2_EPICarray/launch_simulation_script_stratified.sh
            sleep 5
      done
done























