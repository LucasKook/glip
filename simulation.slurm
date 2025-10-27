#!/bin/bash
#SBATCH --job-name=glip # short name for your job
#SBATCH --mail-user=lkook@wu.ac.at
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o %x_%j_%N.out
#SBTACH -e %x_%j_%N.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=short

echo "Running job $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
echo "Allocated hosts: $SLURM_JOB_NODELIST"

. /opt/apps/2023-04-11_lmod.bash
module purge
module load r-4.3.2-gcc-12.2.0-cx3fkbv
module load gmp-6.2.1-gcc-12.2.0-72i3jt5
module load proj-8.2.1-gcc-12.2.0-aoxmoc3
module load freetype-2.11.1-gcc-12.2.0-3fhiifd
module load mpfr-4.1.0-gcc-12.2.0-blxjshv
module load sqlite-3.39.4-gcc-12.2.0-2ja6wxv
module load libjpeg-9e-gcc-12.2.0-zipplqc
module load miniconda3-23.11.0-gcc-12.2.0-whk7wun
module load geos-3.9.1-gcc-12.2.0-y43p3ma
module load libtiff-4.4.0-gcc-12.2.0-sdt4m5m
module load udunits-2.2.28-gcc-12.2.0-ihgdna5
module load icu4c-73.2-gcc-12.2.0-hroicyi
module load cmake-3.24.3-gcc-12.2.0-uzbdae4
module load gdal-3.5.3-gcc-12.2.0-fduwdss
module load cuda-11.8.0-gcc-12.2.0-w5yxkzl

export R_LIBS=/home/lkook/R_LIBS
export GUROBI_HOME=/home/lkook/gurobi1202/
export PATH=$PATH:/home/lkook/gurobi1202/linux64/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/lkook/gurobi1202/linux64/lib

ms=-1
deg=2
admgadd=3
alpha=0.001
uot=0
simname=simulation
wtype=const
nsim=50
reg=lrm
wt=300
for d in 9; do
  for mode in admg; do
    for n in 400 10000; do
      echo "Running with d=$d, mode=$mode, ms=$ms, deg=$deg, alpha=$alpha, wtype=$wtype, n=$n, reg=$reg"
      Rscript --vanilla ./inst/code/run-simulation.R $mode $d $ms $deg $n $nsim $alpha $uot $simname $wtype $reg $wt $admgadd
    done
  done
done
