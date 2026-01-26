#!/bin/bash

### Shared paramters
deg=2
admgadd=3
alpha=0.001
uot=0
wtype=const
nsim=30
reg=lrm
wt=300

### Full
simname=full
for d in {3..9}; do
	for mode in dag admg; do
		for n in 400 10000; do
			echo "Running with d=$d, mode=$mode, ms=-1, deg=$deg, alpha=$alpha, wtype=$wtype, n=$n, reg=$reg"
			Rscript --vanilla ./inst/code/run-simulation.R $mode $d $ms $deg $n $nsim $alpha $uot $simname $wtype $reg $wt $admgadd
		done
	done
done

### Weak
simname=weak
for d in {3..12}; do
	for ms in {1..2}; do
		for mode in dag admg; do
			for n in 400 10000; do
				echo "Running with d=$d, mode=$mode, ms=$ms, deg=$deg, alpha=$alpha, wtype=$wtype, n=$n, reg=$reg"
				Rscript --vanilla ./inst/code/run-simulation.R $mode $d $ms $deg $n $nsim $alpha $uot $simname $wtype $reg $wt $admgadd
			done
		done
	done
done
