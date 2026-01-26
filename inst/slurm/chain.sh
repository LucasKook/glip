#!/bin/bash
ms=-1
deg=4
alpha=0.001
uot=0
simname=full
wtype=const
nsim=30
reg=lrm
wt=300
for d in {3..10}; do
	for n in 1000 100000; do
		echo "Running with d=$d, ms=$ms, deg=$deg, alpha=$alpha, wtype=$wtype, n=$n"
		Rscript --vanilla ./inst/code/run-chain-graph-simulation.R $d $ms $deg $n $nsim $alpha $uot $simname $wtype $reg $wt
	done
done
