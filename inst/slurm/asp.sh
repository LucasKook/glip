#!/bin/bash

### Shared parameters
nsim=20
uot=0
walltime=600
wtype=const
alpha=0.001
n=10000

### Full
simname=full
ms=-1
ds={3..9}

for d in ds; do
	for mode in dag admg; do
		echo "Running with d=$d, mode=$mode, nsim=$nsim, n=$n"
		Rscript --vanilla ./inst/code/run-asp-comparison.R $mode $d $ms $n $nsim $uot $simname $walltime $wtype $alpha
	done
done

### Weak
simname=weak
ms=1
ds={3..14}

for d in ds; do
	for mode in dag admg; do
		echo "Running with d=$d, mode=$mode, nsim=$nsim, n=$n"
		Rscript --vanilla ./inst/code/run-asp-comparison.R $mode $d $ms $n $nsim $uot $simname $walltime $wtype $alpha
	done
done
