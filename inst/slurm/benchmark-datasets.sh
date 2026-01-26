#!/bin/bash
alpha=0.001
uot=0
wtype=const
wt=300
reg=lrm
nsim=30
for mode in admg; do
	for ds in alarm asia child hepar2 sachs; do
		for dmax in 6 8; do
			for ms in 1 2 $((dmax - 2)); do
				echo "Running with ms=$ms, alpha=$alpha, wtype=$wtype, dataset=$ds, dmax=$dmax, nsim=$nsim, reg=$reg"
				Rscript --vanilla ./inst/code/run-benchmark-datasets.R $mode $ds $ms $alpha $uot $wtype $wt $dmax $reg $nsim
			done
		done
	done
done
