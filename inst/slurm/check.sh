#!/bin/bash
for mode in dag admg; do
	for d in {3..10}; do
		for ms in -1; do # for ((ms = 1; ms <= 1; ms++)); do
			echo "Running with d=$d, ms=$ms, mode=$mode"
			Rscript --vanilla ./inst/code/check-random-graphs.R $d 0.5 $mode 1 $ms 1
		done
	done
done
