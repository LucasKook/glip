#!/usr/bin/env bash
# Example 1 of the paper (Sec. 1.1, Fig. 1): small, runs on a laptop.
# Usage (repository root): bash inst/scripts/intro-example.sh [outdir]
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
exec Rscript inst/code/intro-example.R "${1:-results/paper/intro-example}"
