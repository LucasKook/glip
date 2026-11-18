#!/usr/bin/env bash
# Experiment "benchmark" of the paper (see inst/code/README.md); wrapper of reproduce.sh.
# Usage (repository root): bash inst/scripts/benchmark.sh [out]   (default results/paper)
exec bash "$(dirname "$0")/reproduce.sh" "${1:-results/paper}" benchmark
