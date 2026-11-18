#!/usr/bin/env bash
# Experiment "asp" of the paper (see inst/code/README.md); wrapper of reproduce.sh.
# Usage (repository root): bash inst/scripts/asp.sh [out]   (default results/paper)
exec bash "$(dirname "$0")/reproduce.sh" "${1:-results/paper}" asp
