#!/bin/bash
# Install lcd 0.7-3 (archived on CRAN; rnorm.cg() and get.normal.dist() for
# the chain graph simulation). R >= 4.5 removed Calloc/Realloc/Free, so the
# C sources are patched to R_Calloc/R_Realloc/R_Free first, and the LAPACK
# calls get the hidden string length arguments (FCONE, USE_FC_LEN_T).
# Usage: install-lcd.sh <library> [workdir]
set -euo pipefail
lib="$1"
work="${2:-$(mktemp -d)}"
mkdir -p "$lib" "$work"
cd "$work"
curl -sSfLO https://cran.r-project.org/src/contrib/Archive/lcd/lcd_0.7-3.tar.gz
tar -xzf lcd_0.7-3.tar.gz
perl -pi -e 's/\b(Calloc|Realloc|Free)\(/R_$1(/g' lcd/src/*.c
perl -0pi -e 's/^/#define USE_FC_LEN_T\n/; s/(#include <R_ext\/RS.h>\n)/$1#ifndef FCONE\n#define FCONE\n#endif\n/' lcd/src/citest.c
perl -pi -e 's/z, &nrx\);/z, &nrx FCONE FCONE);/; s/\(double\*\) NULL\);/(double*) NULL FCONE);/; s/(dgecon\)\(.*&info)\);/$1 FCONE);/' lcd/src/citest.c
R CMD INSTALL -l "$lib" lcd
