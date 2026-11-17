#!/usr/bin/env bash
# Server (small job, /work/glip mounted): R, Gurobi 13.0.3
# for R next to the image's Gurobi, licence check, R packages, lcd, Python
# environments, clingo check, package tests, task plan, min-length cache
# prebuild and a smoke run. Idempotent: finished steps are skipped.
# Usage: bash /work/glip/setup.sh <code_dir_name>
set -uo pipefail
base=/work/glip
code="$base/$1"
step() { echo "== $1 $(date)"; }
fail() { echo "FAILED: $1"; exit 1; }
[ -d "$code" ] || fail "code $code missing"
[ -x "$base/install-r.sh" ] || fail "$base/install-r.sh missing"
[ -f "$base/secrets/gurobi.lic" ] || fail "$base/secrets/gurobi.lic missing"

step "R"
"$base/install-r.sh" > /tmp/install-r.log 2>&1 || { tail /tmp/install-r.log; fail "install-r"; }
nproc_q=$(awk '{ if ($1 != "max") printf "%d", $1 / $2 }' /sys/fs/cgroup/cpu.max 2>/dev/null)
echo "vCPUs (cgroup): ${nproc_q:-unknown}; memory: $(free -g | awk '/Mem/ {print $2}') GB"

step "Gurobi 13.0.3"
if [ ! -d "$base/gurobi1303/linux64" ]; then
  mkdir -p "$base/gurobi-dl" && cd "$base/gurobi-dl"
  curl -sSfLO https://packages.gurobi.com/13.0/gurobi13.0.3_linux64.tar.gz || fail "download"
  tar -xzf gurobi13.0.3_linux64.tar.gz && mv gurobi1303 "$base/gurobi1303" || fail "unpack"
fi
cat > "$base/env.sh" <<ENV
# GLIP rerun environment (written by setup.sh)
export GUROBI_HOME=$base/gurobi1303/linux64
export PATH=$base/gurobi1303/linux64/bin:\$PATH
export LD_LIBRARY_PATH=$base/gurobi1303/linux64/lib
export GRB_LICENSE_FILE=$base/secrets/gurobi.lic
export R_LIBS=$base/R-lib
export GLIP_PYTHON=$base/venv-glipcut/bin/python
export GLIP_COMP_VENV=$base/venv-comp
export RETICULATE_PYTHON=$base/venv-comp/bin/python
export GLIP_CACHE=$base/cache
export GLIP_ASP_HOME=$base/asp-home
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ENV
. "$base/env.sh"
mkdir -p "$R_LIBS"
if ! Rscript -e 'stopifnot(packageVersion("gurobi") == "13.0.3")' > /dev/null 2>&1; then
  R CMD INSTALL -l "$R_LIBS" "$GUROBI_HOME"/R/gurobi_*.tar.gz > /tmp/gurobi-r.log 2>&1 || { tail /tmp/gurobi-r.log; fail "gurobi R"; }
fi

step "R packages"
Rscript -e '
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"),
        HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))
d <- read.dcf(file.path(commandArgs(TRUE)[1], "DESCRIPTION"))
dep <- function(f) if (f %in% colnames(d)) trimws(sub("\\(.*", "", strsplit(d[, f], ",")[[1]])) else character()
need <- unique(c(dep("Imports"), dep("Depends"), "pcalg", "bnlearn", "dagitty", "igraph", "mvtnorm",
  "reticulate", "deal", "combinat", "hash", "gRbase", "ggm", "pkgload", "testthat", "jsonlite", "slam"))
need <- setdiff(need, c("R", "gurobi", "graph", "RBGL", "Rgraphviz"))
miss <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
cat("missing:", miss, "\n")
bioc <- c("graph", "RBGL")[!vapply(c("graph", "RBGL"), requireNamespace, TRUE, quietly = TRUE)]
if (length(bioc)) { if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = Sys.getenv("R_LIBS")); BiocManager::install(bioc, lib = Sys.getenv("R_LIBS"), update = FALSE, ask = FALSE, Ncpus = parallel::detectCores()) }
if (length(miss)) install.packages(miss, lib = Sys.getenv("R_LIBS"), Ncpus = parallel::detectCores())
still <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
if (length(still)) stop("still missing: ", paste(still, collapse = " "))
' "$code" > /tmp/r-pkgs.log 2>&1 || { tail -20 /tmp/r-pkgs.log; fail "R packages"; }
grep "missing:" /tmp/r-pkgs.log
if ! Rscript -e 'library(lcd)' > /dev/null 2>&1; then
  bash "$code/inst/code/rerun/install-lcd.sh" "$R_LIBS" /tmp/lcd > /tmp/lcd.log 2>&1 || { tail /tmp/lcd.log; fail "lcd"; }
fi

step "licence check (one tiny solve, stop if it fails)"
Rscript -e 'cat("R gurobi", as.character(packageVersion("gurobi")), "\n"); m <- list(A = matrix(1, 1, 1), obj = 1, rhs = 1, sense = ">=", vtype = "B"); r <- gurobi::gurobi(m, list(OutputFlag = 0, Threads = 1)); stopifnot(r$status == "OPTIMAL"); cat("licence OK\n")' || fail "Gurobi 13 licence"

step "Python environments"
if [ ! -x "$GLIP_PYTHON" ]; then
  python3 -m venv "$base/venv-glipcut" && "$GLIP_PYTHON" -m pip install -q --upgrade pip &&
    "$GLIP_PYTHON" -m pip install -q -r "$code/inst/python/requirements.txt" || fail "venv-glipcut"
fi
if [ ! -x "$GLIP_COMP_VENV/bin/python" ]; then
  python3 -m venv "$GLIP_COMP_VENV" && "$GLIP_COMP_VENV/bin/pip" install -q --upgrade pip &&
    "$GLIP_COMP_VENV/bin/pip" install -q torch --index-url https://download.pytorch.org/whl/cpu &&
    "$GLIP_COMP_VENV/bin/pip" install -q dagma==1.1.1 CausalDisco==0.3.0 || fail "venv-comp"
fi
"$GLIP_PYTHON" -c "import gurobipy, networkx; print('gurobipy', gurobipy.gurobi.version(), 'networkx', networkx.__version__)"
"$GLIP_COMP_VENV/bin/python" -c "import dagma, CausalDisco, numpy, torch; print('numpy', numpy.__version__, 'torch', torch.__version__)"

step "ASP code and clingo"
if [ ! -d "$GLIP_ASP_HOME/hyttinen2014uai_ver6" ]; then
  mkdir -p "$GLIP_ASP_HOME" && cd "$GLIP_ASP_HOME" && unzip -q "$code/inst/asp/hyttinen2014uai_ver6.zip" || fail "unzip ASP"
fi
chmod +x "$GLIP_ASP_HOME/hyttinen2014uai_ver6/pkg/ASP/clingo430"
"$GLIP_ASP_HOME/hyttinen2014uai_ver6/pkg/ASP/clingo430" --version | head -1 || fail "clingo"

cd "$code"
step "package tests (one Gurobi thread)"
# as inst/code/test-one-thread.R, without devtools: one Gurobi thread
Rscript -e 'ns <- asNamespace("gurobi"); f <- get("gurobi", ns); unlockBinding("gurobi", ns); assign("gurobi", function(model, params = NULL) { if (is.null(params)) params <- list(); params$Threads <- 1L; f(model, params) }, envir = ns); testthat::test_local(".", stop_on_failure = TRUE)' > /tmp/tests.log 2>&1
tr '\r' '\n' < /tmp/tests.log | grep -E "^\[ FAIL" | tail -1
grep -q "FAIL 0 " <(tr '\r' '\n' < /tmp/tests.log) || { tr '\r' '\n' < /tmp/tests.log | tail -30; fail "package tests"; }

out="$base/runs/rerun-2026-11"
step "plan $out"
if [ ! -f "$out/tasks.csv" ]; then
  Rscript inst/code/rerun/tasks.R plan "$out" 2>&1 | tail -8 || fail "plan"
fi

step "min-length cache prebuild (sequential)"
Rscript inst/code/rerun/tasks.R prebuild "$out" > "$base/prebuild.log" 2>&1 || { tail "$base/prebuild.log"; fail "prebuild"; }
tail -3 "$base/prebuild.log"

step "smoke (one small task per experiment, 10 s)"
smoke="$base/runs/rerun-smoke"
rm -rf "$smoke"
Rscript inst/code/rerun/tasks.R smoke "$smoke" 10 2>&1 | grep -E "start|done|ERROR|exp|asp|bench|chain|simul"
find "$smoke" -name error.txt -exec sh -c 'echo "$1"; cat "$1"' _ {} \;
step "done"
