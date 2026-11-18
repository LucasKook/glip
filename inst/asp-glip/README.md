# ASP timing comparison: modified pipeline

The timing comparison against ASP uses the code of Hyttinen, Eberhardt and
Järvisalo (UAI 2014), contained in `inst/asp/hyttinen2014uai_ver6.zip`
(R code, ASP encodings and the clingo binary `pkg/ASP/clingo430` for
Linux, x86-64, which reports version 4.5.4).

`pipeline.R` replaces `pkg/R/pipeline.R` of that code: after clingo has
run, it calls `GLIP_FIT()`, which is defined in `inst/code/exp-asp.R` and
fits GLIP on the same data, tests, PC/FCI warm start and edge hints. The
rest of the file is unchanged; the oracle branch is not supported.

`inst/code/exp-asp.R` unzips the code once (to `GLIP_ASP_HOME`, by default
a temporary directory), copies it to a separate directory per task (the
code writes to `../tmp`) and adds this `pipeline.R`. Set `GLIP_CLINGO` to
use another clingo executable, e.g. on systems where the included binary
does not run.
