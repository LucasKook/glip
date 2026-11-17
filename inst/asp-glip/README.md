# ASP timing comparison: modified pipeline

`pipeline.R` replaces `pkg/R/pipeline.R` of the code package of Hyttinen,
Eberhardt and Järvisalo (UAI 2014), which is tracked as
`inst/asp/hyttinen2014uai_ver6.zip` (R code, ASP encodings and the Linux
clingo 4.3.0 binary `pkg/ASP/clingo430` used in the paper).

Changes against the zip: after clingo has run, the pipeline calls
`GLIP_FIT()`, which is defined by the driver `inst/code/rerun/exp-asp.R` and
fits GLIP with every requested encoding (min-length and cuts) on the same
data, tests, PC/FCI warm start and edge hints. It returns one row per
encoding with the clingo time and objective. The oracle branch is not
supported. The rest of the file is unchanged.

The driver unzips the package once (`GLIP_ASP_HOME`), copies `pkg/R` and
`pkg/ASP` into a private directory per task (the code writes to `../tmp`)
and overlays this `pipeline.R`. `GLIP_CLINGO` can point to another clingo
binary (e.g. for a local test on macOS).
