## Two-line shim for the cluster: lets nb-r1 reuse the ICDT code already deployed in ~/multiDT/TCV-ICDT/common/ICDT
## (same run_icdt(), functions.R and bgfm.cpp; that copy expects ROOT instead of the ~/vsdt mirror).
## Copy this file to <parent of simu-r1 or nb>/ICDT/icdt_wrapper.R, i.e. the path run_icdt_nb.R resolves as ../../ICDT/icdt_wrapper.R.
ROOT <- path.expand("~/multiDT/TCV-ICDT")
source(file.path(ROOT, "common/ICDT/icdt_wrapper.R"))
