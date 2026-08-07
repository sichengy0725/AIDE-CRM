# Source this file from any working directory.  In nested source() calls,
# locate this bootstrap file rather than the caller's file.
`%||%` <- function(x, y) if (is.null(x)) y else x
.tite_aide_source_files <- vapply(sys.frames(), function(frame) {
  value <- frame$ofile %||% ""
  if (length(value)) as.character(value)[1L] else ""
}, character(1))
.tite_aide_this_file <- tail(.tite_aide_source_files[basename(.tite_aide_source_files) == "TITE-AIDE.R"], 1L)
if (!length(.tite_aide_this_file) || !nzchar(.tite_aide_this_file)) .tite_aide_this_file <- "TITE-AIDE.R"
.tite_aide_root <- normalizePath(dirname(.tite_aide_this_file), winslash = "/", mustWork = TRUE)
options(tite_aide_root = .tite_aide_root)
for (.f in c("aide_phase12_config.R", "aide_phase12_state.R", "aide_phase12_endpoints.R", "aide_phase12_models.R", "aide_phase12_design.R", "aide_phase12_events.R", "aide_phase12_final.R", "aide_phase12_oc.R", "aide_phase12_oc_setup.R")) source(file.path(.tite_aide_root, "R", .f))
