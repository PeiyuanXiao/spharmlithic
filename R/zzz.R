# ==============================================================================
# zzz.R
# Package onLoad hook: lazily imports the bundled Python module.
# ==============================================================================

# Module-level reference to the Python submodule.
# Populated by .onLoad(); kept NULL until reticulate finishes initialization.
sh_py <- NULL

.onLoad <- function(libname, pkgname) {
  # Make the bundled python/ directory importable.
  py_path <- system.file("python", package = pkgname)

  # Delay-load: do NOT initialize Python here. The actual import happens the
  # first time sh_py$<something> is accessed by user-facing code, at which
  # point reticulate will resolve the active Python environment.
  sh_py <<- reticulate::import_from_path(
    "lithicscarpattern_py",
    path = py_path,
    delay_load = list(
      environment = "r-lithicscarpattern",
      on_error = function(e) {
        packageStartupMessage(
          "lithicscarpattern: Python backend not initialized.\n",
          "  Run install_lithicscarpattern_python() to set it up,\n",
          "  or use_lithicscarpattern_python('my_env') to point at an existing env.\n",
          "  Original error: ", conditionMessage(e)
        )
      }
    )
  )

  invisible()
}

# Re-export the magrittr pipe so users don't need library(magrittr).
#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`
