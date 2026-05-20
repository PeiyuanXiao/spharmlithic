# ==============================================================================
# zzz.R
# Package onLoad hook: lazily imports the bundled Python module.
# ==============================================================================

#' @noRd
sh_py <- NULL

.onLoad <- function(libname, pkgname) {
  # Make the bundled python/ directory importable.
  py_path <- system.file("python", package = pkgname)
  
  # Delay-load: do NOT initialize Python here. The actual import happens
  # the first time sh_py$<something> is accessed by user-facing code, at
  # which point reticulate will resolve the active Python environment.
  #
  # Wrapped in tryCatch so that .onLoad() always completes successfully.
  # If Python is already initialized (e.g. after a devtools::load_all()
  # reload in the same session), reticulate skips delay_load and attempts
  # an immediate import; that may fail if the path isn't yet resolvable.
  # In that case sh_py stays NULL and an informative message is shown;
  # the error surfaces only when a Python-backed function is actually called.
  sh_py <<- tryCatch(
    reticulate::import_from_path(
      "spharmlithic_py",
      path = py_path,
      delay_load = list(
        environment = "r-spharmlithic",
        on_error = function(e) {
          packageStartupMessage(
            "spharmlithic: Python backend not initialized.\n",
            "  Run install_spharmlithic_python() to set it up,\n",
            "  or use_spharmlithic_python('my_env') to point at an existing env.\n",
            "  Original error: ", conditionMessage(e)
          )
        }
      )
    ),
    error = function(e) {
      packageStartupMessage(
        "spharmlithic: could not load Python backend (sh_py = NULL).\n",
        "  If you are reloading within an active session, this is expected.\n",
        "  Call use_spharmlithic_python() before using spharm_* functions.\n",
        "  Original error: ", conditionMessage(e)
      )
      NULL
    }
  )
  
  invisible()
}

# Re-export the magrittr pipe so users don't need library(magrittr).
#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`