# ==============================================================================
# zzz.R
# Package hooks:
#   .onLoad   — silent setup: lazily imports the bundled Python module.
#   .onAttach — user-facing startup message if the backend is unavailable.
# ==============================================================================

#' @noRd
sh_py <- NULL

.onLoad <- function(libname, pkgname) {
  py_path <- system.file("python", package = pkgname)
  
  # Attempt a delay-loaded import of the bundled Python module.
  # Wrapped in tryCatch so that .onLoad() always completes successfully:
  # if Python is already initialised (e.g. after a devtools::load_all()
  # reload in the same session), reticulate skips delay_load and attempts
  # an immediate import; that may fail if the path is not yet resolvable.
  # In that case sh_py stays NULL; the error surfaces only when a
  # Python-backed function is actually called.
  sh_py <<- tryCatch(
    reticulate::import_from_path(
      "spharmlithic_py",
      path = py_path,
      delay_load = list(
        environment = "r-spharmlithic",
        on_error    = function(e) NULL
      )
    ),
    error = function(e) NULL
  )
  
  invisible()
}

# Guard called by every Python-backed function. Turns the cryptic
# "attempt to apply non-function" (raised when sh_py is NULL because the
# backend was never set up) into an actionable message.
#' @noRd
.ensure_backend <- function() {
  if (is.null(sh_py)) {
    stop(
      "spharmlithic: Python backend not initialised.\n",
      "  Run install_spharmlithic_python() to set it up, or\n",
      "  use_spharmlithic_python('r-spharmlithic') to point at an existing env.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.onAttach <- function(libname, pkgname) {
  if (is.null(sh_py)) {
    packageStartupMessage(
      "spharmlithic: Python backend not initialised.\n",
      "  Run install_spharmlithic_python() to set it up,\n",
      "  or use_spharmlithic_python('my_env') to point at an existing env.\n",
      "  If you are reloading within an active session, call ",
      "use_spharmlithic_python() to restore the backend."
    )
  }
}

# Re-export the magrittr pipe so users don't need library(magrittr).
#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`