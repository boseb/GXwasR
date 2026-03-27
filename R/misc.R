## Package Startup Message
.onAttach <- function(libname, pkgname) {
  if (interactive()) {
    packageStartupMessage(
      sprintf("To cite %s, run: citation('%s')", pkgname, pkgname)
    )
  }
}