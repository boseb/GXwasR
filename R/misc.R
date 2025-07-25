.onAttach <-
    function(lib, pkg, ...) {
        pkgDescription <- utils::packageDescription(pkg)
        pkgVersion <- pkgDescription$Version
        pkgDate <- pkgDescription$Date
        pkgName <- pkgDescription$Package
        pkgTitle <- pkgDescription$Title
        pkgAuthor <- pkgDescription$`Authors@R`
        pkgMaintainer <- pkgDescription$Maintainer
        packageStartupMessage(paste("\n", pkgName, ": ", pkgTitle, sep = ""))
        packageStartupMessage(paste("Version ", pkgVersion, " (", pkgDate, ") installed", sep = ""))
        packageStartupMessage(paste("Author: ", pkgAuthor, sep = ""))
        packageStartupMessage(paste("Maintainer: ", pkgMaintainer, "\n", sep = ""))
        packageStartupMessage("Tutorial: https://github.com")
        packageStartupMessage('Use citation("GXwasR") to know how to cite this work.\n')
    }
