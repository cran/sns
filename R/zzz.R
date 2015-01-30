.onAttach <- function(libname, pkgname) {
  RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  packageStartupMessage(paste0("Package: ", pkgname, ", Version: ", RFver))
  packageStartupMessage("Stochastic Newton Sampler for log-concave probability densities.")
  #packageStartupMessage("Metropolis-Hastings sampling of log-concave probability densities using the Stochastic Newton Sampler.")
  packageStartupMessage("Scientific Computing Group, Sentrana Inc. &")
  packageStartupMessage("Imperial College London")
}
