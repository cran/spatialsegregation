
.First.lib <- function(lib, pkg) {
    library.dynam("spatialsegregation", pkg, lib)
    v <- read.dcf(file=system.file("DESCRIPTION", package="spatialsegregation"),
                  fields="Version")
    cat(paste("\nLoaded Spatial Segregation measures (\"spatialsegregation\")", v, "\n"))
}
