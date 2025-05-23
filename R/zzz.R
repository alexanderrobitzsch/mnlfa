## File Name: zzz.R
## File Version: 0.043


#  zzz.R
#
# This function is simply copied from the mice package.


# on attach sirt
.onAttach <- function(libname,pkgname)
{
    d <- utils::packageDescription('mnlfa')
    d1 <- d$Version
    packageStartupMessage( paste('- ', d$Package,' ', d1,' (',d$Date,')',sep='') )
}


version <- function(pkg='mnlfa')
{
    lib <- dirname( system.file(package=pkg))
    d <- utils::packageDescription(pkg)
    return( paste(d$Package,d$Version,d$Date,lib))
}

# .First.lib <- function(lib, pkg){
#          library.dynam('sirt', package=pkg, lib.loc=lib)
#          return(invisible(0))
#        }

xx <- function(f1, f2)
{
    v1 <- paste0( rep(' ',f1), collapse='' )
    v2 <- paste0( rep(' ',f2), collapse='' )
    res <- paste0( v1, '=', v2)
    return(res)
}
