## File Name: mnlfa_postproc_trait_parameter_table.R
## File Version: 0.031


mnlfa_postproc_trait_parameter_table <- function(x, type, par0, est0=0)
{
    v1 <- x
    if ( length(v1) > 0 ){
        trait <- data.frame(type=type, par=names(v1), est=v1, estim=1)
    } else {
        trait <- data.frame(type=type, par=par0, est=est0, estim=0)
    }
    return(trait)
}
