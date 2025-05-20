## File Name: mnlfa_postproc_trait.R
## File Version: 0.026

mnlfa_postproc_trait <- function(parm_trait)
{
    trait <- mnlfa_postproc_trait_parameter_table(x=parm_trait$mu, type='mu',
                        par0='mu0', est0=0)
    dfr1 <- mnlfa_postproc_trait_parameter_table(x=parm_trait$sigma, type='sigma',
                        par0='sigma0', est0=0)
    trait <- rbind( trait, dfr1)
    rownames(trait) <- NULL
    return(trait)
}
