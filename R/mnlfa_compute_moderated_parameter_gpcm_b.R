## File Name: mnlfa_compute_moderated_parameter_gpcm_b.R
## File Version: 0.05


mnlfa_compute_moderated_parameter_gpcm_b <- function(Xdes, parm, value=0)
{
    K <- parm$K
    N <- nrow(Xdes)
    NX <- ncol(Xdes)
    b <- sirt::sirt_matrix2(parm$b[1L:K], nrow=N )
    if (NX>1){
        parm_b <- (parm$b)[K + 1L:(NX-1) ]
        Xdes_b <- Xdes[,-1,drop=FALSE]
        b0 <- mnlfa_compute_moderated_parameter(Xdes=Xdes_b, parm=parm_b, value=value)
        KM <- sirt::sirt_matrix2(1L:K, nrow=N)
        b0 <- KM * matrix(b0, nrow=N, ncol=K)
        b <- b + b0
    }
    return(b)
}
